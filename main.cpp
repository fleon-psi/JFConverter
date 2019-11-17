#include <iostream>
#include <fstream>
#include <cstdlib>
#include <string>
#include <sched.h>
#include <unistd.h>

#include "jungfrau.h"

#if defined(AC922)
#define NCPU 160
#elif defined(HPEDL580)
#define NCPU 48
#endif

int main(int argc, char ** argv) {
    std::cout << "Jungfrau 4M image converter" << std::endl;
    std::cout << "Paul Scherrer Institute" << std::endl;
    std::cout << "Villigen, Switzerland" << std::endl << std::endl;
    std::cout << "Authors: S. Redford, F. Leonarski, A. Mozzanica" << std::endl;
    std::cout << "====================================================" << std::endl;
    int ret;

    parseMetadata(px, argc, argv);

    wrt = new JungfrauHDF5Writer(px);
    images_to_do = px.intended_number_of_images;

    number_of_tasks = images_to_do / px.nimages_per_file;

    if (images_to_do <= 0) {
        std::cerr << "Negative number of images (" << images_to_do << ")." << std::endl;
        std::cerr <<  "Most likely since trigger event happened after data collection finished." << std::endl;
        exit(EXIT_FAILURE);
    }

    if (images_to_do % px.nimages_per_file != 0) number_of_tasks += 1;

    std::cout << "Frames to calculate: " << images_to_do << " (Ignored: " << px.frames_delay << ")" << std::endl;
 
    pixel_mask = (uint32_t *) calloc(YPIXEL*XPIXEL,sizeof(uint32_t));

    for (int i = 0; i < YPIXEL*XPIXEL; i++) {
        int col = i % XPIXEL;
        int row = i / XPIXEL;
        if ((col > 1029) && (col < 1038)) pixel_mask[i] = 1;
        if ((row > 0*(514+36)+514-1) && (row < 1*(514+36))) pixel_mask[i] = 1;
        if ((row > 1*(514+36)+514-1) && (row < 2*(514+36))) pixel_mask[i] = 1;
        if ((row > 2*(514+36)+514-1) && (row < 3*(514+36))) pixel_mask[i] = 1;
    }

    pthread_t worker[NUM_WORKERS];
    ConverterThreadArg worker_arg[NUM_WORKERS];
    
    pthread_t writer[px.number_of_writers];
    WriterThreadArg writer_arg[px.number_of_writers];

    task_ready           = (int *)             malloc(number_of_tasks*sizeof(int));
    task_ready_semaphore = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t)*number_of_tasks);
    task_ready_cond      = (pthread_cond_t *)  malloc(sizeof(pthread_cond_t)*number_of_tasks);

    converter_done           = (int *)             calloc(number_of_tasks,sizeof(int));
    converter_done_semaphore = (pthread_mutex_t *) malloc(sizeof(pthread_mutex_t)*number_of_tasks);
    converter_done_cond      = (pthread_cond_t *)  malloc(sizeof(pthread_cond_t)*number_of_tasks);

    for (int i = 0; i < number_of_tasks; i++) {
       pthread_mutex_init(&(task_ready_semaphore[i]), NULL);
       pthread_mutex_init(&(converter_done_semaphore[i]), NULL);
       pthread_cond_init(&(task_ready_cond[i]), NULL);
       pthread_cond_init(&(converter_done_cond[i]), NULL);
       task_ready[i] = -1;
    }
    
    for (int i = 0; i < NMODULES; i++) {
        pixels_in_G1[i] = (int *) calloc(images_to_do,sizeof(int));
        pixels_in_G2[i] = (int *) calloc(images_to_do,sizeof(int));
    }
    
    // Initialize worker threads
    for(int t=0; t<NUM_WORKERS; t++){
        // initialize thread attribute
        pthread_attr_t threadattr;
        ret = pthread_attr_init(&threadattr);
        
        // set cpu
#if defined(HPEDL580) || defined(AC922)

        cpu_set_t *cpusetp;
        cpusetp = CPU_ALLOC(NCPU);
        
        if (cpusetp == NULL) exit(EXIT_FAILURE);
        size_t size = CPU_ALLOC_SIZE(NCPU);
        CPU_ZERO_S(size,cpusetp);
        int num_cpu;
        
        CPU_SET(worker_thread_cpu[t], cpusetp);
        
        pthread_attr_setaffinity_np(&threadattr, size ,cpusetp);
#endif
        worker_arg[t].threadid = t;
        worker_arg[t].module = t%NMODULES;
        worker_arg[t].worker_in_module = t/NMODULES;
        ret = pthread_create(&worker[t],&threadattr, ConverterThread, &worker_arg[t]);
        PTHREAD_ERROR(ret,pthread_create);
        
        // Destroy thread attribute
        ret = pthread_attr_destroy(&threadattr);
#if defined(HPEDL580) || defined(AC922)
        CPU_FREE(cpusetp);
#endif
    }
    
    // Initialize writer threads
    for(int t=0; t< px.number_of_writers; t++){
        writer_arg[t].threadid = t;
        
        // initialize thread attribute
        pthread_attr_t threadattr;
        ret = pthread_attr_init(&threadattr);

#if defined(HPEDL580) || defined(AC922)
        // set cpu
        cpu_set_t *cpusetp;
        cpusetp = CPU_ALLOC(NCPU);
        
        if (cpusetp == NULL) exit(EXIT_FAILURE);
        size_t size = CPU_ALLOC_SIZE(NCPU);
        CPU_ZERO_S(size,cpusetp);
        int num_cpu;
        
        CPU_SET(writer_thread_cpu[t], cpusetp);
        
        pthread_attr_setaffinity_np(&threadattr, size, cpusetp);
#endif
        
        ret = pthread_create(&writer[t],&threadattr, WriterThread, &writer_arg[t]);
        PTHREAD_ERROR(ret,pthread_create);
        
        // Destroy thread attribute
        ret = pthread_attr_destroy(&threadattr);
#if defined(HPEDL580) || defined(AC922)
        CPU_FREE(cpusetp);        
#endif
    }
    

    // Wait till all threads finish
    for(int t=0; t< NUM_WORKERS; t++){
        ret = pthread_join(worker[t], NULL);
        PTHREAD_ERROR(ret,pthread_join);
    }


    for(int t=0; t< px.number_of_writers; t++){
        ret = pthread_join(writer[t], NULL);
        PTHREAD_ERROR(ret,pthread_join);
    } 

    wrt->SavePixelMask(pixel_mask);

            int *gain1 = (int *) calloc(images_to_do,sizeof(int));
            int *gain2 = (int *) calloc(images_to_do,sizeof(int));

            int tmp1 = 0, tmp2 = 0;
            for (int i = 0; i < NMODULES; i++) 
                for (int j = 0; j < images_to_do; j++) {
                    tmp1 += pixels_in_G1[i][j];
                    tmp2 += pixels_in_G2[i][j];

                    gain1[j] += pixels_in_G1[i][j];
                    gain2[j] += pixels_in_G2[i][j];
                }

            std::cout << "pixels G1:" << tmp1 << " G2:" << tmp2 << std::endl;

            int low = 0;
            while ((low < images_to_do) && (gain1[low] < 20)) low++;

            int high = images_to_do - 1;
            while ((high > 0) && (gain1[high] < 20)) high--;

            std::cout << std::endl;
            std::cout << "DATA_RANGE= " << low << " " << high << std::endl;

	    wrt->SaveDetectorSpecificData("pixels_in_G1",gain1,images_to_do);
	    wrt->SaveDetectorSpecificData("pixels_in_G2",gain2,images_to_do);

            for (int i = 0; i < NMODULES; i++) {
   	        wrt->SaveDetectorSpecificData("pixels_in_G1",pixels_in_G1[i],images_to_do, i);
	        wrt->SaveDetectorSpecificData("pixels_in_G2",pixels_in_G2[i],images_to_do, i);
            }
            free(gain1);
            free(gain2);
    

    for (int i = 0; i < NMODULES; i++) {
        wrt->SaveDetectorSpecificData("gain_map_file_mod",gainMapFiles[i],i);
        wrt->SaveDetectorSpecificData("init_mask_file_mod",pixelMaskFiles[i],i);
        wrt->SaveDetectorSpecificData("temp_fpga",px.tempFPGA[i],i);
    }

    for (int i = 0; i < NMODULES; i++) {
        free(pixels_in_G1[i]);
        free(pixels_in_G2[i]);              
    }

    delete wrt;

    free(task_ready);
    free(task_ready_semaphore);
    free(task_ready_cond);

    free(converter_done);
    free(converter_done_semaphore);
    free(converter_done_cond);

    std::cout << "Done." << std::endl;
};

