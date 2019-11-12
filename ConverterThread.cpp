 #include "jungfrau.h"

void *ConverterThread(void *in_threadarg)
{
    int ret;
    ConverterThreadArg *arg = (ConverterThreadArg *) in_threadarg;

#ifdef JF_1M_2kHz
    JungfrauConverter *converter = new JungfrauConverter(px.dataFile[arg->module], gainMapDir + gainMapFiles[arg->module], "", pixel_mask, px.photon_energy_for_normalization/1000.0, px.frames_pedeG0, px.frames_pedeG1, px.frames_pedeG2, startPos[arg->module], px.summation, firstChannel[arg->worker_in_module], arg->module % 2);
#else
    JungfrauConverter *converter = new JungfrauConverter(px.dataFile[arg->module], gainMapDir + gainMapFiles[arg->module], gainMapDir + pixelMaskFiles[arg->module], pixel_mask, px.photon_energy_for_normalization/1000.0, px.frames_pedeG0, px.frames_pedeG1, px.frames_pedeG2, startPos[arg->module], px.summation, firstChannel[arg->worker_in_module], 0);
#endif

    for (int task = 0; task < number_of_tasks; task++) {
        // check if task is ready (i.e. buffer is cleared)
        ret = pthread_mutex_lock(&task_ready_semaphore[task]);
        PTHREAD_ERROR(ret,pthread_mutex_lock);

        while (task_ready[task] == -1) {
            ret = pthread_cond_wait(&task_ready_cond[task], &task_ready_semaphore[task]);
            PTHREAD_ERROR(ret,pthread_cond_wait);
        }

        int curr_buffer = task_ready[task];

        ret = pthread_mutex_unlock(&task_ready_semaphore[task]);
        PTHREAD_ERROR(ret,pthread_mutex_unlock);


        int frames = px.nimages_per_file;
        if (task == number_of_tasks - 1)
            frames = images_to_do - (number_of_tasks-1)*px.nimages_per_file;

        // Do work!
        // If each worker is writing to not overlapping memory locations in the buffer,
        // each thread can operate concurrently.
        std::vector<int> local_pixels_in_G1(frames);
        std::vector<int> local_pixels_in_G2(frames);

        switch (mode) {
        case SIMPLE_MODE:
        	if (px.summation == 1)
        		converter->ConvertSum1<int32_t>((int32_t *) buffer[curr_buffer], local_pixels_in_G1, local_pixels_in_G2, frames);
        	else
        		converter->Convert<int32_t>((int32_t *) buffer[curr_buffer], local_pixels_in_G1, local_pixels_in_G2, frames);
        	break;

        case FLOAT_MODE:
        	if (px.summation == 1)
        		converter->ConvertSum1<float>((float *) buffer[curr_buffer], local_pixels_in_G1, local_pixels_in_G2, frames);
        	else
        		converter->Convert<float>((float *) buffer[curr_buffer], local_pixels_in_G1, local_pixels_in_G2, frames);
        	break;

        case UINT16_MODE:
        	if (px.summation == 1)
        		converter->ConvertSum1<uint16_t>((uint16_t *) buffer[curr_buffer], local_pixels_in_G1, local_pixels_in_G2, frames);
        	else
        		converter->Convert<uint16_t>((uint16_t *) buffer[curr_buffer], local_pixels_in_G1, local_pixels_in_G2, frames);
        	break;
        case UINT32_MODE:
        	if (px.summation == 1)
        		converter->ConvertSum1<uint32_t>((uint32_t *) buffer[curr_buffer], local_pixels_in_G1, local_pixels_in_G2, frames);
        	else
        		converter->Convert<uint32_t>((uint32_t *) buffer[curr_buffer], local_pixels_in_G1, local_pixels_in_G2, frames);
        }

        // Acquire lock to finalize writing.
        ret = pthread_mutex_lock(&semaphore);
        PTHREAD_ERROR(ret,pthread_mutex_lock);

        for (int i = 0; i < frames; i++) {
            pixels_in_G1[arg->module][task*px.nimages_per_file+i] += local_pixels_in_G1[i];
            pixels_in_G2[arg->module][task*px.nimages_per_file+i] += local_pixels_in_G2[i];
        }

        ret = pthread_mutex_unlock(&semaphore);
        PTHREAD_ERROR(ret,pthread_mutex_unlock);

        // Mark as done + wake-up HDF5 writer if all tasks were finished
        ret = pthread_mutex_lock(&converter_done_semaphore[task]);
        PTHREAD_ERROR(ret,pthread_mutex_lock);

        converter_done[task]++;
        if (converter_done[task] == NUM_WORKERS) {
            ret = pthread_cond_broadcast(&converter_done_cond[task]);
            PTHREAD_ERROR(ret,pthread_cond_broadcast);
        }

        ret = pthread_mutex_unlock(&converter_done_semaphore[task]);
        PTHREAD_ERROR(ret,pthread_mutex_unlock);
    }

    delete converter;
    pthread_exit(0);
}

