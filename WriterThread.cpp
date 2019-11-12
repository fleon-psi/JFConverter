#include "jungfrau.h"

// Writer thread
void *WriterThread(void *in_threadarg) {
	int ret;

	WriterThreadArg *in = (WriterThreadArg *)in_threadarg;
	if ((mode == UINT16_MODE) && (px.summation == 1)) {
		buffer[in->threadid] = malloc (px.nimages_per_file*YPIXEL*XPIXEL*sizeof(uint16_t));
		MALLOC_ERROR(buffer[in->threadid]);
	} else {
		buffer[in->threadid] = malloc (px.nimages_per_file*YPIXEL*XPIXEL*sizeof(int32_t));
		MALLOC_ERROR(buffer[in->threadid]);
	}

	char *compression_buffer;
	if (px.use_direct_chunk_writer) {
		compression_buffer = (char *) malloc(wrt->GetOutputBufferMaxSize());
	} else compression_buffer = NULL;

	ret = pthread_mutex_lock(&tasks_assigned_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	while (tasks_assigned < number_of_tasks) {

		int curr_task = tasks_assigned;
		tasks_assigned += 1;

		ret = pthread_mutex_unlock(&tasks_assigned_semaphore);
		PTHREAD_ERROR(ret,pthread_mutex_unlock);

		ret = pthread_mutex_lock(&task_ready_semaphore[curr_task]);
		PTHREAD_ERROR(ret,pthread_mutex_lock);

		task_ready[curr_task] = in->threadid;

		ret = pthread_cond_broadcast(&task_ready_cond[curr_task]);
		PTHREAD_ERROR(ret,pthread_cond_broadcast);

		ret = pthread_mutex_unlock(&task_ready_semaphore[curr_task]);
		PTHREAD_ERROR(ret,pthread_mutex_unlock);

		// Converters are cleared to do their job!

		ret = pthread_mutex_lock(&converter_done_semaphore[curr_task]);
		PTHREAD_ERROR(ret,pthread_mutex_lock);

		while (converter_done[curr_task] != NUM_WORKERS) {
			ret = pthread_cond_wait(&converter_done_cond[curr_task], &converter_done_semaphore[curr_task]);
			PTHREAD_ERROR(ret,pthread_cond_wait);
		}

		ret = pthread_mutex_unlock(&converter_done_semaphore[curr_task]);
		PTHREAD_ERROR(ret,pthread_mutex_unlock);

		wrt->SaveData((char *)buffer[in->threadid], curr_task, compression_buffer);

		ret = pthread_mutex_lock(&tasks_assigned_semaphore);
		PTHREAD_ERROR(ret,pthread_mutex_lock);

	}

	ret = pthread_mutex_unlock(&tasks_assigned_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

	if (px.use_direct_chunk_writer) free(compression_buffer);

	pthread_exit(0);
};
