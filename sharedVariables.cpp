#include "jungfrau.h"

#ifdef AC922
// NUMA: 0 1
int worker_thread_cpu[] =  {
    0,  4,  8, 12,
   16, 20, 24, 28,
   32, 36, 40, 44,
   48, 52, 56, 60,
   80, 84, 88, 92,
   96,100,104,108,
  112,116,120,124,
  128,132,136,140
};

int writer_thread_cpu[] =  {
   64,144, 68,148,
   72,152, 76,154,
    1, 81,  5, 85,
    9, 89, 13, 93,
   17, 97, 21,101,
   25,105, 29,109,
   33,113, 37,117,
   41,121, 45,125,
   49,129, 53,133,
   57,137, 61,141,
   65,145, 69,149,
   73,153, 77,157
};

#elif HPEDL580

// NUMA: 2 0 3 1
int worker_thread_cpu[] =  {
    28, 29, 30, 31,
    32, 33, 34, 35,
     4,  5,  6,  7,
     8,  9, 10, 11,
    40, 41, 42, 43,
    44, 45, 46, 47,
    16, 17, 18, 19,
    20, 21, 22, 23
 };
 
 int writer_thread_cpu[] =  {
    0, 12, 24, 36,
    1, 13, 25, 37,
    2, 14, 26, 38,
    3, 15, 27, 39
 };

#endif

#ifdef JF_1M_2kHz

int firstChannel[] = {0, 128*1024, 64*1024, 192*1024};

int startPos[] = {0, 1030*258, 1030*(514+GAPY), 1030*(514+258+GAPY)};

#ifdef HPEDL580
std::string gainMapDir = "/home/jungfrau/daq/JF1M_DLS_Oct2019/calibration/";
std::vector<std::string> gainMapFiles = {"gainMaps_M264.bin",
    "gainMaps_M264.bin",
    "gainMaps_M253.bin",
    "gainMaps_M253.bin"};
#elif AC922
std::string gainMapDir = "/home/jungfrau/daq/JF1M_DLS_Oct2019/calibration/";
std::vector<std::string> gainMapFiles = {"gainMaps_M264.bin",
    "gainMaps_M264.bin",
    "gainMaps_M253.bin",
    "gainMaps_M253.bin"};
#else
std::string gainMapDir = "/mnt/das-gpfs/home/leonarski_f/p16371/JF1M_2kHz_test/";
std::vector<std::string> gainMapFiles = {"gainMaps_M264_2019-07-29.bin",
    "gainMaps_M264_2019-07-29.bin",
    "gainMaps_M253_2019-07-29.bin",
    "gainMaps_M253_2019-07-29.bin"};
#endif

std::vector<std::string> pixelMaskFiles = {"",
    "",
    "",
    "",
    "",
    "",
    "",
    ""};

double corner[2][2] = {{0.0,0.0},{0.0,514+36}};
double fast[2][3] = {{1.0,0.0,0.0},{1.0,0.0,0.0}};
double slow[2][3] = {{0.0,1.0,0.0},{0.0,1.0,0.0}};
#else

int firstChannel[] = {0, 256*1024, 128*1024, 384*1024};

std::string gainMapDir = "/home/jungfrau/daq/JF4M_PX_Aug2018/calibration/";

std::vector<std::string> gainMapFiles = {"gainMaps_M120.bin",
    "gainMaps_M210.bin",
    "gainMaps_M202.bin",
    "gainMaps_M115.bin",
    "gainMaps_M049.bin",
    "gainMaps_M043.bin",
    "gainMaps_M060.bin",
    "gainMaps_M232.bin"};

std::vector<std::string> pixelMaskFiles = {"trap_mask_v2_0.txt",
    "trap_mask_v2_1.txt",
    "trap_mask_v2_2.txt",
    "trap_mask_v2_3.txt",
    "trap_mask_v2_4.txt",
    "trap_mask_v2_5.txt",
    "trap_mask_v2_6.txt",
    "trap_mask_v2_7.txt"};

int startPos[] = {0, 1038,
    1137400, 1138438,
    2274800, 2275838,
    3412200, 3413238};

double corner[8][2] = {{0.0,    0.0   },{1038,    0.0   },
                       {0.0, 514+36   },{1038, 514+36   },
                       {0.0,(514+36)*2},{1038,(514+36)*2},
                       {0.0,(514+36)*3},{1038,(514+36)*3}};

double fast[8][3] = {{1.0,0.0,0.0},
		     {1.0,0.0,0.0},
		     {1.0,0.0,0.0},
		     {1.0,0.0,0.0},
		     {1.0,0.0,0.0},
		     {1.0,0.0,0.0},
		     {1.0,0.0,0.0},
		     {1.0,0.0,0.0}};

double slow[8][3] = {{0.0,1.0,0.0},
		     {0.0,1.0,0.0},
		     {0.0,1.0,0.0},
		     {0.0,1.0,0.0},
		     {0.0,1.0,0.0},
		     {0.0,1.0,0.0},
		     {0.0,1.0,0.0},
                     {0.0,1.0,0.0}};
#endif


pthread_mutex_t semaphore = PTHREAD_MUTEX_INITIALIZER;
pthread_mutex_t hdf5_semaphore = PTHREAD_MUTEX_INITIALIZER;

long images_done;
long images_to_do;
int threads_done;
int frame_low;

uint32_t *pixel_mask;
void *buffer[MAX_NUM_WRITERS];

int curr_writer;
bool all_done;

PXMetadata px;

int *pixels_in_G1[NMODULES];
int *pixels_in_G2[NMODULES];

JungfrauHDF5Writer *wrt;

ConverterMode mode = SIMPLE_MODE;

bool keep_negative_numbers = false;

pthread_mutex_t tasks_assigned_semaphore = PTHREAD_MUTEX_INITIALIZER;
int tasks_assigned = 0;

// How many writers are ready
pthread_mutex_t *task_ready_semaphore;
pthread_cond_t *task_ready_cond;
int *task_ready;

// How many readers are done
pthread_mutex_t *converter_done_semaphore;
pthread_cond_t *converter_done_cond;
int *converter_done;

size_t number_of_tasks;

