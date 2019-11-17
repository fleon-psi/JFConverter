#ifndef _JUNGFRAU_H
#define _JUNGFRAU_H

#define PEDESTAL_DRIFT
#define ALLOW_NEGATIVE

#include <string>
#include <algorithm> // std::count
#include <iostream>  // std::cout
#include <stdint.h>
#include <ctime>     // clock
#include <iomanip>   // std::setprecision
#include <fstream>
#include <vector>

#include <time.h>
#include <sys/time.h>
#include <hdf5.h>

#include <pthread.h>

#define PIXEL_SIZE_IN_MM 0.075

#ifdef JF_1M_2kHz

#define NC 1024
#define NR 256

#define MODULESX 1
#define MODULESY 2
#define NMODULES (MODULESX*MODULESY*2)
#define DETECTOR_NAME "JF1M"

#define NEW_RECEIVER

#else

#define NC 1024
#define NR 512

#define MODULESX 2
#define MODULESY 4

#define NMODULES (MODULESX*MODULESY)
#define DETECTOR_NAME "JF4M"

#endif

#define NMODULES_GEOM (MODULESX*MODULESY)

#define NCH (NC*NR)

#define IMAGES_IN_FILE 10000

#define GAPX	 8
#define GAPY     36

#define XPIXEL (1030*MODULESX + GAPX*(MODULESX-1))
#define YPIXEL (514*MODULESY + GAPY*(MODULESY-1))

#define PEDE_G0 1000

#define PEDE_WINDOW 100

#ifdef AC922
#define NUM_WORKERS 32
#define MAX_NUM_WRITERS 48
#else
#define NUM_WORKERS 32
#define MAX_NUM_WRITERS 16
#endif

#define CHANNELS_PER_WORKER (NCH/(NUM_WORKERS/NMODULES))

#define PTHREAD_ERROR(ret,func) if (ret) printf("%s(%d) %s: err = %d\n",__FILE__,__LINE__, #func, ret), exit(ret)
#define HDF5_ERROR(ret,func) if (ret) printf("%s(%d) %s: err = %d\n",__FILE__,__LINE__, #func, ret), exit(ret)
#define MALLOC_ERROR(ret) if (ret == NULL) printf("%s(%d): memory allocation error\n",__FILE__,__LINE__), exit(EXIT_FAILURE)

// for bitshuffle H5 filter
#define BSHUF_H5_COMPRESS_LZ4  2
#define BSHUF_H5_COMPRESS_ZSTD 3

#define BSHUF_H5FILTER    32008
#define LZ4_H5FILTER      32004
#define ZSTD_H5FILTER 	  32015
#define SZ_H5FILTER       32017

#define FRAME_BUF_SIZE 200

typedef enum
{
    INT32_MODE, INT16_MODE, UINT16_MODE, UINT32_MODE, FLOAT_MODE
} ConverterMode;

typedef enum
{
    COMPRESSION_NONE,
	COMPRESSION_LZ4,
	COMPRESSION_BSHUF_LZ4,
	COMPRESSION_GZIP,
	COMPRESSION_SZ,
	COMPRESSION_ZSTD,
	COMPRESSION_BSHUF_ZSTD
} ConverterCompression;

/* Metadata for an experiment, necessary for data processing
 * In future should be imported directly from DAQ handler
 */
struct PXMetadata {
    ConverterMode mode;

    int number_of_writers;

    int nimages_collected; // total number of images collected during the experiment (raw)

    int intended_number_of_images;

    int summation; // summation amount
    int nimages_per_file; // # of summed images to save in a single HDF5 file

    int frames_pedeG0; // Frames collected for G0 pedestal
    int frames_pedeG1; // Frames collected for G1 pedestal
    int frames_pedeG2; // Frames collected for G2 pedestal

    int frames_delay;

    std::string hdf5_prefix; // prefix for HDF5 files
    std::string jf_raw_prefix; // prefix for JF files

    double detector_distance; // in mm
    double beamx,beamy; // in pixel
    double omega_start; // in degree

    double omega_increment; // in degree
    double integration_time; // in us
    double frame_time; // in us
    double photon_energy; // in eV
    double photon_energy_for_normalization; // in eV
    double sensor_thickness; // in um
    double x_pixel_size, y_pixel_size; // in micron

    std::vector<std::string> dataFile;

    bool sumFloat; // use float summation
    double tempFPGA[NMODULES];

    ConverterCompression compressionAlgorithm;
    bool use_direct_chunk_writer;
    size_t bshuf_block_size;
    double absolute_err_bound_sz;
    int gzip_level;

    bool stillImage;

    double transmission;

};

// Representation of a frame in the file.
struct JF_rawheader {
#ifdef JF_1M_2kHz
    uint64_t framenumber;
    uint64_t packetnumber;
    uint64_t bunchid;
    uint64_t thing4;
    uint64_t thing5;
    uint64_t thing6;
    uint64_t thing7;
    uint64_t thing8;
    uint64_t thing9;
    uint64_t thing10;
    uint64_t thing11;
    uint64_t thing12;
    uint64_t thing13;
    uint64_t thing14;
#else
    uint64_t framenumber;
    uint32_t subframenumber;
    uint32_t packetnumber;
    uint64_t bunchid;
    uint64_t timestamp;
    uint16_t moduleid;
    uint16_t xcoord;
    uint16_t ycoord;
    uint16_t zcoord;
    uint32_t debug;
    uint16_t rrnumber;
    uint8_t detectortype;
    uint8_t headerversion;
#endif
};

struct JF_frame {
    uint16_t imagedata[NCH];
};

/*
 * Class implementing pixel mask operations
 */
class JungfrauPixelMask {
uint32_t* mask;

public:
    JungfrauPixelMask(uint32_t* mask);
    void MaskIfPedestalStep(uint16_t *imagedata, float* pedestalG0, int step);
    void MaskIfPedeRMSG0GreaterThan(float* pedestalG0RMSs, int val);
    void MaskIfGainNot(int expected, uint16_t *imagedata);
};


class JungfrauFile {
    std::ifstream filebin; // file descriptor
public:
    JungfrauFile(char *fformat, int fileindex);
    ~JungfrauFile();
    bool ReadHeader(JF_rawheader *header);
    bool IgnoreHeader();
    bool IgnoreFrames(uint32_t nframes);
    bool IgnoreFramesFromBeg(uint32_t nframes);
    bool ReadFrame(uint16_t *frame, uint32_t store, uint32_t ignore_in_front, uint32_t ignore_afterwards);
    bool ReadFrame(uint16_t *frame, uint32_t start_channel);
};

class JungfrauConverter {
    float gainG0[CHANNELS_PER_WORKER] __attribute__((aligned(64)));
    float gainG1[CHANNELS_PER_WORKER] __attribute__((aligned(64)));
    float gainG2[CHANNELS_PER_WORKER] __attribute__((aligned(64)));

    float pedeG0[CHANNELS_PER_WORKER] __attribute__((aligned(64)));
    float pedeG1[CHANNELS_PER_WORKER] __attribute__((aligned(64)));
    float pedeG2[CHANNELS_PER_WORKER] __attribute__((aligned(64)));
    float pedeG0_RMS[CHANNELS_PER_WORKER] __attribute__((aligned(64)));

    float localBuffer[CHANNELS_PER_WORKER] __attribute__((aligned(64)));
    uint32_t channelMask[CHANNELS_PER_WORKER] __attribute__((aligned(64)));

    uint16_t jf_frame[CHANNELS_PER_WORKER] __attribute__((aligned(64)));
    uint16_t jf_frame_buffer[FRAME_BUF_SIZE*CHANNELS_PER_WORKER] __attribute__((aligned(64)));

    int start_pos;
    int start_channel;

    JungfrauFile *thisfile;

    long images_read;
    long frames_written;
    int curr_file;
    int even_module; // even_module == 0: lines   0..255, even_module == 1: lines 256..511

    uint16_t summation;

    void ReloadFile(int summation);

    void CopyLine(uint16_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local);
    void CopyLine(int16_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local);
    void CopyLine(uint32_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local);
    void CopyLine(int32_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local);
    void CopyLine(float *imageBuffer, int start_pos_in_buffer, int start_pos_in_local);

    void CopyLineWithHorizGap(uint16_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap);
    void CopyLineWithHorizGap(int16_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap);
    void CopyLineWithHorizGap(uint32_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap);
    void CopyLineWithHorizGap(int32_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap);
    void CopyLineWithHorizGap(float *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap);

    void LoadGain(std::string gainMap, double energy_in_keV);
    void LoadPixelMask(std::string pixelMaskFile);

    void CalculatePede(float *gain, uint16_t expected, uint32_t size);
    void CalculatePede(float *gain, float *gainRMS, uint16_t expected, uint32_t size);
    void CalculatePede(float *gain, float *gainRMS, float *maxRMS, uint16_t expected, uint32_t size);

    std::string dataFileName;
public:
    JungfrauConverter(std::string data, std::string gainMapFile, std::string pixelMapFile, uint32_t *pixel_mask, double energy_in_keV, int nimg_pedeG0, int nimg_pedeG1,
    		          int nimg_pedeG2, int in_start_pos, uint16_t summation, int in_start_channel = 0, int in_even_module = 0);
    ~JungfrauConverter();

    template<typename T> void Convert(T *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
    template<typename T> void ConvertSum1(T *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
};

/* HDF5 Writer
 */
class JungfrauHDF5Writer {
    std::string prefix;
    int images_per_file; // images in single hdf5 file
    int nimages; // total images to save in hdf5 (after summation)

    int nchunks_used; // number of chunks already saved
    int nimages_saved; // # of images saved
    int nchunks; // Total number of chunks
    hid_t master_file_id; // HDF5 handle for master file

    size_t elem_size;
    hid_t data_type;

    ConverterCompression compressionAlgorithm;
    bool use_direct_chunk_writer;
    size_t bshuf_block_size;
    int gzip_level;
    size_t output_size;
    char *compression_buffer;

    void InitCompression();
    void SetCompressionFilter(hid_t data_dcpl_id);

    // Save metadata according to NeXuS specification
    void SaveMetadata(const PXMetadata& in_metadata);
    // Create "angle" container (phi/omega/chi/...) and fill it with incrementing numbers.
    int SaveAngleContainer(hid_t location, std::string name, double start, double increment, std::string units);

    // Maximal output size
    size_t CompressionOutputMaxSize();
    // Compress for direct chunk writing
    size_t Compress(char *output, char *input);

public:
    // Create JungfrauHDF5writer object
    // Creates a HDF5 master and writes metadata
    JungfrauHDF5Writer(const PXMetadata& in_metadata);

    // Close file
    ~JungfrauHDF5Writer();

    // Writes modules metrology
    void WriteMetrology(const PXMetadata& in_metadata);

    // Saves data under /entry/instrument/detector/detectorSpecific
    // with name as a container name.
    void SaveDetectorSpecificData(std::string name, std::vector<double> &data);
    void SaveDetectorSpecificData(std::string name, std::vector<int> &data);
    void SaveDetectorSpecificData(std::string name, std::string data);
    void SaveDetectorSpecificData(std::string name, int *data, size_t size);

    void SaveDetectorSpecificData(std::string name, std::vector<double> &data, int module);
    void SaveDetectorSpecificData(std::string name, std::vector<int> &data, int module);
    void SaveDetectorSpecificData(std::string name, std::string data, int module);
    void SaveDetectorSpecificData(std::string name, double data, int module);
    void SaveDetectorSpecificData(std::string name, int *data, size_t size, int module);

    // Save pixel mask
    void SavePixelMask(uint32_t *pixelMask);

    // Save images
    void SaveData(char* data, int chunk_number, char *buffer = NULL);

    // Return maximal size of compression buffer
    size_t GetOutputBufferMaxSize();
};

int addStringAttribute(hid_t location, std::string name, std::string val);
int addDoubleAttribute(hid_t location, std::string name, const double *val, int dim);
hid_t createGroup(hid_t master_file_id, std::string group, std::string nxattr);
int saveDouble(hid_t location, std::string name, double val, std::string units = "");
int saveDouble1D(hid_t location, std::string name, const double *val, std::string units = "", int dim = 1);

int saveInt(hid_t location, std::string name, int val, std::string units = "");
int saveInt1D(hid_t location, std::string name, const int *val, std::string units = "", int dim = 1);
int saveUInt2D(hid_t location, std::string name, const uint32_t *val, std::string units, int dim1, int dim2);

int saveString(hid_t location, std::string name, std::string val, std::string units = "");
int saveString1D(hid_t location, std::string name, char *val, std::string units, int dim, int len);

// int createDataChunk(std::string filename, int frame_low, int *data, int nframes);
int createDataChunkLink(std::string filename, hid_t location, std::string name);

int parseMetadata(PXMetadata &px, int argc, char **argv);

struct ConverterThreadArg {
    int threadid;
    int module;
    int worker_in_module;
};

struct WriterThreadArg {
    int threadid;
    bool ready;
    bool finished;
    int chunk;
};

void *WriterThread(void *in_threadarg);
void *ConverterThread(void *in_threadarg);

size_t findFirstFrameWithBeam();
size_t findLastFrameWithBeam();

// constants
extern int worker_thread_cpu[]; // CPU assigned to converter threads
extern int writer_thread_cpu[]; // CPU assigned to writer threads
extern std::string gainMapDir; // directory with gain data

extern std::vector<std::string> gainMapFiles; // filename with gain for each module
extern std::vector<std::string> pixelMaskFiles; // filename with gain for each module

extern int firstChannel[]; // first channel
extern int startPos[];     // start position of module in full JF4M coordinate system

extern double corner[NMODULES_GEOM][2];
extern double fast[NMODULES_GEOM][3];
extern double slow[NMODULES_GEOM][3];

// set up startup and not modified
extern PXMetadata px;
extern size_t number_of_tasks;
extern long images_to_do;

// HDF5 writer
extern JungfrauHDF5Writer *wrt;

// How many writers are ready
extern pthread_mutex_t *task_ready_semaphore;
extern pthread_cond_t *task_ready_cond;
extern int *task_ready;

// Synchronization of c threads
extern pthread_mutex_t *converter_done_semaphore;
extern pthread_cond_t *converter_done_cond;
extern int *converter_done;

extern pthread_mutex_t tasks_assigned_semaphore;
extern int tasks_assigned;

// Synchronization of converter threads
extern pthread_mutex_t semaphore;

// protected with sempahore
extern uint32_t *pixel_mask;
extern int *pixels_in_G1[NMODULES];
extern int *pixels_in_G2[NMODULES];

// protected with writer_sempahore[]
extern void *buffer[MAX_NUM_WRITERS];

// Synchronization of HDF5 usage
extern pthread_mutex_t hdf5_semaphore;

#endif /* _JUNGFRAU_H */
