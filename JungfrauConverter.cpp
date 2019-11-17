#include <cmath>
#include <fstream>
#include <cstdlib>
#include "stdint.h"

#include "jungfrau.h"

#define RMS_THRESHOLD 100
#define TIMES_RMS_THRESHOLD_FOR_PEDE 2
#define ROUND_LOC_BUFFER(x) (x)

inline uint16_t jf_round_u16(float x) {
	if (x >= (float) UINT16_MAX) return (uint16_t) (UINT16_MAX-2);
	if (x <= (float) INT16_MIN)   return (uint16_t) (UINT16_MAX-1);
	else if (x > 0.0f) return (uint16_t) (x + 0.5f);
	else return (uint16_t) 0;
}

inline uint16_t jf_round_i16(float x) {
	if (x >= (float) INT16_MAX) return (uint16_t) (INT16_MAX-1);
	if (x <= (float) INT16_MIN) return (uint16_t) (INT16_MIN+1);
	else if (x > 0.0f) return (uint16_t) (x + 0.5f);
	else return (uint16_t) (x-0.5f);
}

inline uint32_t jf_round_u32(float x) {
	if (x >= (float) INT32_MAX) return (uint32_t) (INT32_MAX / 2 - 1);
	if (x <= (float) INT32_MIN) return (uint32_t) (INT32_MAX / 2 + 1);
	else if (x >= 0.0f) return (uint32_t) (x + 0.5f);
	else return (int32_t) 0;
}

inline int32_t jf_round_i32(float x) {
	if (x >= (float) INT32_MAX) return (int32_t) (INT32_MAX / 2 - 1);
	if (x <= (float) INT32_MIN) return (int32_t) (INT32_MIN / 2 + 1);
	else if (x >= 0.0f) return (int32_t) (x + 0.5f);
	else return (int32_t) (x - 0.5f);
}

// Assumption that pedestal always fits into a single file
JungfrauConverter::JungfrauConverter (std::string data, std::string gainMapFile, std::string pixelMaskFile,
		uint32_t * pixelMask, double energy_in_keV,
		int nimg_pedeG0, int nimg_pedeG1, int nimg_pedeG2,
		int in_start_pos, uint16_t in_summation, int in_start_channel, int in_even_module)
{

	start_pos = in_start_pos;
	start_channel = in_start_channel;
	even_module = in_even_module;
	summation = in_summation;
	dataFileName = data;

	// Create data file
	char filename[1000];
#ifdef NEW_RECEIVER
	sprintf (filename, "%s_f0_0.raw", data.c_str ());
#else
	sprintf (filename, "%s_0.raw", data.c_str ());
#endif
	thisfile = new JungfrauFile (filename, 0);

	// Read Pixel Mask from file
	LoadPixelMask(pixelMaskFile);

	// Read Gain from file
	LoadGain(gainMapFile,energy_in_keV);

	// Calculate Pedestal
	CalculatePede(pedeG1, 1, nimg_pedeG1);
	CalculatePede(pedeG2, 3, nimg_pedeG2);
	CalculatePede(pedeG0, pedeG0_RMS, 0, nimg_pedeG0);

	// mask pixels with high RMS in G0
	for (int i = 0; i < CHANNELS_PER_WORKER; i++) {
		if (pedeG0_RMS[i] > RMS_THRESHOLD) channelMask[i] |= 4;
	}

	images_read = nimg_pedeG0 + nimg_pedeG1 + nimg_pedeG2;
	frames_written = 0;
	curr_file = 0;

	// Ignore frames with no photons
	for (int i = 0; i < px.frames_delay; i ++) {
		if (thisfile->ReadFrame (jf_frame, start_channel)) {
#ifdef PEDESTAL_DRIFT
			for (int line = 0 ; line < CHANNELS_PER_WORKER/1024 ; line ++) {
#pragma ivdep
				for (int j = 0; j < 1024; j++)
				{
					int ch = j + line*1024;

//					uint16_t gain = (jf_frame[ch] & 0xc000);
//					uint16_t adc = jf_frame[ch] & 0x3fff;

					float adcfloat = jf_frame[ch]; // only gain0, so jf_frame[ch] == adc
					if ((jf_frame[ch] < 0xc000) && (fabs (adcfloat - pedeG0[ch]) < TIMES_RMS_THRESHOLD_FOR_PEDE * pedeG0_RMS[ch]))
						pedeG0[ch] = (pedeG0[ch] * (PEDE_WINDOW - 1.0f) + adcfloat) / PEDE_WINDOW;
				}
			}
#endif
		} else std::cout << "\nFile error!\n";
		ReloadFile (1);
	}

	// Copy mask to common buffer
	// Accounts for multi-size pixels
	for (int i = 0; i < CHANNELS_PER_WORKER; i++) {
		int ch = start_channel + i;
		int pixelid = start_pos + (ch / 1024) * XPIXEL + ((ch % 1024) / 256) * 258 + (ch % 256);

		// Mirror flip
		pixelid = (YPIXEL - 1 - pixelid / XPIXEL) * XPIXEL +  pixelid % XPIXEL;

		pixelMask[pixelid] = channelMask[i];

		if ((ch%256 == 0) && (ch%1024 != 0)) pixelMask[pixelid - 1 ] = channelMask[i];
		if ((ch%256 == 255) && (ch%1024 != 1023)) pixelMask[pixelid + 1 ] = channelMask[i];
		if (NCH == 256*1024) { // 2 kHz
                    if (even_module) {
		        if (ch/1024 == 0) pixelMask[pixelid - YPIXEL ] = channelMask[i];
		    } else {
		        if (ch/1024 == 255) pixelMask[pixelid + YPIXEL ] = channelMask[i];
                    }
		} else {
		    if (ch/1024 == 255) pixelMask[pixelid + YPIXEL ] = channelMask[i];
		    if (ch/1024 == 256) pixelMask[pixelid - YPIXEL ] = channelMask[i];
		}
	}
};

JungfrauConverter::~JungfrauConverter ()
{
	delete thisfile;
}

void JungfrauConverter::LoadPixelMask(std::string pixelMaskFile) {

	std::fill_n(channelMask, CHANNELS_PER_WORKER, 0); // initialise with zeros

	if (pixelMaskFile != "") {
		std::ifstream inPixelMask(pixelMaskFile.c_str());
		if (inPixelMask.is_open())
		{
			for (int i = 0; i < start_channel; i++) {
				int tmp;
				inPixelMask >> tmp;
			}
			for (int i = 0; i < CHANNELS_PER_WORKER; i++) {
				int tmp;
				inPixelMask >> tmp;
				if (tmp == 0) channelMask[i] |= 8;
			}
			inPixelMask.close();
		}
		else std::cout << "Unable to open file " << pixelMaskFile << std::endl;
	}
}

void JungfrauConverter::LoadGain(std::string gainMapFile, double energy_in_keV) {
	// Gain maps (double is found in the original file)
	double *tmp_gainG0 = (double *) malloc(CHANNELS_PER_WORKER*sizeof(double));
	MALLOC_ERROR(tmp_gainG0);

	double *tmp_gainG1 = (double *) malloc(CHANNELS_PER_WORKER*sizeof(double));
	MALLOC_ERROR(tmp_gainG1);

	double *tmp_gainG2 = (double *) malloc(CHANNELS_PER_WORKER*sizeof(double));
	MALLOC_ERROR(tmp_gainG2);

	std::fstream infile;
	infile.open (gainMapFile.c_str (), infile.in | infile.binary);
	if (!infile.is_open ())
	{
		std::cout << "File " << gainMapFile << " not found!" << std::endl;
		return;
	}

	size_t byte_ignore0;
	size_t byte_ignore1;

	if (NCH == 256*1024) {
		byte_ignore0 = (even_module * NCH + start_channel)*sizeof(double);
		byte_ignore1 = byte_ignore0 + (2*NCH - even_module * NCH - CHANNELS_PER_WORKER - start_channel)*sizeof(double);
	} else {
		byte_ignore0 = start_channel*sizeof(double);
		byte_ignore1 = byte_ignore0 + (NCH - CHANNELS_PER_WORKER - start_channel)*sizeof(double);
	}

	infile.seekg(byte_ignore0, infile.cur);
	infile.read ((char *) tmp_gainG0, CHANNELS_PER_WORKER * sizeof (double));

	infile.seekg(byte_ignore1, infile.cur);
	infile.read ((char *) tmp_gainG1, CHANNELS_PER_WORKER * sizeof (double));

	infile.seekg(byte_ignore1, infile.cur);
	infile.read ((char *) tmp_gainG2, CHANNELS_PER_WORKER * sizeof (double));

	infile.close ();

	for (int i = 0; i < CHANNELS_PER_WORKER; i++)
	{
		gainG0[i] = 1 / (tmp_gainG0[i] * energy_in_keV);
		gainG1[i] = 1 / (tmp_gainG1[i] * energy_in_keV);
		gainG2[i] = 1 / (tmp_gainG2[i] * energy_in_keV);
	}

	free(tmp_gainG0);
	free(tmp_gainG1);
	free(tmp_gainG2);
}

// Calculates pedestal
void JungfrauConverter::CalculatePede(float *pede, uint16_t expected, uint32_t size) {
	double *tmp_pede  = (double *) malloc(CHANNELS_PER_WORKER*sizeof(double));
	MALLOC_ERROR(tmp_pede);

	for (int j = 0; j < CHANNELS_PER_WORKER; j++)
		tmp_pede[j] = 0.0;

	// plain average for 0...PEDE_WINDOW
	for (int i = 0; i < PEDE_WINDOW; i++) {
		thisfile->ReadFrame(jf_frame, start_channel);

		for (int j = 0; j < CHANNELS_PER_WORKER; j++)
		{
			uint16_t gain = (jf_frame[j] & 0xc000) >> 14;
			if (gain != expected) channelMask[j] |= 2;
			else tmp_pede[j] += (double) (jf_frame[j] & 0x3fff);
		}
	}

	// moving average for PEDE_WINDOW...
	for (int i = 0; i < size-PEDE_WINDOW; i++) {
		thisfile->ReadFrame(jf_frame, start_channel);
		for (int j = 0; j < CHANNELS_PER_WORKER; j++)
		{
			uint16_t gain = (jf_frame[j] & 0xc000) >> 14;
			if (gain != expected) channelMask[j] |= 2;
			else tmp_pede[j] += (double) (jf_frame[j] & 0x3fff)  - tmp_pede[j]/PEDE_WINDOW;
		}
	}

	for (int j = 0; j < CHANNELS_PER_WORKER; j++) {
		pede[j] = tmp_pede[j]/PEDE_WINDOW;
	}
	free(tmp_pede);
}

// Calculates pedestal and it standard deviation
void JungfrauConverter::CalculatePede(float *pede, float *pedeRMS, uint16_t expected, uint32_t size) {
	double *tmp_pede  = (double *) malloc(CHANNELS_PER_WORKER*sizeof(double));
	MALLOC_ERROR(tmp_pede);

	double *tmp_pede2  = (double *) malloc(CHANNELS_PER_WORKER*sizeof(double));
	MALLOC_ERROR(tmp_pede2);

	for (int j = 0; j < CHANNELS_PER_WORKER; j++) {
		tmp_pede[j] = 0.0;
		tmp_pede2[j] = 0.0;
	}

	for (int i = 0; i < PEDE_WINDOW; i++) {
		thisfile->ReadFrame(jf_frame, start_channel);
		for (int j = 0; j < CHANNELS_PER_WORKER; j++)
		{
			uint16_t gain = (jf_frame[j] & 0xc000) >> 14;
			if (gain != expected) channelMask[j] |= 2;
			else {
				double adc = jf_frame[j] & 0x3fff;
				tmp_pede[j] += adc;
				tmp_pede2[j] += adc*adc;
			}
		}
	}

	for (int i = 0; i < size-PEDE_WINDOW; i++) {
		thisfile->ReadFrame(jf_frame, start_channel);
		for (int j = 0; j < CHANNELS_PER_WORKER; j++)
		{
			uint16_t gain = (jf_frame[j] & 0xc000) >> 14;
			if (gain != expected) channelMask[j] |= 2;
			else {
				double adc = jf_frame[j] & 0x3fff;
				tmp_pede[j] += adc - tmp_pede[j]/PEDE_WINDOW;
				tmp_pede2[j] += adc*adc - tmp_pede2[j]/PEDE_WINDOW;
			}
		}
	}

	for (int j = 0; j < CHANNELS_PER_WORKER; j++) {
		pede[j] = tmp_pede[j]/PEDE_WINDOW;
		// variance = <x^2> - <x>^2
		double tmp_variance = (tmp_pede2[j]/PEDE_WINDOW - pede[j]*pede[j]);
		pedeRMS[j] = sqrt(tmp_variance);
	}

	free(tmp_pede);
	free(tmp_pede2);
}

// Copies single line
void JungfrauConverter::CopyLine(int32_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_i32(localBuffer[start_pos_in_local+i] / 2.0f);
			imageBuffer[start_pos_in_buffer+idx+1] = jf_round_i32(localBuffer[start_pos_in_local+i] / 2.0f);
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_i32(localBuffer[start_pos_in_local+i]);
			idx++;
		}
	}
}

void JungfrauConverter::CopyLine(uint32_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_u32(localBuffer[start_pos_in_local+i] / 2.0f);
			imageBuffer[start_pos_in_buffer+idx+1] = jf_round_u32(localBuffer[start_pos_in_local+i] / 2.0f);
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_u32(localBuffer[start_pos_in_local+i]);
			idx++;
		}
	}
}

void JungfrauConverter::CopyLine(uint16_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_u16(localBuffer[start_pos_in_local+i] / 2.0f);
			imageBuffer[start_pos_in_buffer+idx+1] = jf_round_u16(localBuffer[start_pos_in_local+i] / 2.0f);
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_u16(localBuffer[start_pos_in_local+i]);
			idx++;
		}
	}
}

void JungfrauConverter::CopyLine(int16_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_i16(localBuffer[start_pos_in_local+i] / 2.0f);
			imageBuffer[start_pos_in_buffer+idx+1] = jf_round_i16(localBuffer[start_pos_in_local+i] / 2.0f);
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_i16(localBuffer[start_pos_in_local+i]);
			idx++;
		}
	}
}

// Copies single line
void JungfrauConverter::CopyLine(float *imageBuffer, int start_pos_in_buffer, int start_pos_in_local) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = localBuffer[start_pos_in_local+i] / 2.0f;
			imageBuffer[start_pos_in_buffer+idx+1] = localBuffer[start_pos_in_local+i] / 2.0f;
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = localBuffer[start_pos_in_local+i];
			idx++;
		}
	}
}

// Copies single line
void JungfrauConverter::CopyLineWithHorizGap(int32_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_i32(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+1] = jf_round_i32(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = jf_round_i32(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL+1] = jf_round_i32(localBuffer[start_pos_in_local+i] / 4.0f);
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_i32(localBuffer[start_pos_in_local+i] / 2.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = jf_round_i32(localBuffer[start_pos_in_local+i] / 2.0f);
			idx++;
		}
	}
}

// Copies single line
void JungfrauConverter::CopyLineWithHorizGap(uint32_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_u32(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+1] = jf_round_u32(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = jf_round_u32(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL+1] = jf_round_u32(localBuffer[start_pos_in_local+i] / 4.0f);
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_u32(localBuffer[start_pos_in_local+i] / 2.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = jf_round_u32(localBuffer[start_pos_in_local+i] / 2.0f);
			idx++;
		}
	}
}

// Copies single line
void JungfrauConverter::CopyLineWithHorizGap(uint16_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_u16(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+1] = jf_round_u16(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = jf_round_u16(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL+1] = jf_round_u16(localBuffer[start_pos_in_local+i] / 4.0f);
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_u16(localBuffer[start_pos_in_local+i] / 2.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = jf_round_u16(localBuffer[start_pos_in_local+i] / 2.0f);
			idx++;
		}
	}

}

// Copies single line
void JungfrauConverter::CopyLineWithHorizGap(int16_t *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_i16(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+1] = jf_round_i16(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = jf_round_i16(localBuffer[start_pos_in_local+i] / 4.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL+1] = jf_round_i16(localBuffer[start_pos_in_local+i] / 4.0f);
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = jf_round_i16(localBuffer[start_pos_in_local+i] / 2.0f);
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = jf_round_i16(localBuffer[start_pos_in_local+i] / 2.0f);
			idx++;
		}
	}

}

// Copies single line
void JungfrauConverter::CopyLineWithHorizGap(float *imageBuffer, int start_pos_in_buffer, int start_pos_in_local, int gap) {
	int idx = 0;
#pragma ivdep
	for (int i = 0; i < 1024; i++) {
		if (((i % 256 == 255) && (i / 256 < 3)) || (i % 256 == 0) && (i / 256 > 0)) {
			imageBuffer[start_pos_in_buffer+idx] = localBuffer[start_pos_in_local+i] / 4.0f;
			imageBuffer[start_pos_in_buffer+idx+1] = localBuffer[start_pos_in_local+i] / 4.0f;
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = localBuffer[start_pos_in_local+i] / 4.0f;
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL+1] = localBuffer[start_pos_in_local+i] / 4.0f;
			idx += 2;
		} else {
			imageBuffer[start_pos_in_buffer+idx] = localBuffer[start_pos_in_local+i] / 2.0f;
			imageBuffer[start_pos_in_buffer+idx+gap*XPIXEL] = localBuffer[start_pos_in_local+i] / 2.0f;
			idx++;
		}
	}

}

void JungfrauConverter::ReloadFile (int summation)
{
	images_read++;
	if (images_read % summation == 0)
		frames_written++;
	if (images_read % IMAGES_IN_FILE == 0)
	{
		curr_file++;
		char filename[1000];
		//TODO: Keep old settings for 1 kHz
#ifdef NEW_RECEIVER
		sprintf (filename, "%s_f%d_0.raw",
				dataFileName.c_str(),
				curr_file);
#else
		sprintf (filename, "%s%012d_0.raw",
				dataFileName.substr (0, dataFileName.length () - 12).c_str (),
				10000 * curr_file);

#endif
		//std::cout << filename << std::endl;
		delete thisfile;
		thisfile = new JungfrauFile(filename, curr_file);
	}
}

template<typename T> void JungfrauConverter::Convert (T *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write) {
	float adcfloat;
	//TODO - if summation > FRAME_BUF_SIZE
	for (int i = 0; i < frames_to_write; i ++)
	{
		int local_G1 = 0;
		int local_G2 = 0;

		bool loaded_ok = true;
		for (int k = 0; k < summation; k++) {
			if (! thisfile->ReadFrame (&(jf_frame_buffer[k*CHANNELS_PER_WORKER]), start_channel)) loaded_ok = false;
			ReloadFile (1);
		}

		for (int line = 0 ; line < CHANNELS_PER_WORKER/1024 ; line ++) {
			uint32_t pixelid0 = XPIXEL * YPIXEL * i;
			uint32_t line_in_module = line + start_channel/1024;

			// Mirror image of the detector
			uint32_t pos_img = start_pos + line_in_module*XPIXEL;
			uint32_t pos = pixelid0 + (YPIXEL - pos_img / XPIXEL - 1) * XPIXEL + pos_img % XPIXEL;
#pragma ivdep
			for (int j = 0; j < 1024; j++)
			{
				int ch = j + line*1024;

				uint16_t gain = jf_frame_buffer[ch] & 0xc000;
				uint16_t adc  = jf_frame_buffer[ch] & 0x3fff;

				adcfloat = adc;
				switch (gain)
				{
				case 0:
#ifdef PEDESTAL_DRIFT
					if (fabs (adcfloat - pedeG0[ch]) < TIMES_RMS_THRESHOLD_FOR_PEDE * pedeG0_RMS[ch])
						pedeG0[ch] = (pedeG0[ch] * ((float) (PEDE_WINDOW-1.0)) + adcfloat) / ((float) (PEDE_WINDOW));
#endif
					localBuffer[j] = (adcfloat - pedeG0[ch]) * gainG0[ch];
					break;
				case 0x4000:
					localBuffer[j] = (adcfloat - pedeG1[ch]) * gainG1[ch];
					if (channelMask[ch] == 0) local_G1++;
					break;
				case 0xc000:
					localBuffer[j] = (adcfloat - pedeG2[ch]) * gainG2[ch];
					if (channelMask[ch] == 0) local_G2++;
					if (adc == 0) localBuffer[j] = 1.0f / 0.0f;
					else if (adc == 0x3fff) localBuffer[j] = 1.0f / 0.0f;
					break;
				default:
					localBuffer[j] = 1.0f / 0.0f;
					break;
				}
			}

			for (int k = 1; k < summation; k++) {
#pragma ivdep
				for (int j = 0; j < 1024; j++)
				{

					int ch = j + line*1024;

					uint16_t gain = jf_frame_buffer[ch + k*CHANNELS_PER_WORKER] & 0xc000;
					uint16_t adc = jf_frame_buffer[ch+k*CHANNELS_PER_WORKER] & 0x3fff;

					adcfloat = adc;
					switch (gain)
					{
					case 0:
#ifdef PEDESTAL_DRIFT
						if (fabs (adcfloat - pedeG0[ch]) < TIMES_RMS_THRESHOLD_FOR_PEDE * pedeG0_RMS[ch])
							pedeG0[ch] = (pedeG0[ch] * ((float) (PEDE_WINDOW-1.0)) + adcfloat) / ((float) (PEDE_WINDOW));
#endif
						localBuffer[j] += (adcfloat - pedeG0[ch]) * gainG0[ch];
						break;
					case 0x4000:
						localBuffer[j] += (adcfloat - pedeG1[ch]) * gainG1[ch];
						if (channelMask[ch] == 0) local_G1++;
						break;
					case 0xc000:
						localBuffer[j] += (adcfloat - pedeG2[ch]) * gainG2[ch];
						if (channelMask[ch] == 0) local_G2++;
						if (adc == 0) localBuffer[j] = 1.0f / 0.0f;
						else if (adc == 0x3fff) localBuffer[j] = 1.0f / 0.0f;
						break;
					default:
						localBuffer[j] = 1.0f / 0.0f;
						break;
					}
				}

			}
			if (NCH == 256*1024) { // 2 kHz
				if (even_module) {
					if (line_in_module == 0)      CopyLineWithHorizGap(imageBuffer, pos , 0,+1);
					else                          CopyLine(imageBuffer, pos , 0);
				} else {
					if (line_in_module == 255)    CopyLineWithHorizGap(imageBuffer, pos , 0,-1);
					else                          CopyLine(imageBuffer, pos , 0);
				}
			} else {
				if (line_in_module > 255) pos -= 2*XPIXEL; // account for horizontal line with multi-sized pixels

				if (line_in_module == 255)	 CopyLineWithHorizGap(imageBuffer, pos , 0,-1);
				else if (line_in_module == 256)  CopyLineWithHorizGap(imageBuffer, pos , 0,+1);
				else                             CopyLine(imageBuffer, pos , 0);
			}
		}
		pixels_in_G1[i] = local_G1;
		pixels_in_G2[i] = local_G2;
	}
};

template<typename T> void JungfrauConverter::ConvertSum1 (T *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write) {
	float adcfloat;

	//TODO - if frames_to_write > FRAME_BUF_SIZE
	for (int k = 0; k < frames_to_write; k++) {
		bool loaded_ok = true;
		if (! thisfile->ReadFrame (&(jf_frame_buffer[k*CHANNELS_PER_WORKER]), start_channel)) loaded_ok = false;
		ReloadFile (1);
	}

	for (int line = 0 ; line < CHANNELS_PER_WORKER/1024 ; line ++) {

		uint32_t line_in_module = line + start_channel/1024;
		for (int i = 0; i < frames_to_write; i++) {

			int local_G1 = 0;
			int local_G2 = 0;

			uint32_t pixelid0 = XPIXEL * YPIXEL * i;

			// Mirror image of the detector
			uint32_t pos_img = start_pos + line_in_module*XPIXEL;
			uint32_t pos = pixelid0 + (YPIXEL - pos_img / XPIXEL - 1) * XPIXEL + pos_img % XPIXEL;

#pragma ivdep
			for (int j = 0; j < 1024; j++) {
				int ch = j + line*1024;

				uint16_t gain =  jf_frame_buffer[ch + i*CHANNELS_PER_WORKER] & 0xc000;
				uint16_t adc =   jf_frame_buffer[ch + i*CHANNELS_PER_WORKER] & 0x3fff;

				adcfloat = adc;
				switch (gain)
				{
				case 0:
#ifdef PEDESTAL_DRIFT
					if (fabs (adcfloat - pedeG0[ch]) < TIMES_RMS_THRESHOLD_FOR_PEDE * pedeG0_RMS[ch])
						pedeG0[ch] = (pedeG0[ch] * ((float) (PEDE_WINDOW-1.0)) + adcfloat) / ((float) (PEDE_WINDOW));
#endif
					localBuffer[j] = (adcfloat - pedeG0[ch]) * gainG0[ch];
					break;
				case 0x4000:
					localBuffer[j] = (adcfloat - pedeG1[ch]) * gainG1[ch];
					if (channelMask[ch] == 0) local_G1++;
					break;
				case 0xc000:
					localBuffer[j] = (adcfloat - pedeG2[ch]) * gainG2[ch];
					if (channelMask[ch] == 0) local_G2++;
					if (adc == 0) localBuffer[j] = 1.0f / 0.0f;
					else if (adc == 0x3fff) localBuffer[j] = 1.0f / 0.0f;
					break;
				default:
					localBuffer[j] = 1.0f / 0.0f;
					break;
				}
			}
			if (NCH == 256*1024) { // 2 kHz
				if (even_module) {
					if (line_in_module == 0)      CopyLineWithHorizGap(imageBuffer, pos , 0,+1);
					else                          CopyLine(imageBuffer, pos , 0);
				} else {
					if (line_in_module == 255)    CopyLineWithHorizGap(imageBuffer, pos , 0,-1);
					else                          CopyLine(imageBuffer, pos , 0);
				}
			} else {
				if (line_in_module > 255) pos -= 2*XPIXEL; // account for horizontal line with multi-sized pixels

				if (line_in_module == 255)	 CopyLineWithHorizGap(imageBuffer, pos , 0,-1);
				else if (line_in_module == 256)  CopyLineWithHorizGap(imageBuffer, pos , 0,+1);
				else                             CopyLine(imageBuffer, pos , 0);
			}
			pixels_in_G1[i] += local_G1;
			pixels_in_G2[i] += local_G2;
		}
	}
};

template void JungfrauConverter::ConvertSum1<float> (float *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
template void JungfrauConverter::ConvertSum1<uint16_t> (uint16_t *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
template void JungfrauConverter::ConvertSum1<int16_t> (int16_t *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
template void JungfrauConverter::ConvertSum1<uint32_t> (uint32_t *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
template void JungfrauConverter::ConvertSum1<int32_t> (int32_t *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);

template void JungfrauConverter::Convert<float> (float *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
template void JungfrauConverter::Convert<int32_t> (int32_t *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
template void JungfrauConverter::Convert<uint16_t> (uint16_t *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
template void JungfrauConverter::Convert<int16_t> (int16_t *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
template void JungfrauConverter::Convert<uint32_t> (uint32_t *imageBuffer, std::vector<int> &pixels_in_G1, std::vector<int> &pixels_in_G2, long frames_to_write);
