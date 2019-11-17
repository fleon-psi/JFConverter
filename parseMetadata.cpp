#include <pthread.h>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstdint>
#include <string>
#include <sched.h>
#include "nlohmann/json.hpp"
#include <getopt.h>

#include "jungfrau.h"
#include "zstd.h"

void show_usage() {
	std::cerr << std::endl;
	std::cerr << "Usage ./JFConverter -j<json file> -s|-r|-d {...}" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  Always specify:" << std::endl;
	std::cerr << "    -j<json file>         JSON file" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  And one option among the following:" << std::endl;
	std::cerr << "    -s                    Read raw images from SSD  (/mnt/ssd)" << std::endl;
	std::cerr << "    -r                    Read raw images from RAM  (/mnt/n*ram)" << std::endl;
	std::cerr << "    -d                    Read raw images from disk (/mnt/zfs)" << std::endl;
	std::cerr << "    -d<folder>            Read raw images from disk (<folder>)" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  Optional parameters:" << std::endl;
	std::cerr << "    -l<delay>             Ignore <delay> first frames" << std::endl;
	std::cerr << "    -f<frames>            Generate <frames> frames" << std::endl;
	std::cerr << "    -h<suffix>            Add suffix for the HDF5 filenames" << std::endl;
	std::cerr << "    -S<summation>         Specify summation amount" << std::endl;
	std::cerr << "    -o<folder>            Save results in a <folder> (/mnt/zfs/output as default)" << std::endl;
	std::cerr << "    -c<compression>       Compress data (bshuf_lz4 default)" << std::endl;
	std::cerr << "    -m<multiply>          Round to multiply of photon counts" << std::endl;
	std::cerr << "    -M<multiply>          Round to multiply of energy in eV" << std::endl;
	std::cerr << "    -W<threads>           Run <threads> writers (" << MAX_NUM_WRITERS << ")" << std::endl;
	std::cerr << "    -B<size>              Set BSHUF block size to <size>" << std::endl;
	std::cerr << "    -N                    Don't use direct chunk writer"    << std::endl;
	//std::cerr << "    -a                    Write unconverted ADC (slow and doesn't compress very well!)" << std::endl;
	std::cerr << "    -x                    Output as floating point numbers" << std::endl;
	std::cerr << "    -u                    Output as unsigned numbers" << std::endl;
	std::cerr << std::endl;
	std::cerr << "  Compression algorithms available: gzip, lz4, bshuf_lz4, sz, zstd, bshuf_zstd, none" << std::endl;
	std::cerr << "  Prefix fast_ indicates direct chunk writer" << std::endl;
	std::cerr << std::endl;


}

int parseMetadata(PXMetadata &px, int argc, char **argv) {
	bool parsed_ok = true;

	bool json_opened = false;
	bool use_direct_chunk_writer = true;
	bool use_signed = true;
	bool use_float = false;
	int c;
	int option_index = 0;
	std::string json_filename;
	std::string path_node0, path_node1, path_node2, path_node3;
	std::string suffix = "";
	int number_of_writers = MAX_NUM_WRITERS;
        
#if defined(HPEDL580)
	std::string output = "/mnt/zfs/output";
#else
	std::string output = "";
#endif
	int sum = -1;
	px.sumFloat = false;
	int frames = -1;
	int delay = -1;

	bool still = false;
	float mult = -1.0;
	float multE = -1.0;

	ConverterCompression compressionAlgorithm = COMPRESSION_BSHUF_LZ4;
	double absolute_err_bound_sz = 1.0;
	int gzip_level = -1;

	static struct option long_options[] = {
			{"disk"  ,    optional_argument,  0, 'd'},
			{"ssd"   ,          no_argument,  0, 's'},
			{"ram"   ,          no_argument,  0, 'r'},
			{"json"  ,    required_argument,  0, 'j'},
			{"frames",    required_argument,  0, 'f'},
			{"delay",     required_argument,  0, 'l'},
			{"hdf5suffix",required_argument,  0, 'h'},
			{"sum"   ,    required_argument,  0, 'S'},
			{"output",    required_argument,  0, 'o'},
			{"nopede",          no_argument,  0, 'n'},
			{"mult"  ,    required_argument,  0, 'm'},
			{"multE"  ,   required_argument,  0, 'M'},
			{"debug" ,          no_argument,  0, 'g'},
			{"uint",            no_argument,  0, 'u'},
			{"float",           no_argument,  0, 'x'},
			{"compression",required_argument, 0, 'c'},
			{"gzip",       required_argument, 0, 'Z'},
			{"still",      required_argument, 0, 't'},
			{"abs",        required_argument, 0, 'A'},
			{"nodirect",        no_argument,  0, 'N'},
			{"bshuf_block",required_argument, 0, 'B'},
                        {"writer",     required_argument, 0, 'W'},
			{0       ,                    0,  0,  0 }
	};
	int tmp;

	while((c = getopt_long (argc, argv, "dsrl:f:j:h:t:m:M:h:S:uo:ngc:ixZ:tA:NB:", long_options, &option_index)) != -1) {
		switch(c) {
		case 'W':
			number_of_writers = atoi(optarg);
			if ((number_of_writers <= 0) || (number_of_writers > MAX_NUM_WRITERS)) {
                                std::cerr << "ERROR: Number of writers must be positive number and less than " << MAX_NUM_WRITERS << std::endl;
                                parsed_ok = false;
                        }
			break;	
		case 'x':
			use_float = true;
			break;
		case 'u':
			use_signed = false;
			break;
		case 'B':
			tmp = atoi(optarg);
			if (tmp < 0) {
				std::cerr << "ERROR: Bit-shuffle block size has to be > 0" << std::endl;
				parsed_ok = false;
			}
			px.bshuf_block_size = tmp;
			break;
		case 'f':
			frames = atoi(optarg);
			if (frames <= 0) {
				std::cerr << "ERROR: Frame # must be > 0" << std::endl;
				parsed_ok = false;
			}
			break;
		case 'l':
			delay = atoi(optarg);
			if (delay < 0) {
				std::cerr << "ERROR: Frame delay # must be >= 0" << std::endl;
				parsed_ok = false;
			}
			break;
		case 'c':
			if (std::string(optarg) == "bshuf_lz4") compressionAlgorithm = COMPRESSION_BSHUF_LZ4;
			else if (std::string(optarg) == "gzip") compressionAlgorithm = COMPRESSION_GZIP;
			else if (std::string(optarg) == "lz4") compressionAlgorithm = COMPRESSION_LZ4;
			else if (std::string(optarg) == "sz") compressionAlgorithm = COMPRESSION_SZ;
			else if (std::string(optarg) == "zstd") compressionAlgorithm = COMPRESSION_ZSTD;
			else if (std::string(optarg) == "bshuf_zstd") {
				compressionAlgorithm = COMPRESSION_BSHUF_ZSTD;
				std::cerr << "WARNING: Bitshuffle/zstd fully experimental compression mode and currently not implemented in the mainstream bitshuffle." << std::endl;
			}
			else if (std::string(optarg) == "none") compressionAlgorithm = COMPRESSION_NONE;
			else {
				std::cerr << "ERROR: Unknown compression " << optarg << std::endl;
				parsed_ok = false;
			}
			break;
		case 'N':
			use_direct_chunk_writer = false;
			break;
		case 'd':
			if (optarg == NULL) {
				path_node0 = std::string("/mnt/zfs");
				path_node1 = std::string("/mnt/zfs");
				path_node2 = std::string("/mnt/zfs");
				path_node3 = std::string("/mnt/zfs");
			} else {
				path_node0 = std::string(optarg);
				path_node1 = std::string(optarg);
				path_node2 = std::string(optarg);
				path_node3 = std::string(optarg);
			}
			break;
		case 's':
			path_node0 = std::string("/mnt/ssd");
			path_node1 = std::string("/mnt/ssd");
			path_node2 = std::string("/mnt/ssd");
			path_node3 = std::string("/mnt/ssd");
			break;

		case 'r':
#ifdef AC922
			path_node0 = std::string("/mnt/n0ram");
			path_node1 = std::string("/mnt/n1ram");
			path_node2 = std::string("/mnt/n0ram");
			path_node3 = std::string("/mnt/n1ram");
#else
			path_node0 = std::string("/mnt/n0ram");
			path_node1 = std::string("/mnt/n1ram");
			path_node2 = std::string("/mnt/n2ram");
			path_node3 = std::string("/mnt/n3ram");
#endif
			break;
		case 'j':
			json_filename = std::string(optarg);
			break;
		case 'h':
			suffix = std::string(optarg);
			break;
		case 'o':
			output = std::string(optarg);
			break;
		case 'Z':
			gzip_level = atoi(optarg);
			if ((gzip_level < 0) || (gzip_level > 9)) {
				std::cerr << "ERROR: Gzip level must be > 0" << std::endl;
				parsed_ok = false;
			}
			break;
		case 'S':
			sum = atoi(optarg);
			if ((sum <= 0) || (sum > FRAME_BUF_SIZE)) {
				std::cerr << "ERROR: Summation amount must be > 0 and <= " + std::to_string(FRAME_BUF_SIZE)<< std::endl;
				parsed_ok = false;
			}
			break;
		case 'M':
			multE = atof(optarg);
			if (multE <= 0.0) {
				std::cerr << "ERROR: Rounding multiplier (energy) needs to be > 0 " << std::endl;
				parsed_ok = false;
			}
			break;
		case 'm':
			mult = atof(optarg);
			if (mult <= 0.0) {
				std::cerr << "ERROR: Rounding multiplier (photons) needs to be > 0 " << std::endl;
				parsed_ok = false;
			}
			break;
		case 't':
			still = true;
			break;
		case 'A':
			absolute_err_bound_sz = atof(optarg);
			if (absolute_err_bound_sz <= 0.0) {
				std::cerr << "ERROR: Absolute error bound for SZ needs to be > 0 " << std::endl;
				parsed_ok = false;
			}
		default :
			break;
		}
	}

	if (path_node0 == "") {
		std::cerr << "ERROR: One of -s, -r or -d must be specified for the path." << std::endl;
		parsed_ok = false;
	}


	nlohmann::json j;
	std::ifstream input;
	if (json_filename == "") {
		std::cerr << "ERROR: JSON filename must be provided" << std::endl;
		parsed_ok = false;
	}

	input.open(json_filename.c_str());
	if (input.fail()) {
		std::cerr << "ERROR: Error opening/reading file " << json_filename << std::endl;
		parsed_ok = false;
	}

	px.number_of_writers = number_of_writers;
	input >> j;

	nlohmann::json detector = j.at("detector");
	nlohmann::json filewriter = j.at("filewriter");
	nlohmann::json jfr = j.at("jungfrau");

	if (delay != -1) px.frames_delay = delay;
	else {
		if (jfr.find("bl_delay") != jfr.end()) {
			int json_bl_delay = jfr.at("bl_delay");
			if (json_bl_delay == 0) px.frames_delay = 0;
			else if (json_bl_delay < PEDE_G0) {
				parsed_ok = false;
				std::cerr << "ERROR: Trigger during pedestal data collection" << std::endl;
			}
			else px.frames_delay = json_bl_delay - PEDE_G0;
		} else px.frames_delay = 0;
	}

	px.frames_pedeG0 = PEDE_G0;
	px.frames_pedeG1 = jfr.at("framesG1");
	px.frames_pedeG2 = jfr.at("framesG2");

	px.stillImage = still;

	try {
		px.transmission = detector.at("transmission");
	} catch (nlohmann::json::exception& e) {
		px.transmission = -1.0;
	}

	int json_framesG0 = jfr.at("framesG0");
	px.nimages_collected= json_framesG0 - px.frames_pedeG0;

	if (px.frames_delay > px.nimages_collected) {
		std::cerr << "ERROR: Delay longer than actually collected images, check trigger settings" << std::endl;
		parsed_ok = false;
	}

	double speed = ((double)detector.at("omega_increment")) / ((double)detector.at("count_time"));
	if (sum > 0 ) px.summation = sum;
	else if (speed >= 100.0) px.summation = 1;
	else px.summation = std::lround(100.0/speed);
	int bl_summation = jfr.at("summation");

	if (frames > 0) px.intended_number_of_images = frames;
	else {
   	        int bl_nimages = detector.at("nimages");
	        px.intended_number_of_images = bl_summation*bl_nimages/px.summation;
	}



	px.nimages_per_file = filewriter.at("nimages_per_file");
	if ((px.nimages_per_file <= 0) || (px.nimages_per_file > FRAME_BUF_SIZE)) {
		std::cerr << "ERROR: Number of images per data HDF5 file must be > 0 and <= " + std::to_string(FRAME_BUF_SIZE)<< std::endl;
		parsed_ok = false;
	}

	std::string file = filewriter.at("name_pattern");
	if (output != "")
		px.hdf5_prefix = output +"/" + file; // prefix for HDF5 files
	else 
		px.hdf5_prefix = file; // prefix for HDF5 files

	if (suffix != "") px.hdf5_prefix += "_" + suffix;

	std::cout << "Summation of " << px.summation << std::endl;
	std::cout << "Data should be present in " << px.intended_number_of_images << " images." << std::endl;

	px.detector_distance = detector.at("detector_distance");
	px.beamx = detector.at("beam_center_x");
	px.beamy = detector.at("beam_center_y");

	// px.detector_distance = 0.06;
	// px.beamx = 1034;
	// px.beamy = 603;

	px.omega_start = detector.at("omega_start");
	double period;

	try {	
	    px.integration_time = jfr.at("exptime");
 	    period = jfr.at("periodG0");

	} catch (nlohmann::json::exception& e) {
	    px.integration_time = 0.000840;
            period = 0.000880;

	}

	px.frame_time = period * ((double) px.summation);
	px.omega_increment = speed * px.frame_time;

	px.photon_energy = detector.at("photon_energy");
	px.sensor_thickness = detector.at("sensor_thickness");
	px.x_pixel_size = PIXEL_SIZE_IN_MM / 1000.0;
	px.y_pixel_size = PIXEL_SIZE_IN_MM / 1000.0;

	// Compression settings
	px.use_direct_chunk_writer = use_direct_chunk_writer;
	px.compressionAlgorithm = compressionAlgorithm;
	if ((gzip_level == -1) && (compressionAlgorithm == COMPRESSION_ZSTD)) gzip_level = ZSTD_CLEVEL_DEFAULT;
	if ((gzip_level == -1) && (compressionAlgorithm == COMPRESSION_GZIP)) gzip_level = 6;
	px.absolute_err_bound_sz = absolute_err_bound_sz;
	px.gzip_level = gzip_level;

	// Rounding settings
	if ((mult > 0.0) && (multE > 0.0)) {
		std::cerr << "ERROR: Rounding multipliers in photons and energy cannot be used at the same time." << std::endl;
		parsed_ok = false;
	}

	if (mult > 0.0) px.photon_energy_for_normalization = mult*px.photon_energy;
	else if (multE > 0.0) px.photon_energy_for_normalization = multE;
	else px.photon_energy_for_normalization = px.photon_energy;

#ifdef JF_1M_2kHz
#if defined(HPEDL580) || defined(AC922)
	px.dataFile.push_back(path_node3 + "/" + file + "_d0");
	px.dataFile.push_back(path_node3 + "/" + file + "_d1");
	px.dataFile.push_back(path_node2 + "/" + file + "_d2");
	px.dataFile.push_back(path_node2 + "/" + file + "_d3");
#else
	px.dataFile.push_back(file + "_d0");
	px.dataFile.push_back(file + "_d1");
	px.dataFile.push_back(file + "_d2");
	px.dataFile.push_back(file + "_d3");
#endif
#else
#if defined(HPEDL580) || defined(AC922)
	px.dataFile.push_back(path_node2 + "/" + file + "_d0_f000000000000");
	px.dataFile.push_back(path_node2 + "/" + file + "_d1_f000000000000");
	px.dataFile.push_back(path_node0 + "/" + file + "_d2_f000000000000");
	px.dataFile.push_back(path_node0 + "/" + file + "_d3_f000000000000");
	px.dataFile.push_back(path_node3 + "/" + file + "_d4_f000000000000");
	px.dataFile.push_back(path_node3 + "/" + file + "_d5_f000000000000");
	px.dataFile.push_back(path_node1 + "/" + file + "_d6_f000000000000");
	px.dataFile.push_back(path_node1 + "/" + file + "_d7_f000000000000");
#else
	px.dataFile.push_back(file + "_d0_f000000000000");
	px.dataFile.push_back(file + "_d1_f000000000000");
	px.dataFile.push_back(file + "_d2_f000000000000");
	px.dataFile.push_back(file + "_d3_f000000000000");
	px.dataFile.push_back(file + "_d4_f000000000000");
	px.dataFile.push_back(file + "_d5_f000000000000");
	px.dataFile.push_back(file + "_d6_f000000000000");
	px.dataFile.push_back(file + "_d7_f000000000000");
#endif
#endif

        if (use_float) px.mode = FLOAT_MODE;
        else {
            // Approximate value of the highest possible pixel value
            float approx_saturation = 15000.0 * 12400.0 * px.summation / px.photon_energy_for_normalization;
	    if (use_signed) {
               if (approx_saturation > INT16_MAX/2) px.mode = INT32_MODE;
               else px.mode = INT16_MODE;
            } else {
               if (approx_saturation > UINT16_MAX/2) px.mode = UINT32_MODE;
               else px.mode = UINT16_MODE;
            }
        }

	for (int i = 0; i < NMODULES; i++) {
		try {
			nlohmann::json jmod = jfr.at("mod" +std::to_string(i));
			px.tempFPGA[i] = jmod.at("tempFPGA_disarm");
		} catch (nlohmann::json::exception& e) {
			px.tempFPGA[i] = -273.15;
		}
	}

	if (!parsed_ok) {
		show_usage();
		exit(EXIT_FAILURE);
	}
	return 0;
};

