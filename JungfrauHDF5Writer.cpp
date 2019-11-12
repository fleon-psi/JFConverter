#include <iostream>
#include <cstdlib>
#include <math.h>
#include <cstring>

#include "jungfrau.h"

JungfrauHDF5Writer::JungfrauHDF5Writer(const PXMetadata& in_metadata) {
	nchunks_used = 0;
	nimages_saved = 0;
	nimages = in_metadata.intended_number_of_images;
	prefix = in_metadata.hdf5_prefix;
	bshuf_block_size = in_metadata.bshuf_block_size;

	compressionAlgorithm = in_metadata.compressionAlgorithm;
	gzip_level = in_metadata.gzip_level;

	images_per_file = in_metadata.nimages_per_file;

	nchunks = nimages / images_per_file;
	if (nimages % images_per_file != 0) nchunks++;

	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	InitCompression();

	std::string master_file_name = prefix+"_master.h5";
	master_file_id = H5Fcreate(master_file_name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
	if (master_file_id < 0) {
		std::cout << "File error." << std::endl;
		std::exit(1);
	};

	this->SaveMetadata(in_metadata);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

	switch (mode) {
	case UINT16_MODE:
		elem_size = 2;
		data_type = H5T_STD_U16LE;
		break;
	case UINT32_MODE:
		elem_size = 4;
		data_type = H5T_STD_U32LE;
		break;
	case FLOAT_MODE:
		elem_size = 4;
		data_type = H5T_NATIVE_FLOAT;
		break;
	case SIMPLE_MODE:
		elem_size = 4;
		data_type = H5T_STD_I32LE;
		break;
	}

	output_size = CompressionOutputMaxSize();

}

void JungfrauHDF5Writer::WriteMetrology(const PXMetadata& in_metadata)
{
	double moduleOrigin[NMODULES_GEOM][3];
	double detectorCenter[3] = {0.0, 0.0, 0.0};
	for (int i = 0; i < NMODULES_GEOM; i ++) {
		// 1. Find module corner in lab coordinates, based on pixel coordinates
		moduleOrigin[i][0] = (corner[i][0] - in_metadata.beamx) * PIXEL_SIZE_IN_MM;
		moduleOrigin[i][1] = (corner[i][1] - in_metadata.beamy) * PIXEL_SIZE_IN_MM;
		moduleOrigin[i][2] = 0.0; // No Z offset for modules at the moment

		// 2. Calculate vector from module in mm
		double moduleCenter[3];
		if ((slow != NULL) && (fast != NULL)) {
			double fastNorm = sqrt(fast[i][0] * fast[i][0] + fast[i][1] * fast[i][1] + fast[i][2] * fast[i][2]);
			double slowNorm = sqrt(slow[i][0] * slow[i][0] + slow[i][1] * slow[i][1] + slow[i][2] * slow[i][2]);
			moduleCenter[0] = moduleOrigin[i][0] + PIXEL_SIZE_IN_MM * (1030 / 2 * fast[i][0] / fastNorm + 514 / 2 * slow[i][0] / slowNorm);
			moduleCenter[1] = moduleOrigin[i][1] + PIXEL_SIZE_IN_MM * (1030 / 2 * fast[i][1] / fastNorm + 514 / 2 * slow[i][1] / slowNorm);
			moduleCenter[2] = moduleOrigin[i][2] + PIXEL_SIZE_IN_MM * (1030 / 2 * fast[i][2] / fastNorm + 514 / 2 * slow[i][2] / slowNorm);
		} else {
			moduleCenter[0] = moduleOrigin[i][0] + PIXEL_SIZE_IN_MM * 1030 / 2;
			moduleCenter[1] = moduleOrigin[i][1] + PIXEL_SIZE_IN_MM * 514 / 2;
			moduleCenter[2] = moduleOrigin[i][2];
		}

		// 3. Find detector center
		detectorCenter[0] += moduleCenter[0];
		detectorCenter[1] += moduleCenter[1];
		detectorCenter[2] += moduleCenter[2];
	}

	detectorCenter[0] /= NMODULES_GEOM;
	detectorCenter[1] /= NMODULES_GEOM;
	detectorCenter[2] /= NMODULES_GEOM;

	hid_t grp, dataset;

	grp = createGroup(master_file_id, "/entry/instrument/" DETECTOR_NAME "/transformations","NXtransformations");

	saveDouble(grp, "AXIS_RAIL", in_metadata.detector_distance * 1000, "mm");

	double rail_vector[3] = {0,0,1};

	dataset = H5Dopen2(grp, "AXIS_RAIL", H5P_DEFAULT);
	addStringAttribute(dataset, "depends_on", ".");
	addStringAttribute(dataset, "equipment", "detector");
	addStringAttribute(dataset, "equipment_component", "detector_arm");
	addStringAttribute(dataset, "transformation_type", "translation");
	addDoubleAttribute(dataset, "vector", rail_vector, 3);
	H5Dclose(dataset);

	saveDouble(grp, "AXIS_D0", 0.0, "degrees");

	double d0_vector[3] = {0, 0, -1};

	dataset = H5Dopen2(grp , "AXIS_D0", H5P_DEFAULT);
	addStringAttribute(dataset, "depends_on", "AXIS_RAIL");
	addStringAttribute(dataset, "equipment", "detector");
	addStringAttribute(dataset, "equipment_component", "detector_arm");
	addStringAttribute(dataset, "transformation_type", "rotation");
	addDoubleAttribute(dataset, "vector", d0_vector, 3);
	addDoubleAttribute(dataset, "offset", detectorCenter, 3);
	addStringAttribute(dataset, "offset_units", "mm");
	H5Dclose(dataset);

	for (int i = 0; i < NMODULES_GEOM; i++) {
		double mod_vector[3] = {0, 0, -1};
		double mod_offset[3] = {moduleOrigin[i][0] - detectorCenter[0], moduleOrigin[i][1] - detectorCenter[1], moduleOrigin[i][2] - detectorCenter[2]};
		std::string detModuleAxis = "AXIS_D0M" + std::to_string(i);
		saveDouble(grp, detModuleAxis, 0.0, "degrees");
		dataset = H5Dopen2(grp , detModuleAxis.c_str(), H5P_DEFAULT);
		addStringAttribute(dataset, "depends_on", "AXIS_D0");
		addStringAttribute(dataset, "equipment", "detector");
		addStringAttribute(dataset, "equipment_component", "detector_module");
		addStringAttribute(dataset, "transformation_type", "rotation");
		addDoubleAttribute(dataset, "vector", mod_vector, 3);
		addDoubleAttribute(dataset, "offset", mod_offset, 3);
		addStringAttribute(dataset, "offset_units", "mm");
		H5Dclose(dataset);
	}

	H5Gclose(grp);

	for (int i = 0; i < NMODULES_GEOM; i++) {
		std::string moduleGroup = "/entry/instrument/detector/ARRAY_D0M" + std::to_string(i);
		grp = createGroup(master_file_id, moduleGroup.c_str() ,"NXdetector_module");
		int origin[2];

		if (NCH == 256*1024)
			origin[0] = startPos[2*i] / XPIXEL, origin[1] = startPos[2*i] % XPIXEL;
		else
			origin[0] = startPos[i] / XPIXEL, origin[1] = startPos[i] % XPIXEL;

		int size[2] = {514,1030};
		saveInt1D(grp, "data_origin", origin, "", 2);
		saveInt1D(grp, "data_size", size, "", 2);

		saveDouble(grp, "fast_pixel_direction", PIXEL_SIZE_IN_MM, "mm");

		double offset_fast[3] = {0,0,0};

		dataset = H5Dopen2(grp , "fast_pixel_direction", H5P_DEFAULT);
		addStringAttribute(dataset, "transformation_type","translation");
		addStringAttribute(dataset, "depends_on","/entry/instrument/" DETECTOR_NAME "/transformations/AXIS_D0M" + std::to_string(i));
		addDoubleAttribute(dataset, "offset", offset_fast, 3);
		if (fast == NULL) {
			double vector_fast[3] = {1,0,0};
			addDoubleAttribute(dataset, "vector", vector_fast, 3);
		} else {
			double vector_fast[3] = {fast[i][0],fast[i][1],fast[i][2]};
			addDoubleAttribute(dataset, "vector", vector_fast, 3);
		}

		H5Dclose(dataset);

		saveDouble(grp, "slow_pixel_direction", PIXEL_SIZE_IN_MM, "mm");

		double offset_slow[3] = {0,0,0};

		dataset = H5Dopen2(grp , "slow_pixel_direction", H5P_DEFAULT);
		addStringAttribute(dataset, "transformation_type","translation");
		addStringAttribute(dataset, "depends_on","/entry/instrument/" DETECTOR_NAME "/transformations/AXIS_D0M" + std::to_string(i));
		addDoubleAttribute(dataset, "offset", offset_slow, 3);
		if (slow == NULL) {
			double vector_slow[3] = {0,-1,0};
			addDoubleAttribute(dataset, "vector", vector_slow, 3);
		} else {
			double vector_slow[3] = {slow[i][0], slow[i][1], slow[i][2]};
			addDoubleAttribute(dataset, "vector", vector_slow, 3);
		}
		H5Dclose(dataset);

		H5Gclose(grp);
	}
}


void JungfrauHDF5Writer::SaveMetadata(const PXMetadata& in_metadata) {

	hid_t grp;
	grp = createGroup(master_file_id, "/entry", "NXentry");
	saveString(grp,"definition", "NXmx");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/data","NXdata");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument","NXinstrument");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/filter", "NXattenuator");
	if (px.transmission > 0.0) saveDouble(grp,"attenuator_transmission", px.transmission,"");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/" DETECTOR_NAME,"NXdetector_group");

	char *group_names = (char *) malloc(24 *sizeof(char *));
	strcpy(group_names, DETECTOR_NAME);
	strcpy(group_names+16,"detector");

	int32_t group_index[2] = {1,2};
	int32_t group_parent[2] = {-1, 1};
	int32_t group_type[2] = {1, 2};

	saveString1D(grp,"group_names", group_names, "", 2, 16);
	saveInt1D(grp, "group_parent", group_parent, "", 2);
	saveInt1D(grp, "group_index", group_index, "", 2);
	saveInt1D(grp, "group_type", group_type, "", 2);
	H5Gclose(grp);
	free(group_names);

	grp = createGroup(master_file_id, "/entry/instrument/beam","NXbeam");
	saveDouble(grp, "incident_wavelength",12400.0/in_metadata.photon_energy,"angstrom");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector","NXdetector");
	saveDouble(grp,"beam_center_x",in_metadata.beamx,"pixel");
	saveDouble(grp,"beam_center_y",in_metadata.beamy,"pixel");
	saveDouble(grp,"count_time",in_metadata.integration_time,"s");
	saveDouble(grp,"frame_time",in_metadata.frame_time,"s");
	saveDouble(grp,"detector_distance",in_metadata.detector_distance,"m");
	saveDouble(grp,"sensor_thickness",in_metadata.sensor_thickness,"m");
	saveDouble(grp,"x_pixel_size",PIXEL_SIZE_IN_MM/1000,"m");
	saveDouble(grp,"y_pixel_size",PIXEL_SIZE_IN_MM/1000,"m");
	saveString(grp,"sensor_material","Si");
	saveString(grp,"description","PSI Jungfrau 1M");
	if (mode == UINT16_MODE) {
		saveInt(grp,"saturation_value", UINT16_MAX-2);
	} else {
		saveInt(grp,"underload_value", INT32_MIN/2+2);
		saveInt(grp,"saturation_value", INT32_MAX/2-2);
	}
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/detectorSpecific","NXcollection");
	saveDouble(grp,"photon_energy",in_metadata.photon_energy,"eV");
	saveInt(grp,"nimages",nimages);
	saveInt(grp,"ntrigger",1);
	saveInt(grp,"x_pixels_in_detector",XPIXEL);
	saveInt(grp,"y_pixels_in_detector",YPIXEL);

	switch(compressionAlgorithm) {
	case COMPRESSION_BSHUF_LZ4:
		saveString(grp,"compression","bslz4");
		break;
	case COMPRESSION_GZIP:
		saveString(grp,"compression","gzip");
		break;
	case COMPRESSION_SZ:
		saveString(grp,"compression","sz");
		break;
	case COMPRESSION_LZ4:
		saveString(grp,"compression","lz4");
		break;
	case COMPRESSION_ZSTD:
		saveString(grp,"compression","zstd");
		break;
	case COMPRESSION_BSHUF_ZSTD:
		saveString(grp,"compression","bszstd");
		break;
	default:
		saveString(grp,"compression","");
		break;
	}

	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/detectorSpecific/detectorModule_000","");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/detectorSpecific/detectorModule_001","");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/detectorSpecific/detectorModule_002","");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/detectorSpecific/detectorModule_003","");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/detectorSpecific/detectorModule_004","");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/detectorSpecific/detectorModule_005","");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/detectorSpecific/detectorModule_006","");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/detectorSpecific/detectorModule_007","");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/instrument/detector/geometry","NXgeometry");
	H5Gclose(grp);

	//grp = createGroup(master_file_id, "/entry/instrument/detector/geometry/orientation","NXorientation");
	//double orient[] = {-1,0,0,0,-1,0};
	//saveDouble1D(grp, "value", orient, "", 6);
	//H5Gclose(grp);

	//grp = createGroup(master_file_id, "/entry/instrument/detector/geometry/translation","NXtranslation");
	//double trans[] = {in_metadata.beamx*in_metadata.x_pixel_size,in_metadata.beamy*in_metadata.y_pixel_size,-in_metadata.detector_distance};
	//saveDouble1D(grp, "distances", trans, "", 3);

	//H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/sample","NXsample");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/sample/beam","NXbeam");
	saveDouble(grp, "incident_wavelength",12400.0/in_metadata.photon_energy,"angstrom");
	H5Gclose(grp);

	grp = createGroup(master_file_id, "/entry/sample/goniometer","NXtransformations");


	if(px.stillImage) {
		SaveAngleContainer(grp,"omega", 0.0, 0.0, "degree");
		SaveAngleContainer(grp,"omega_end", 0.0, 0.0, "degree");
		saveDouble(grp,"omega_increment", 0.0,"degree");
		saveDouble(grp,"omega_range_average", 0.0,"degree");
		saveDouble(grp,"omega_range_total", 0.0,"degree");
		saveDouble(grp,"omega_start", 0.0,"degree");
	} else {
		double omega_range = in_metadata.omega_start + in_metadata.omega_increment * nimages;
		SaveAngleContainer(grp,"omega",in_metadata.omega_start, in_metadata.omega_increment, "degree");
		SaveAngleContainer(grp,"omega_end",in_metadata.omega_start+in_metadata.omega_increment, in_metadata.omega_increment, "degree");
		saveDouble(grp,"omega_increment", in_metadata.omega_increment,"degree");
		saveDouble(grp,"omega_range_average", in_metadata.omega_increment,"degree");
		saveDouble(grp,"omega_range_total", omega_range,"degree");
		saveDouble(grp,"omega_start", in_metadata.omega_start,"degree");
	}

	SaveAngleContainer(grp,"chi", 0.0, 0.0, "degree");
	SaveAngleContainer(grp,"chi_end", 0.0, 0.0, "degree");
	saveDouble(grp,"chi_increment", 0.0,"degree");
	saveDouble(grp,"chi_range_average", 0.0,"degree");
	saveDouble(grp,"chi_range_total", 0.0,"degree");
	saveDouble(grp,"chi_start", 0.0,"degree");

	SaveAngleContainer(grp,"kappa", 0.0, 0.0, "degree");
	SaveAngleContainer(grp,"kappa_end", 0.0, 0.0, "degree");
	saveDouble(grp,"kappa_increment", 0.0,"degree");
	saveDouble(grp,"kappa_range_average", 0.0,"degree");
	saveDouble(grp,"kappa_range_total", 0.0,"degree");
	saveDouble(grp,"kappa_start", 0.0,"degree");

	SaveAngleContainer(grp,"phi", 0.0, 0.0, "degree");
	SaveAngleContainer(grp,"phi_end", 0.0, 0.0, "degree");
	saveDouble(grp,"phi_increment", 0.0,"degree");
	saveDouble(grp,"phi_range_average", 0.0,"degree");
	saveDouble(grp,"phi_range_total", 0.0,"degree");
	saveDouble(grp,"phi_start", 0.0,"degree");

	H5Gclose(grp);

	WriteMetrology(in_metadata);
}

JungfrauHDF5Writer::~JungfrauHDF5Writer() {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	H5Fclose(master_file_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);
}


void JungfrauHDF5Writer::SaveData(char* data, int chunk_number, char *buffer) {

	int frames = images_per_file;
	if (chunk_number == nchunks - 1)
		frames = nimages-(nchunks-1)*images_per_file;

	// generate filename for data file
	char buff[12];
	snprintf(buff,12,"data_%06d",chunk_number+1);
	std::string filename = prefix+"_"+std::string(buff)+".h5";

	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	// Make link to particular chunk in the master dataset
	hid_t grp = H5Gopen2(master_file_id, "/entry/data", H5P_DEFAULT);
	createDataChunkLink(filename, grp, std::string(buff));
	H5Gclose(grp);

	// Add link to data in detector
	std::string entry_str = "/entry/data/" + std::string(buff);
	grp = H5Gopen2(master_file_id, "/entry/instrument/detector", H5P_DEFAULT);
	herr_t status = H5Lcreate_soft(entry_str.c_str(), grp, buff, H5P_DEFAULT, H5P_DEFAULT);
	H5Gclose(grp);

	// Create separate data file
	hid_t data_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

	// if (data_file_id < 0) return 1;

	grp = createGroup(data_file_id, "/entry", "NXentry");
	H5Gclose(grp);

	hid_t data_grp_id = createGroup(data_file_id, "/entry/data","NXdata");

	herr_t h5ret;

	// https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_crtdat.c
	hsize_t dims[3], chunk[3];

	dims[0] = frames;
	dims[1] = YPIXEL;
	dims[2] = XPIXEL;

	chunk[0] = 1;
	chunk[1] = YPIXEL;
	chunk[2] = XPIXEL;

	// Create the data space for the dataset.
	hid_t data_dataspace_id = H5Screate_simple(3, dims, NULL);

	hid_t data_dcpl_id = H5Pcreate (H5P_DATASET_CREATE);
	h5ret = H5Pset_chunk (data_dcpl_id, 3, chunk);

	// Set appropriate compression filter
	SetCompressionFilter(data_dcpl_id);

	// Create the dataset.
	hid_t data_dataset_id = H5Dcreate2(data_grp_id, "data", data_type, data_dataspace_id,
			H5P_DEFAULT, data_dcpl_id, H5P_DEFAULT);

	// Add attributes
	int tmp = chunk_number * images_per_file + 1;
	hid_t aid = H5Screate(H5S_SCALAR);
	hid_t attr = H5Acreate2(data_dataset_id, "image_nr_low", H5T_STD_I32LE, aid, H5P_DEFAULT, H5P_DEFAULT);
	h5ret = H5Awrite(attr, H5T_NATIVE_INT, &tmp);
	h5ret = H5Sclose(aid);
	h5ret = H5Aclose(attr);

	tmp = tmp+frames-1;
	aid = H5Screate(H5S_SCALAR);
	attr = H5Acreate2(data_dataset_id, "image_nr_high", H5T_STD_I32LE, aid, H5P_DEFAULT, H5P_DEFAULT);
	h5ret = H5Awrite(attr, H5T_NATIVE_INT, &tmp);
	h5ret = H5Sclose(aid);
	h5ret = H5Aclose(attr);

	if (use_direct_chunk_writer && (compressionAlgorithm != COMPRESSION_NONE)) {
		ret = pthread_mutex_unlock(&hdf5_semaphore);
		PTHREAD_ERROR(ret,pthread_mutex_unlock);

		for (unsigned int i = 0; i < frames; i ++) {
			size_t compressed_size = Compress(buffer, data+XPIXEL*YPIXEL*i*elem_size);

			ret = pthread_mutex_lock(&hdf5_semaphore);
			PTHREAD_ERROR(ret,pthread_mutex_lock);

			hsize_t offset[] = {i,0,0};
			h5ret = H5Dwrite_chunk(data_dataset_id, H5P_DEFAULT, 0, offset, compressed_size, buffer);
			HDF5_ERROR(h5ret,H5Dwrite_chunk);

			ret = pthread_mutex_unlock(&hdf5_semaphore);
			PTHREAD_ERROR(ret,pthread_mutex_unlock);
		}

		ret = pthread_mutex_lock(&hdf5_semaphore);
		PTHREAD_ERROR(ret,pthread_mutex_lock);

	} else {
		h5ret = H5Dwrite(data_dataset_id, data_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
		HDF5_ERROR(h5ret,H5Dwrite);
	}

	h5ret = H5Pclose (data_dcpl_id);

	/* End access to the dataset and release resources used by it. */
	h5ret = H5Dclose(data_dataset_id);

	/* Terminate access to the data space. */
	h5ret = H5Sclose(data_dataspace_id);

	h5ret = H5Gclose(data_grp_id);
	h5ret = H5Fclose(data_file_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);
}

void JungfrauHDF5Writer::SavePixelMask(uint32_t *pixelMask) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	hid_t group_id = H5Gopen2(master_file_id, "/entry/instrument/detector/detectorSpecific", H5P_DEFAULT);
	saveUInt2D(group_id, "pixel_mask", pixelMask, "", YPIXEL, XPIXEL);
	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

}

void JungfrauHDF5Writer::SaveDetectorSpecificData(std::string name, int *data, size_t size) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	hid_t group_id = H5Gopen2(master_file_id, "/entry/instrument/detector/detectorSpecific", H5P_DEFAULT);
	saveInt1D(group_id, name, data, "", size);
	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

}

void JungfrauHDF5Writer::SaveDetectorSpecificData(std::string name, int *data, size_t size, int module) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	char strbuf[100];
	sprintf(strbuf,"/entry/instrument/detector/detectorSpecific/detectorModule_%03d",module);

	hid_t group_id = H5Gopen2(master_file_id, strbuf, H5P_DEFAULT);
	saveInt1D(group_id, name, data, "", size);
	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

}

void JungfrauHDF5Writer::SaveDetectorSpecificData(std::string name, std::vector<double> &data) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	hid_t group_id = H5Gopen2(master_file_id, "/entry/instrument/detector/detectorSpecific", H5P_DEFAULT);
	saveDouble1D(group_id, name, &data.front(), "", data.size());
	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

}

void JungfrauHDF5Writer::SaveDetectorSpecificData(std::string name, std::vector<double> &data, int module) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	char strbuf[100];
	sprintf(strbuf,"/entry/instrument/detector/detectorSpecific/detectorModule_%03d",module);

	hid_t group_id = H5Gopen2(master_file_id, strbuf, H5P_DEFAULT);
	saveDouble1D(group_id, name, &data.front(), "", data.size());

	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

}

void JungfrauHDF5Writer::SaveDetectorSpecificData(std::string name, std::vector<int> &data) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	hid_t group_id = H5Gopen2(master_file_id, "/entry/instrument/detector/detectorSpecific", H5P_DEFAULT);
	saveInt1D(group_id, name, &data.front(), "", data.size());
	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);
}

void JungfrauHDF5Writer::SaveDetectorSpecificData(std::string name, std::vector<int> &data, int module) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	char strbuf[100];
	sprintf(strbuf,"/entry/instrument/detector/detectorSpecific/detectorModule_%03d",module);

	hid_t group_id = H5Gopen2(master_file_id, strbuf, H5P_DEFAULT);
	saveInt1D(group_id, name, &data.front(), "", data.size());

	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);
}

void JungfrauHDF5Writer::SaveDetectorSpecificData(std::string name, std::string data) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	hid_t group_id = H5Gopen2(master_file_id, "/entry/instrument/detector/detectorSpecific", H5P_DEFAULT);
	saveString(group_id, name, data);
	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

}

void JungfrauHDF5Writer::SaveDetectorSpecificData(std::string name, std::string data, int module) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	char strbuf[100];
	sprintf(strbuf,"/entry/instrument/detector/detectorSpecific/detectorModule_%03d",module);

	hid_t group_id = H5Gopen2(master_file_id, strbuf, H5P_DEFAULT);
	saveString(group_id, name, data);
	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

}

void JungfrauHDF5Writer::SaveDetectorSpecificData(std::string name, double data, int module) {
	int ret = pthread_mutex_lock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_lock);

	char strbuf[100];
	sprintf(strbuf,"/entry/instrument/detector/detectorSpecific/detectorModule_%03d",module);

	hid_t group_id = H5Gopen2(master_file_id, strbuf, H5P_DEFAULT);
	saveDouble(group_id, name, data);
	H5Gclose(group_id);

	ret = pthread_mutex_unlock(&hdf5_semaphore);
	PTHREAD_ERROR(ret,pthread_mutex_unlock);

}

int JungfrauHDF5Writer::SaveAngleContainer(hid_t location, std::string name, double start, double increment, std::string units) {
	double val[nimages];
	for (int i = 0; i < nimages; i++) {
		val[i] = start + i * increment;
	}
	return saveDouble1D(location, name, val, units, nimages);
}

size_t JungfrauHDF5Writer::GetOutputBufferMaxSize() {
	return output_size;
}
