#include <iostream>
#include <pthread.h>

#include "jungfrau.h"

#define H5Z_FILTER_LZ4        32004

int addStringAttribute(hid_t location, std::string name, std::string val) {
    /* https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_attribute.c */
    hid_t aid = H5Screate(H5S_SCALAR);
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, val.length());
    H5Tset_strpad(atype,H5T_STR_NULLTERM);
    hid_t attr = H5Acreate2(location, name.c_str(), atype, aid, H5P_DEFAULT, H5P_DEFAULT);
    herr_t ret = H5Awrite(attr, atype, val.c_str());
    ret = H5Sclose(aid);
    ret = H5Tclose(atype);
    ret = H5Aclose(attr);
    
    return 0;
    
}

int addDoubleAttribute(hid_t location, std::string name, const double *val, int dim) {
    // https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_crtdat.c
    hsize_t dims[1];
    dims[0] = dim;

    /* Create the data space for the dataset. */
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

    hid_t attr = H5Acreate2(location, name.c_str(), H5T_IEEE_F64LE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT);

    herr_t ret = H5Awrite(attr, H5T_NATIVE_DOUBLE, val);
    HDF5_ERROR(ret, H5Awrite);

    ret = H5Sclose(dataspace_id);
    ret = H5Aclose(attr);

    return 0;
}

hid_t createGroup(hid_t master_file_id, std::string group, std::string nxattr) {
    hid_t group_id = H5Gcreate(master_file_id, group.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    if (nxattr != "") addStringAttribute(group_id, "NX_class", nxattr);
    return group_id;
}

int saveString1D(hid_t location, std::string name, char *val, std::string units, int dim, int len) {
    herr_t status;

    // https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_crtdat.c
    hsize_t dims[1];
    dims[0] = dim;

    /* Create the data space for the dataset. */
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);

    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, len);
    H5Tset_strpad(atype,H5T_STR_NULLTERM);
    /* Create the dataset. */
    hid_t dataset_id = H5Dcreate2(location, name.c_str(), atype, dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    /* Write the dataset. */
    status = H5Dwrite(dataset_id, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      val);
    if (units != "") addStringAttribute(dataset_id, "units", units);

    /* End access to the dataset and release resources used by it. */
    status = H5Dclose(dataset_id);

    /* Terminate access to the data space. */
    status = H5Sclose(dataspace_id);
    return 0;
}

int saveDouble1D(hid_t location, std::string name, const double *val, std::string units, int dim) {
    herr_t status;
    
    // https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_crtdat.c
    hsize_t dims[1];
    dims[0] = dim;
    
    /* Create the data space for the dataset. */
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    /* Create the dataset. */
    hid_t dataset_id = H5Dcreate2(location, name.c_str(), H5T_IEEE_F32LE, dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    /* Write the dataset. */
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      val);
    
    if (units != "") addStringAttribute(dataset_id, "units", units);
    
    /* End access to the dataset and release resources used by it. */
    status = H5Dclose(dataset_id);
    
    /* Terminate access to the data space. */
    status = H5Sclose(dataspace_id);
    
    return 0;
}

int saveUInt2D(hid_t location, std::string name, const uint32_t *val, std::string units, int dim1, int dim2) {
    herr_t status;
    
    // https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_crtdat.c
    hsize_t dims[2];
    dims[0] = dim1;
    dims[1] = dim2;
    
    /* Create the data space for the dataset. */
    hid_t dataspace_id = H5Screate_simple(2, dims, NULL);
    
    /* Create the dataset. */
    hid_t dataset_id = H5Dcreate2(location, name.c_str(), H5T_STD_U32LE, dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    /* Write the dataset. */
    status = H5Dwrite(dataset_id, H5T_NATIVE_UINT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      val);
    
    if (units != "") addStringAttribute(dataset_id, "units", units);
    
    /* End access to the dataset and release resources used by it. */
    status = H5Dclose(dataset_id);
    
    /* Terminate access to the data space. */
    status = H5Sclose(dataspace_id);
    
    return 0;
}

int saveInt1D(hid_t location, std::string name, const int *val, std::string units, int dim) {
    herr_t status;
    
    // https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_crtdat.c
    hsize_t dims[1];
    dims[0] = dim;
    
    /* Create the data space for the dataset. */
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    /* Create the dataset. */
    hid_t dataset_id = H5Dcreate2(location, name.c_str(), H5T_STD_I32LE, dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    /* Write the dataset. */
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      val);
    
    if (units != "") addStringAttribute(dataset_id, "units", units);
    
    /* End access to the dataset and release resources used by it. */
    status = H5Dclose(dataset_id);
    
    /* Terminate access to the data space. */
    status = H5Sclose(dataspace_id);
    
    return 0;
}

int saveString(hid_t location, std::string name, std::string val, std::string units) {
    herr_t status;
    
    // https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_crtdat.c
    hsize_t dims[1];
    dims[0] = 1;
    
    /* Create the data space for the dataset. */
    hid_t dataspace_id = H5Screate_simple(1, dims, NULL);
    
    hid_t atype = H5Tcopy(H5T_C_S1);
    H5Tset_size(atype, val.length()+1);
    H5Tset_strpad(atype,H5T_STR_NULLTERM);
    /* Create the dataset. */
    hid_t dataset_id = H5Dcreate2(location, name.c_str(), atype, dataspace_id,
                                  H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    
    /* Write the dataset. */
    status = H5Dwrite(dataset_id, atype, H5S_ALL, H5S_ALL, H5P_DEFAULT,
                      val.c_str());
    if (units != "") addStringAttribute(dataset_id, "units", units);
    
    /* End access to the dataset and release resources used by it. */
    status = H5Dclose(dataset_id);
    
    /* Terminate access to the data space. */
    status = H5Sclose(dataspace_id);
    return 0;
}


int saveDouble(hid_t location, std::string name, double val, std::string units) {
    double tmp = val;
    return saveDouble1D(location, name, &tmp, units, 1);
}

int saveInt(hid_t location, std::string name, int val, std::string units) {
    int tmp = val;
    return saveInt1D(location, name, &tmp, units, 1);
}



// X and Y are swapped in Dectris HDF5 !!!
int createDataChunk(std::string filename, int frame_low, int *data, int frames) {
    /* Create file */
    hid_t data_file_id = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    
    if (data_file_id < 0) return 1;
    
    hid_t grp;
    grp = createGroup(data_file_id, "/entry", "NXentry");
    H5Gclose(grp);
    
    grp = createGroup(data_file_id, "/entry/data","NXdata");

    herr_t status;
    
    // https://support.hdfgroup.org/ftp/HDF5/current/src/unpacked/examples/h5_crtdat.c
    hsize_t dims[3], chunk[3];
    
    dims[0] = frames;
    dims[1] = YPIXEL;
    dims[2] = XPIXEL;
    
    chunk[0] = 1;
    chunk[1] = YPIXEL;
    chunk[2] = XPIXEL;
    
    /* Create the data space for the dataset. */
    hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
    
    hid_t dcpl = H5Pcreate (H5P_DATASET_CREATE);
    status = H5Pset_chunk (dcpl, 3, chunk);
    unsigned int params[] = {0,2};

    status = H5Pset_filter(dcpl,(H5Z_filter_t)32008,H5Z_FLAG_MANDATORY,(size_t)2,params);
    
    /* Create the dataset. */
    hid_t dataset_id = H5Dcreate2(grp, "data", H5T_STD_I32LE, dataspace_id,
                                  H5P_DEFAULT, dcpl, H5P_DEFAULT);
    
    /* Write dataset. */
    H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data);
    
    /* Add attributes */
    int tmp = frame_low+1;
    hid_t aid = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate2(dataset_id, "image_nr_low", H5T_STD_I32LE, aid, H5P_DEFAULT, H5P_DEFAULT);
    herr_t ret = H5Awrite(attr, H5T_NATIVE_INT, &tmp);
    
    ret = H5Sclose(aid);
    ret = H5Aclose(attr);
    
    tmp = frame_low+frames;
    aid = H5Screate(H5S_SCALAR);
    attr = H5Acreate2(dataset_id, "image_nr_high", H5T_STD_I32LE, aid, H5P_DEFAULT, H5P_DEFAULT);
    ret = H5Awrite(attr, H5T_NATIVE_INT, &tmp);
    
    ret = H5Sclose(aid);
    ret = H5Aclose(attr);
    
    status = H5Pclose (dcpl);
    
    /* End access to the dataset and release resources used by it. */
    status = H5Dclose(dataset_id);
    
    /* Terminate access to the data space. */
    status = H5Sclose(dataspace_id);
    std::cout << "Written " << frames << " frames ("<<frame_low <<"-"<<frames+frame_low-1 <<") to "<< filename <<"." << std::endl;

    H5Gclose(grp);
    H5Fclose(data_file_id);
    return 0;
}

int createDataChunkLink(std::string filename, hid_t location, std::string name) {
    std::string pure_filename;
    int pos = filename.rfind("/");
    if (pos != std::string::npos) {
        herr_t status = H5Lcreate_external(filename.substr(pos+1).c_str(),"/entry/data/data",location,name.c_str(), H5P_DEFAULT, H5P_DEFAULT);
    } else {
        herr_t status = H5Lcreate_external(filename.c_str(),"/entry/data/data",location,name.c_str(), H5P_DEFAULT, H5P_DEFAULT);

    }
    return 0;
}


