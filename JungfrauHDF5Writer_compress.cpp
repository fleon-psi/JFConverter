#include <cstring>

#include "jungfrau.h"

#include "bitshuffle.h"
#include "bitshuffle_core.h"

#include "lz4.h"

#include "zlib.h"
#include "zstd.h"

//#include "sz.h"

#define DEFLATE_SIZE_ADJUST(s) (ceil(((double)(s))*1.001)+12)

extern "C" {
// prototypes for bitshuffle internal functions to write big-endian integers into a buffer
// necessary for bshuf header in H5
void bshuf_write_uint64_BE(void* buf, uint64_t num);
void bshuf_write_uint32_BE(void* buf, uint32_t num);

extern H5Z_class_t bshuf_H5Filter[1];
extern H5Z_class_t H5Z_LZ4[1];
extern H5Z_class_t zstd_H5Filter[1];
}

void JungfrauHDF5Writer::InitCompression() {
	H5Zregister(&bshuf_H5Filter);
	H5Zregister(&H5Z_LZ4);
	H5Zregister(&zstd_H5Filter);
}

size_t JungfrauHDF5Writer::CompressionOutputMaxSize() {
	switch(compressionAlgorithm) {
	case COMPRESSION_BSHUF_LZ4:
		return bshuf_compress_lz4_bound(XPIXEL*YPIXEL,elem_size,bshuf_block_size) + 12;
	case COMPRESSION_BSHUF_ZSTD:
		return bshuf_compress_zstd_bound(XPIXEL*YPIXEL,elem_size,bshuf_block_size) + 12;
	case COMPRESSION_GZIP:
		return DEFLATE_SIZE_ADJUST(XPIXEL*YPIXEL*elem_size);
	case COMPRESSION_LZ4:
		return LZ4_COMPRESSBOUND(XPIXEL*YPIXEL*elem_size) + 16;
	case COMPRESSION_ZSTD:
		return ZSTD_compressBound(XPIXEL*YPIXEL*elem_size);
	default:
		return XPIXEL*YPIXEL*elem_size;
	}
}

size_t JungfrauHDF5Writer::Compress(char *output,  char *input) {
	int gzip_ret;
	uLongf z_dst_nbytes;
	size_t res;
	size_t block_size;
	size_t input_size = XPIXEL*YPIXEL*elem_size;

	switch(compressionAlgorithm) {
	case COMPRESSION_BSHUF_LZ4:
		if (bshuf_block_size == 0)
			// Default block size determined by bitshuffle internals
			block_size = bshuf_default_block_size(elem_size);
		else block_size = bshuf_block_size;
		// write header - taken from bshuf_h5filter.c
		// uncompressed size
		bshuf_write_uint64_BE(output, input_size);

		// block_size * element_size
		bshuf_write_uint32_BE(output + 8, block_size * elem_size );

		return bshuf_compress_lz4(input, output+12, XPIXEL*YPIXEL,elem_size,block_size) + 12;

	case COMPRESSION_BSHUF_ZSTD:
		if (bshuf_block_size == 0)
			// Default block size determined by bitshuffle internals
			block_size = bshuf_default_block_size(elem_size);
		else block_size = bshuf_block_size;

		// write header - taken from bshuf_h5filter.c
		// uncompressed size
		bshuf_write_uint64_BE(output, input_size);

		// block_size * element_size
		bshuf_write_uint32_BE(output + 8, block_size * elem_size );

		return bshuf_compress_zstd(input, output+12, XPIXEL*YPIXEL,elem_size,block_size) + 12;
	case COMPRESSION_GZIP:
		z_dst_nbytes = (uLongf) output_size;

		gzip_ret = compress2((Bytef *) output, &z_dst_nbytes, (const Bytef *) input, (uLong) input_size, gzip_level);

		/* Check for various zlib errors */
		if(gzip_ret == Z_BUF_ERROR)
			std::cerr << "gzip overflow" << std::endl;
		else if(gzip_ret == Z_MEM_ERROR)
			std::cerr << "gzip mem error" << std::endl;
		else if(gzip_ret != Z_OK)
			std::cerr << "gzip other error" << std::endl;
		return z_dst_nbytes;
	case COMPRESSION_LZ4:
		// Currently whole image is saved as a single LZ4 block
		// Write uncompressed size
		bshuf_write_uint64_BE(output, input_size);

		// Write block size
		bshuf_write_uint32_BE((char*) output + 8, input_size);

		// Compress
		res = LZ4_compress_default(input, output + 16, input_size, output_size - 16);

		// Write size after compression
		bshuf_write_uint32_BE(output + 12, res);
		return res + 16;
	case COMPRESSION_ZSTD:
		return ZSTD_compress(output, output_size, input, input_size, ZSTD_CLEVEL_DEFAULT);
	default:
		memcpy(output, input,XPIXEL*YPIXEL*elem_size);
		return XPIXEL*YPIXEL*elem_size;
	}
}

void JungfrauHDF5Writer::SetCompressionFilter(hid_t data_dcpl_id) {
	herr_t h5ret;
	switch(compressionAlgorithm) {
	case COMPRESSION_BSHUF_LZ4:
	{
		unsigned int params[] = {bshuf_block_size, BSHUF_H5_COMPRESS_LZ4};
		h5ret = H5Pset_filter(data_dcpl_id, (H5Z_filter_t)BSHUF_H5FILTER, H5Z_FLAG_MANDATORY, (size_t)2, params);
		HDF5_ERROR(h5ret,H5Pset_filter);
		break;
	}
	case COMPRESSION_BSHUF_ZSTD:
	{
		unsigned int params[] = {bshuf_block_size, BSHUF_H5_COMPRESS_ZSTD};
		h5ret = H5Pset_filter(data_dcpl_id, (H5Z_filter_t)BSHUF_H5FILTER, H5Z_FLAG_MANDATORY, (size_t)2, params);
		HDF5_ERROR(h5ret,H5Pset_filter);
		break;
	}
	case COMPRESSION_GZIP:
		h5ret = H5Pset_deflate(data_dcpl_id, gzip_level);
		HDF5_ERROR(h5ret,H5Pset_deflate);
		break;
	case COMPRESSION_LZ4:
	{
		const unsigned int lz4_values[1] = {3}; // taken from LZ4 HDF5 Plugin example
		h5ret = H5Pset_filter(data_dcpl_id, LZ4_H5FILTER, H5Z_FLAG_MANDATORY, (size_t)1, lz4_values);
		HDF5_ERROR(h5ret,H5Pset_filter);
		break;
	}
	case COMPRESSION_SZ:
		std::cout << "SZ not yet implemented, no compression applied" << std::endl;
		break;
	case COMPRESSION_ZSTD:
		h5ret = H5Pset_filter(data_dcpl_id, ZSTD_H5FILTER, H5Z_FLAG_MANDATORY, (size_t)0, NULL);
		break;
	case COMPRESSION_NONE:
		break;
	}
}
