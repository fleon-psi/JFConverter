/*
 * This file is an example of an HDF5 filter plugin.
 * The plugin can be used with the HDF5 library vesrion 1.8.11+ to read
 * HDF5 datasets compressed with lz4.
 */

#include <sys/types.h>
#include <stdlib.h>
#include <string.h>
#include <assert.h>
#include <stdio.h>
#if defined(_WIN32)
#include <Winsock2.h>
#endif
#include "H5PLextern.h"
#include "lz4.h"

static size_t H5Z_filter_lz4(unsigned int flags, size_t cd_nelmts,
        const unsigned int cd_values[], size_t nbytes,
        size_t *buf_size, void **buf);

#define H5Z_FILTER_LZ4 32004

#define htonll(x) ( ( (uint64_t)(htonl( (uint32_t)((x << 32) >> 32)))<< 32) | htonl( ((uint32_t)(x >> 32)) ))
#define ntohll(x) htonll(x)

#define htobe16t(x) htons(x)
#define htobe32t(x) htonl(x)
#define htobe64t(x) htonll(x)
#define be16toht(x) ntohs(x)
#define be32toht(x) ntohl(x)
#define be64toht(x) ntohll(x)


#define DEFAULT_ACCEL 1; /* 1GB. LZ4 needs blocks < 1.9GB. */

const H5Z_class2_t H5Z_LZ4[1] = {{
        H5Z_CLASS_T_VERS,       /* H5Z_class_t version */
        (H5Z_filter_t)H5Z_FILTER_LZ4,         /* Filter id number             */
        1,              /* encoder_present flag (set to true) */
        1,              /* decoder_present flag (set to true) */
        "HDF5 lz4 filter; see  http://github.com/RGLab/h5lz",
        /* Filter name for debugging    */
        NULL,                       /* The "can apply" callback     */
        NULL,                       /* The "set local" callback     */
        (H5Z_func_t)H5Z_filter_lz4,         /* The actual filter function   */
}};

//H5PL_type_t   H5PLget_plugin_type(void) {return H5PL_TYPE_FILTER;}
//const void *H5PLget_plugin_info(void) {return H5Z_LZ4;}

static size_t H5Z_filter_lz4(unsigned int flags, size_t cd_nelmts,
        const unsigned int cd_values[], size_t nbytes,
        size_t *buf_size, void **buf)
{
    void * outBuf = NULL;
    size_t ret_value;

    if (flags & H5Z_FLAG_REVERSE)
    {
        uint32_t *i32Buf;
        uint32_t compressedSize;
        char *roBuf;   /* pointer to current write position */
        uint64_t decompSize;
        const char* rpos = (char*)*buf; /* pointer to current read position */
        const uint64_t * const i64Buf = (uint64_t *) rpos;
        const uint64_t origSize = (uint64_t)(be64toht(*i64Buf));/* is saved in be format */
        rpos += 8; /* advance the pointer */

        i32Buf = (uint32_t*)rpos;
        compressedSize = (uint32_t)(be32toht(*i32Buf));/// is saved in be format

        rpos += 4;

        if (NULL==(outBuf = malloc(origSize)))
        {
            printf("cannot malloc\n");
            goto error;
        }
        roBuf = (char*)outBuf;   /* pointer to current write position */
        decompSize     = 0;
        /// start with the first block ///
        uint32_t nDest = LZ4_decompress_fast(rpos, roBuf, origSize);
		if(nDest != compressedSize)
		{
			 printf("decompressed size not the same: %d, != %d\n", nDest, compressedSize);
			 goto error;
		}
        free(*buf);
        *buf = outBuf;
        outBuf = NULL;
        ret_value = (size_t)origSize;  // should always work, as orig_size cannot be > 2GB (sizeof(size_t) < 4GB)
    }
    else /* forward filter */
    {
    	size_t acceleration;
        size_t outSize; /* size of the output buffer. Header size (12 bytes) is included */
        uint64_t *i64Buf;
        uint32_t *i32Buf;
        char *rpos;      /* pointer to current read position */
        char *roBuf;    /* pointer to current write position */

        if (nbytes > INT32_MAX)
        {
            /* can only compress chunks up to 2GB */
            goto error;
        }

        if(cd_nelmts > 0 && cd_values[0] > 0)
        {
        	acceleration = cd_values[0];
        }
        else
        {
        	acceleration = DEFAULT_ACCEL;
        }
        uint32_t nMaxSize = LZ4_COMPRESSBOUND(nbytes);

        if (NULL==(outBuf = malloc(nMaxSize)))
        {
            goto error;
        }

        rpos  = (char*)*buf;      /* pointer to current read position */
        roBuf = (char*)outBuf;    /* pointer to current write position */
        /* header */
        i64Buf = (uint64_t *) (roBuf);
        i64Buf[0] = htobe64t((uint64_t)nbytes); /* Store decompressed size in be format */
        roBuf += 8;

        uint32_t nCompressed = LZ4_compress_fast(rpos, roBuf+4, nbytes, nMaxSize, acceleration);
        if(!nCompressed)
        	goto error;
        outSize =  nCompressed + 12; /* size of the output buffer. Header size (12 bytes) is included */
        i32Buf =  (uint32_t *) (roBuf);
        i32Buf[0] = htobe32t((uint32_t)nCompressed); /* Store the compressed size in be format */

        free(*buf);
        *buf = outBuf;
        *buf_size = outSize;
        outBuf = NULL;
        ret_value = outSize;

    }

    if(outBuf)
        free(outBuf);
    return ret_value;


    error:
    if(outBuf)
        free(outBuf);
    outBuf = NULL;
    return 0;
}
