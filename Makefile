HOSTNAME := $(shell hostname)
UNAME_P := $(shell uname -p)

ifndef HDF5_PATH
    HDF5_PATH=/opt/hdf5
endif

SZ_PATH=SZ/sz
LZ_PATH=lz4
GZIP_LIB=-lz

CPPFLAGS= -I$(HDF5_PATH)/include -I$(LZ_PATH) -I$(SZ_PATH)/include -Ijson/single_include -Izstd/lib
JF_LDLIBS=$(HDF5_PATH)/lib/libhdf5.a zstd/lib/libzstd.a $(GZIP_LIB) -ldl

SRCS= parseMetadata.o JungfrauFile.o JungfrauPixelMask.o JungfrauHDF5Writer.o JungfrauHDF5Writer_h5functions.o\
      JungfrauHDF5Writer_compress.o JungfrauConverter.o sharedVariables.o WriterThread.o ConverterThread.o main.o \
      bitshuffle/bitshuffle.o bitshuffle/bitshuffle_core.o bitshuffle/iochain.o bitshuffle/bshuf_h5filter.o\
      lz4/lz4.o lz4_plugin/H5Zlz4.o zstd_plugin/zstd_h5plugin.o

ifeq ($(UNAME_P), ppc64le)
    CC=xlc_r
    CXX=xlC_r
    LDFLAGS= -g -O5 -qarch=pwr9 -qsimd -qhot
    CXXFLAGS= -g -std=gnu++11 -O5 -qarch=pwr9 -qsimd -qhot
    CFLAGS= -g -std=c99 -O5 -qarch=pwr9 -qsimd -qhot
endif

ifeq ($(UNAME_P), x86_64)
	ifeq ($(HOSTNAME), mx-jungfrau-1.psi.ch)
	    CXX=icpc
	endif

	ifeq ($(CXX),g++)
	    CC=gcc
	    CXX=g++
	    LDFLAGS=-Ofast -pthread -march=native
	    CXXFLAGS=-std=c++11 -Ofast -pthread -march=native
	    CFLAGS=-std=c99 -Ofast -pthread -march=native
	else
            JF_LDLIBS += $(IPPROOT)/lib/intel64/libippdc.a $(IPPROOT)/lib/intel64/libipps.a $(IPPROOT)/lib/intel64/libippcore.a
            CPPFLAGS+=-DWITH_IPP
	    CC=icc
	    CXX=icpc
	    LDFLAGS= -Ofast -ipo -pthread -fno-alias -xHost
	    CXXFLAGS= -std=c++11 -Ofast -ipo -pthread -xHost -qopt-zmm-usage=high -fno-alias
	    CFLAGS= -std=c99 -Ofast -ipo -pthread -xHost -qopt-zmm-usage=high -fno-alias

#	    LDFLAGS= -Ofast -ipo -g -shared-intel -debug inline-debug-info -pthread -qopt-report=5 -qopt-report-phase=pgo -fno-alias
#	    CXXFLAGS= -std=c++11 -Ofast -ipo -g -shared-intel -debug inline-debug-info -pthread -xHost -qopt-zmm-usage=high -fno-alias -prof-gen-sampling
#	    CFLAGS= -std=c99 -Ofast -ipo -g -shared-intel -debug inline-debug-info -pthread -xHost -qopt-zmm-usage=high -fno-alias -prof-gen-sampling

	endif
endif

ifeq ($(HOSTNAME), mx-ac922-test.psi.ch)
   CPPFLAGS += -DAC922   
endif

ifeq ($(HOSTNAME), mx-jungfrau-1.psi.ch)
   CPPFLAGS += -DHPEDL580
endif

all: JFConverter

JFConverter: $(SRCS) jungfrau.h
	$(CXX) $(LDFLAGS) $(SRCS) -o JFConverter $(JF_LDLIBS)

clean:
	rm -f *.o bitshuffle/*.o lz4/*.o lz4_plugin/*.o JFConverter

install:
	cp JFConverter /usr/local/bin
