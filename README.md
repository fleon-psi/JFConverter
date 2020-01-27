# JFConverter
Conversion from JUNGFRAU detector format to photons.

Currently converts data collected for synhcrotron MX experiments with JF 1M (new firmware) and JF 4M (old firmware), using SLS Detector Package and with PSI MX control scripts.

# Requirements
* Recent HDF5 library (e.g. 1.10.5)
* Zlib - gzip library - part of standard Linux distribution
* Compression libraries (Zstandard, bitshuffle, lz4) and HDF5 plugins - available as submodules
* JSON C++ library - available as submodule

# Compilation

Might require adjusting paths in the Makefile and manually compiling submodules.

Path to gain files needs also to be manually adjusted within jungfrau.h file.

```
git clone https://github.com/fleon-psi/JFConverter
git submodule update --init
make
```

# Test data

Available from:
https://dx.doi.org/10.16907/808de0df-a9d3-4698-8e9f-d6e091516650

Test data can be processed with the following command:
```
JFConverter -j NA102I_Lyso5_12p4_100dps_tr1p0_720_1.json  -r -o /mnt/ssd/output/ -f2045 -l2500 
