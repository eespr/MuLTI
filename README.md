# MuLTI
Multimodal Layered Transdimensional Inversion of Seismic Dispersion Curves with Depth Constraints.

The MuLTI algorithm was developed in Matlab version 2017a, therefore all Matlab codes supplied with this paper will work with this version or subsequent Matlab versions. This code can be run on Windows or Linux based platforms. This github repository includes: gpdc mex files used to enable the original gpdc code from Wathelet (2005), written in C++, to be called in Matlab within the MuLTI script for both Linux and Windows based platforms; all MuLTI Matlab codes used to run synthetic and real data examples, along with all input dispersion curve data files and all results from these data examples.

To run the gpdc forward modelling calculation in Matlab, within the MuLTI algorithm, the gpdc Geospy software version 2.10 must be installed in the Matlab bin directory. This can be obtained from http://www.geopsy.org/. The mex interface has been mostly tested on CentOS 7, and Windows 7, with Matlab 2017a and geopsy 2.10.0. Further details and instructions on how to compile and use the gpdc mex file can be found in the corresponding folders. If working on a Linux based platform the gpdc mex file created for Linux specifically must be compiled in the Matlab working directory and then the gpdc function can be called in Matlab. If working on a Windows platform the gpdc mex files created for windows specifically does not need compiled and can just be called in the matlab working directory. 

To run the MuLTI Matlab codes, .m files, the corresponding platform based gpdc mex file must be in the active Matlab working directory along with the Matlab functions ‘thicknesses_and_priors’ and ‘whichnuclei’ and the input dispersion curve data files.

The input dispersion curve data files are .mat files with dispersion curve picks saved as a column vector variable called “data” with column 1 being frequency in Hertz and column 2 being phase velocity in m/s. The fitting error is determined from the half width of the dispersion curve image, this is saved as a column vector variable called “half_width” with column 1 being frequency and column 2 being the half width of the dispersion curve in m/s.

DOI 10.5281/zenodo.1489959
