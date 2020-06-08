# MuLTI
Multimodal Layered Transdimensional Inversion of Seismic Dispersion Curves with Depth Constraints.

The MuLTI algorithm was developed in Matlab version 2017a, therefore all Matlab codes provided will work with this version or subsequent Matlab versions. This code can be run on Windows or Linux based platforms. This github repository includes: the main MuLTI matlab codes; gpdc mex files used to enable the original gpdc code from Wathelet (2005), written in C++, to be called in Matlab within the MuLTI script for both Linux and Windows based platforms; all fuctions called within MuLTI; tools to display all results output from MuLTI; example datasets.

## gpdc forward model
To run the gpdc forward modelling calculation in Matlab, within the MuLTI algorithm, the gpdc Geospy software version 2.10 must be installed in the Matlab bin directory. This can be obtained from http://www.geopsy.org/. The mex interface has been mostly tested on CentOS 7, and Windows 7, with Matlab 2017a and geopsy 2.10.0. Further details and instructions on how to compile and use the gpdc mex file can be found in the corresponding folders and at the link: https://github.com/cemac/MEX-gpdc. If working on a Linux based platform the gpdc mex file created for Linux specifically must be compiled in the Matlab working directory and then the gpdc function can be called in Matlab. If working on a Windows platform the gpdc mex files created for windows specifically does not need compiled and can just be called in the matlab working directory. 

## Functions
To run the MuLTI Matlab code ('MuLTI.m') the corresponding platform based gpdc mex file must be in the active Matlab working directory along with the Matlab functions ‘thicknesses_and_priors.m’ and ‘whichnuclei.m’.

## Input data (MuLTI)
The input dispersion curve data files are .mat files with dispersion curve picks saved as a column vector variable called “data” with column 1 frequency in Hertz and column 2 phase velocity in m/s. The fitting error is determined from the half width of the dispersion curve image, this is saved as a column vector variable called “half_width” with column 1 frequency and column 2 the half width of the dispersion curve in m/s.

# MuLTI III
MuLTI III has been developed further to address key limitations in the original MuLTI code, specifically that Vp and density must be fixed. MuLTI III overcomes this limitation by allowing both Vp and density to vary, together with estimates of their curves uncertainty. This tool is useful when the Vp and density structures of the subsurface are known, for example, from seismic refraction investigations and borehole measurements.

## Functions
MuLTI III setup procedures are very similar to the original MuLTI code, described above. However, the function “thickneses_and_priors_III.m” is now needed for MuLTI III ('MuLTI_III.m') instead of the original “thickneses_and_priors.m”. The original function ‘whichnuclei.m’ is still needed also.

## Input data (MuLTI III)
The input dispersion curve and fitting error data format are the same as for the original MuLTI code.
The input Vp profiles are saved as .mat files containing a column vector variable called “vpdata” with column 1 Vp depths in meters, column 2 mean Vp value in m/s and column 3 one standard deviation (estimated error of Vp).
The input density profiles are saved as .mat files containing a column vector variable called “dendata” with column 1 density depths in meters, column 2 mean density value in g/cm3 and column 3 one standard deviation (estimated error of density).


