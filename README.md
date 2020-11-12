# CalibrationCpp
Calibration of snowflakes parameters using Bayesian approach (cpp version)

# Installation

Install static libraries
 cd static_lib/galib274
 # Set the variable "DESTDIR"="absolute_path_to_static_lib" in makevars file
 echo $PWD
 make install
 (make clean) ~ to remove build files (free space)

 cd ../..

 cd static_lib/libconfig-1.5
 echo $PWD
 ./configure --prefix="absolute_path_to_static_lib"
 make install
 (make clean) ~ to remove build files (free space)

 cd ../..


Compile main program
 make
