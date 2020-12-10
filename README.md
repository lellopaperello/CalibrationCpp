# CalibrationCpp
Calibration of snowflakes parameters using Bayesian approach (cpp version)

# Installation

Install static libraries ("~" commands are suggestions)

 cd static_lib/galib247
 # Set the variable "DESTDIR"="absolute_path_to_static_lib" in makevars file
 ~ echo $(dirname "$PWD")
 ~ nano makevars
 make install
 (make clean) ~ to remove build files (free space)

 cd ../..

 cd static_lib/libconfig-1.5
 ~ echo $(dirname "$PWD")
 ./configure --prefix="absolute_path_to_static_lib"
 make install
 (make clean) ~ to remove build files (free space)

 cd ../..


Compile main program
 make
