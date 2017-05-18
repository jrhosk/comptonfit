# comptonfit

# Required libraries
In order to run this software you will need the follow libraries. They should be included already if you are running on a farm machine.

- BOOST C++ libraries http://www.boost.org/doc/libs/1_48_0/index.html
- ROOT https://root.cern.ch/

# Configuration Files

The fitting runs using parameters from external cofiguration files found in config/. You will need to provide atleast four config files
all of which musy end in '.config'. Two magnet files are required for the EIC chicane setup. The files needed and an explaination of 
what each parameter is can be found below.

  parameter.config
   beam : beam is the beam energy given in GeV
   laser: laser energy given in GeV

  detector
   width  : strip width in meters
   spacing: strip spacing in meters
   theta  : detector tilts in radians 
   x : detector x position in meters
   y : detector y position in meters
   z : detector z position in meters
   
  magnet
   theta : dipole angular tilt in radians
   length: dipole length (full) in meters
   dipole: dipole strength in T
   x : magnet x position in meters
   y : magnet y position in meters
   z : magnet z position in meters


# Installation
Simply unzip the code in your favorite directory and make sure the ROOT and BOOST environment variables are set properly. The current
default is /usr/lib64/ and /usr/include/. If your libraries are in a different place you will need to either make a symbolic link
or change the directories in the Makefile.

   make clean
   make

# Running the code

   ./comptonfit

Command-line options can be found using ./comptonfit --help