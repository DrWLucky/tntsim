:: TEXANSIM ::

   A GEANT4 simulation package for the Texas Array for Neutrons.
   Authors:
     G. Christian
     gchristian@tamu.edu


Installation Instructions:
  Prerequisites:
     GEANT4
        I use version 4.10.01.p02 and have not tested with any other version of GEANT.
       	Geant must be installed with the GDML option: -DGEANT4_USE_GDML=ON when running
	cmake this requires installation of the xerces-c package which is freely available
	online. See the GEANT manual for more info.

     xerces-c (see above)

     ROOT
        Output is written into ROOT files (http://root.cern.ch). I use v5.34/26. Older versions
	should be fine, but at the moment I would not recomment using anything in ROOT 6.
	
          

     Following the protocol of GEANT4, this package uses cmake to set up the build environment:
          cd build
	  edit cmake.sh to reflect your GEANT4 installation location
	  ./cmake.sh
	  make
	  make install (optional)



Coding Conventions:
  - Use texansim namespace to enclose project classes (often alias to txs in .cc files)
  - Prefix file names with Texan, i.e. TexanDetectorConstruction.hh