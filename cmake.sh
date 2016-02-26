#!/bin/bash

## Edit to reflect your GEANT4 library directory
MY_G4LIBDIR=$G4INSTALL/lib64/Geant4-10.1.2


mkdir -p build
cd build/
cmake -DGeant4_DIR=$MY_G4LIBDIR $PWD/../ -DCMAKE_INSTALL_PREFIX=${HOME}/install
