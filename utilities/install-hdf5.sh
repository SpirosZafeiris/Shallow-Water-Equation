#/bin/bash
#Download HDF5 repository
git clone https://github.com/HDFGroup/hdf5.git

# enter cloned repository dir
cd ..
SWE=$PWD
cd utilities
TARGET=$SWE/lib/hdf5

## Set Compilers
##-- Intel
export CC=ifc
export CXX=ifc
export F90=ifx


cd hdf5/
mkdir build
cd build/
#configure installation
cmake -DBUILD_TESTING=OFF -DCMAKE_BUILD_TYPE=Release -DHDF5_BUILD_EXAMPLES=OFF -DHDF5_DISABLE_COMPILER_WARNINGS=ON \
      -DCMAKE_INSTALL_PREFIX=$TARGET ../

make -j4 install

# append library search path to include all hdf5 libs
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$TARGET/lib


