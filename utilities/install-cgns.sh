#/bin/bash
# Downloads CGNS master branch
git clone  https://github.com/CGNS/CGNS.git

cd ..
SWE=$PWD
cd utilities
HDF5_DIR=$SWE/lib/hdf5/cmake
TARGET=$SWE/lib/cgns
# Intel fortran
export CC=ifc
export FC=ifx

# Set paths for installation
# Configure and make
cd CGNS
mkdir build
cd build

cmake  -DCMAKE_BUILD_TYPE=Release -DCGNS_ENABLE_HDF5=ON -DHDF5_NEED_MPI=OFF \
       -DCGNS_ENABLE_PARALLEL=OFF -DCGNS_ENABLE_FORTRAN=ON -DHDF5_DIR=$HDF5_DIR -DCMAKE_INSTALL_PREFIX=$TARGET ../

make -j4 install
