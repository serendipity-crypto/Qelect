# MPS

```
# If permission required, please add sudo before the commands as needed

sudo apt-get update && sudo apt-get install build-essential # if needed
sudo apt-get install autoconf # if no autoconf
sudo apt-get install cmake # if no cmake
sudo apt-get install libgmp3-dev # if no gmp
sudo apt-get install libntl-dev=11.4.3-1build1 # if no ntl

# change build_path to where you want the dependency libraries installed
MPSDIR=~
BUILDDIR=$MPSDIR/MPS/build

# used for data type definitions and distribution, generators, etc
cd $MPSDIR && git clone -b v1.11.3 https://gitlab.com/palisade/palisade-release
cd palisade-release
mkdir build
cd build
cmake .. -DCMAKE_INSTALL_PREFIX=$BUILDDIR
make
make install


# a separate fork of SEAL library that overwrite some private functions, used in prior works
# we also depend on this library which is used by the prior work we are based on
cd $MPSDIR && git clone https://github.com/wyunhao/SEAL
cd SEAL
cmake -S . -B build -DCMAKE_INSTALL_PREFIX=$BUILDDIR -DSEAL_USE_INTEL_HEXL=ON 
cmake --build build
cmake --install build


cd $MPSDIR/MPS/build
cmake .. -DCMAKE_PREFIX_PATH=$BUILDDIR
mkdir ../data
mkdir ../data/perm
make
```
