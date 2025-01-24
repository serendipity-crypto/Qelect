



# MPS 
## The instructions to run the project via Dockerfile:
```
# install docker
sudo apt-get update
sudo apt-get install -y ca-certificates curl
# might not needed, ignore for now: sudo apt-get install -y -m 0755 -d /etc/apt/keyrings
sudo curl -fsSL https://download.docker.com/linux/ubuntu/gpg -o /etc/apt/keyrings/docker.asc
sudo chmod a+r /etc/apt/keyrings/docker.asc

echo \
  "deb [arch=$(dpkg --print-architecture) signed-by=/etc/apt/keyrings/docker.asc] https://download.docker.com/linux/ubuntu \
  $(. /etc/os-release && echo "$UBUNTU_CODENAME") stable" | \
  sudo tee /etc/apt/sources.list.d/docker.list > /dev/null

sudo apt-get update
sudo apt-get install -y docker-ce docker-ce-cli containerd.io docker-buildx-plugin docker-compose-plugin


# execute the MPE
# download the zip file and unzip the folder under MPS
cd MPS
sudo docker build --no-cache -t mps_project .
sudo docker run mps_project
```

An example output is as follows:
```
After param generation.
Initial noise budget: 613
Lefted noise budget: 47
Total broadcast communication size (for each party): 2147 MB.
---------------------------------
Total number of parties: 1024
Preprocessed time      : 136201928 us.
Total time             : 5921790 us.
```
The local computation time for a single party, with party size = 1024, is ~5.9 seconds.

## The instructions to run the project via cmake:
```
# If permission required, please add sudo before the commands as needed

sudo apt-get update && sudo apt-get install build-essential # if needed
sudo apt-get install autoconf # if no autoconf
sudo apt-get install cmake # if no cmake
sudo apt-get install libgmp3-dev # if no gmp
sudo apt-get install libntl-dev=11.4.3-1build1 # if no ntl, or sudo apt-get install libntl-dev is version not found

# change build_path to where you want the dependency libraries installed
MPSDIR=~
BUILDDIR=$MPSDIR/MPS/build

# used for data type definitions and distribution, generators, etc
cd $MPSDIR && git clone -b v1.11.9 https://gitlab.com/palisade/palisade-release
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
