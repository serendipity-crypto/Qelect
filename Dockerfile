FROM ubuntu:22.04
# Install dependencies
RUN apt-get update && apt-get install -y \
    build-essential \
    autoconf \
    cmake \
    libgmp3-dev \
    libntl-dev \
    graphviz \
    git
# Set working directory
WORKDIR /home/ubuntu/MPS

# Clone the MPS repository
RUN git clone https://github.com/wyunhao/MPS /home/ubuntu/MPS
# Install PALISADE library
RUN git clone -b v1.11.9 https://gitlab.com/palisade/palisade-release \
    && cd palisade-release && mkdir build && cd build \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/home/ubuntu/MPS/build \
    && make -j8 && make install
# Install SEAL library
RUN git clone https://github.com/wyunhao/SEAL \
    && cd SEAL && cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/home/ubuntu/MPS/build -DSEAL_USE_INTEL_HEXL=ON \
    && cmake --build build -- -j8 && cmake --install build
# Build the MPS project
RUN cd /home/ubuntu/MPS && cd build \
    && mkdir ../data && mkdir ../data/perm \
    && cmake .. -DCMAKE_PREFIX_PATH=/home/ubuntu/MPS/build \
    && make -j8
# Set entrypoint
ENTRYPOINT ["/home/ubuntu/MPS/build/mps"]
