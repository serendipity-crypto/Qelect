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
WORKDIR /home/ubuntu/qelect

RUN git config --global http.proxy "http://10.210.3.114:20172"
RUN git config --global https.proxy "http://10.210.3.114:20172"

# Clone the qelect repository
RUN git clone https://github.com/serendipity-crypto/Qelect /home/ubuntu/qelect
# Install PALISADE library
RUN git clone -b v1.11.9 https://gitlab.com/palisade/palisade-release \
    && cd palisade-release && mkdir build && cd build \
    && cmake .. -DCMAKE_INSTALL_PREFIX=/home/ubuntu/qelect/build \
    && make -j && make install
# Install SEAL library
RUN git clone https://github.com/wyunhao/SEAL \
    && cd SEAL && cmake -S . -B build -DCMAKE_INSTALL_PREFIX=/home/ubuntu/qelect/build -DSEAL_USE_INTEL_HEXL=ON \
    && cmake --build build -- -j && cmake --install build
# Build the qelect project
RUN cd /home/ubuntu/qelect && cd build \
    && mkdir ../data && mkdir ../data/perm \
    && cmake .. -DCMAKE_PREFIX_PATH=/home/ubuntu/qelect/build \
    && make -j
# Set entrypoint
ENTRYPOINT ["/home/ubuntu/qelect/build/mps"]
