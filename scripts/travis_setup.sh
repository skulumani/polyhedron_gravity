#!/bin/bash

# setup script to install dependencies on Travis
# If it does exist then only update the conda environment rather than create it

TEMP_DIR="$(mktemp -d)"
WORK_DIR="$(pwd)"
DIR="$(pwd)"

#ANACONDA
ANACONDA_ENV="asteroid"
ANACONDA_PATH="$HOME/anaconda3/envs"

# EIGEN
EIGEN_VER=3.3.4
EIGEN_INSTALL_DIR="/usr/local/include"
EIGEN_RELEASE_URL="https://github.com/eigenteam/eigen-git-mirror/archive/${EIGEN_VER}.tar.gz"

#HDF5
HDF5_VER=1.10.2
HDF5_RELEASE_URL="https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.2/src/hdf5-${HDF5_VER}.tar.gz"
HDF5_INSTALL_DIR="/usr/local/include"

# BOOST
BOOST_VER=1.67.0
BOOST_SHA256_SUM="8aa4e330c870ef50a896634c931adf468b21f8a69b77007e45c444151229f665"
BOOST_URL="https://dl.bintray.com/boostorg/release/${BOOST_VER}/source/boost_1_67_0.tar.gz"
BOOST_INSTALL_DIR="/usr/local"

# CGAL
CGAL_VER='CGAL-4.12'
CGAL_RELEASE_URL='https://github.com/CGAL/cgal/releases/download/releases%2F'${CGAL_VER}'/'${CGAL_VER}'.tar.xz'

# CMAKE
CMAKE_VERSION=3.11
CMAKE_BUILD="3"

if [[ ! "$TEMP_DIR" || ! -d "$TEMP_DIR" ]]; then
	echo "Could not create temp dir"
	exit 1
fi

# delete the temp directory on cleanup
function cleanup {
    rm -rf "$TEMP_DIR"
    echo "Deleted temp working directory $TEMP_DIR"
}

# trap cleanup EXIT

echo "Installing some dependencies"
sudo apt-get -qq update
sudo apt-get install -y build-essential g++ libicu-dev libbz2-dev autotools-dev
sudo apt-get install -y build-essential libgmp-dev libmpfr-dev zlib1g-dev libgmp-dev libmpc3 libmpfr4
sudo apt-get install -y libstdc++6 libgcc1 libc6 libntl-dev 
sudo apt-get purge -y cmake
sudo apt-get remove -y cmake

echo "Downloading and Installing Miniconda"
## get miniconda installed
## check if anaconda3 directory exists
if [ -f "$HOME/anaconda3/bin/conda" ]; then
    echo "Anaconda already installed"
else
    echo "Anaconda not installed"
    wget https://repo.continuum.io/miniconda/Miniconda3-latest-Linux-x86_64.sh -O $HOME/miniconda.sh
    bash $HOME/miniconda.sh -b -u -p $HOME/anaconda3
fi

export PATH="$HOME/anaconda3/bin:$PATH"
hash -r

echo "Configuring Conda settings"
## configure conda
conda config --set always_yes yes --set changeps1 no
conda update conda
echo ". /home/travis/anaconda3/etc/profile.d/conda.sh" >> ~/.bashrc
export PATH="$HOME/anaconda3/bin:$PATH"
hash -r 
source $HOME/anaconda3/bin/activate

echo "Downloading the shape models"
## download the shape model
wget https://github.com/skulumani/asteroid_dumbbell/releases/download/v0.3/shape_model.tar.gz -O ./data/shape_model/shape_model.tar.gz
tar xf ./data/shape_model/shape_model.tar.gz -C ./data/shape_model

echo "Creating the asteroid environment"
# setup development enviornment
if [ -d "$HOME/anaconda3/envs/asteroid" ]; then
    echo "asteroid enviornment exists. Just update"
    conda env update --name asteroid --file ./utilities/asteroid.yml
else
    echo "No asteroid enviornment"
    conda env create --name asteroid --file ./utilities/asteroid.yml

fi

conda activate asteroid

echo "Anaconda Setup is complete"

#CMAKE
sudo apt-get -y purge cmake
sudo apt-get -y remove cmake

echo "Now going to download cmake v$CMAKE_VERSION.$CMAKE_BUILD"
cd ${TEMP_DIR}
wget https://cmake.org/files/v${CMAKE_VERSION}/cmake-${CMAKE_VERSION}.${CMAKE_BUILD}-Linux-x86_64.sh -O install_cmake.sh
chmod +x install_cmake.sh
sudo ./install_cmake.sh --prefix=/usr/local --skip-license
## EIGEN
echo "We're going to download Eigen ${EIGEN_VER} and install to ${EIGEN_INSTALL_DIR}"
cd ${TEMP_DIR}
mkdir ${EIGEN_VER}
wget ${EIGEN_RELEASE_URL} -O ${TEMP_DIR}/${EIGEN_VER}.tar.gz
tar -xvzf ${EIGEN_VER}.tar.gz -C ./${EIGEN_VER} --strip-components=1

echo "Going to install Eigen using CMake"
cd ${EIGEN_VER}
mkdir build
cd build
cmake ..
sudo make install

echo "Eigen installed"

## HDF5
echo "We're going to download HDF5 ${HDF5_VER} and install to ${HDF5_INSTALL_DIR}"
cd ${TEMP_DIR}
mkdir ${HDF5_VER}
wget ${HDF5_RELEASE_URL} -O ${TEMP_DIR}/${HDF5_VER}.tar.gz
tar -xvzf ${HDF5_VER}.tar.gz -C ./${HDF5_VER} --strip-components=1

echo "Going to install HDF5 using the configure script"
cd ${HDF5_VER}
bash ./configure --prefix=/usr/local/hdf5 --enable-cxx 
make -j5
# make check
sudo make install

echo "Finished installing HDF5"

## BOOST
echo "Now downloading Boost"
cd ${TEMP_DIR}
mkdir boost
wget ${BOOST_URL} -O ${TEMP_DIR}/boost.tar.gz
# verify sha256 sum
if ! echo "${BOOST_SHA256_SUM} boost.tar.gz" | sha256sum -c; then
    echo "Checksum does not match. Aborting!!"
    exit 1
fi

tar -xzf boost.tar.gz -C ./boost --strip-components=1

echo "Remove old boost"
sudo rm -rf /usr/include/boost
sudo rm /usr/local/lib/libboost*

echo "Now installing Boost and compiled libraries"
cd boost
./bootstrap.sh --prefix=${BOOST_INSTALL_DIR} --with-libraries=all --with-python=$HOME/anaconda3/bin/python3

sudo ./b2 -j 4 install

echo "Boost and Boost-Python are installed to $BOOST_INSTALL_DIR"

## CGAL
echo "Now installing CGAL"
cd ${TEMP_DIR}
wget ${CGAL_RELEASE_URL} $TEMP_DIR
tar xf ${CGAL_VER}.tar.* 
cd $CGAL_VER
cmake -DWITH_examples=ON -DWITH_demos=ON -DWITH_CGAL_Qt5=OFF .
make -j 4
sudo make install

cd $DIR
