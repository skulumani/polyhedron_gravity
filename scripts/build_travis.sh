#!/bin/bash

TEMP_DIR="$(mktemp -d)"
if [[ ! "$TEMP_DIR" || ! -d "$TEMP_DIR" ]]; then
	echo "Could not create temp dir"
	exit 1
fi

# delete the temp directory on cleanup
function cleanup {
    rm -rf "$TEMP_DIR"
    echo "Deleted temp working directory $TEMP_DIR"
}

trap cleanup EXIT

#############################################################################
# INSTALL CMAKE
############################################################################
# 3.10.2 has an issue with the Finding boost 1.66
# VERSION=3.11
# BUILD="3"

# echo "This will download and install the latest stable version of CMAKE"
# echo "First we'll remove cmake"

# sudo apt-get purge cmake

# echo "Now going to download cmake v$VERSION.$BUILD"

# cd ${TEMP_DIR}
# wget https://cmake.org/files/v$VERSION/cmake-$VERSION.$BUILD.tar.gz
# tar -xzvf cmake-$VERSION.$BUILD.tar.gz
# cd cmake-$VERSION.$BUILD/

# echo "Now going to configure cmake"

# ./bootstrap

# echo "Now build cmake"
# make -j 4
# # TODO Look up checkinstall flags to make this automatic
# sudo make install

# echo "CMake $VERSION.$BUILD installed"

# sudo ldconfig

##############################################################################
# INSTALL EIGEN
###############################################################################
# EIGEN_VER=3.3.4 INSTALL_DIR="/usr/local/include"
# EIGEN_RELEASE_URL="https://github.com/eigenteam/eigen-git-mirror/archive/${EIGEN_VER}.tar.gz"
# TEMP_DIR="$(mktemp -d)"
# WORK_DIR="$(pwd)"

# # download  latest release of eigen

# echo "We're going to download Eigen ${EIGEN_VER} and install to ${INSTALL_DIR}"
# cd ${TEMP_DIR}
# mkdir ${EIGEN_VER}
# wget ${EIGEN_RELEASE_URL} -O ${TEMP_DIR}/${EIGEN_VER}.tar.gz
# tar -xvzf ${EIGEN_VER}.tar.gz -C ./${EIGEN_VER} --strip-components=1

# echo "Going to install Eigen using CMake"
# cd ${EIGEN_VER}
# mkdir build
# cd build
# cmake ..
# sudo make install

# echo "Eigen installed"

# sudo ldconfig

##########################################################################
# INSTALL BOOST
#########################################################################
# BOOST_VER=1.68.0
# BOOST_SHA256_SUM="da3411ea45622579d419bfda66f45cd0f8c32a181d84adfa936f5688388995cf"
# BOOST_URL="https://dl.bintray.com/boostorg/release/${BOOST_VER}/source/boost_1_68_0.tar.gz"
# TEMP_DIR="$(mktemp -d)"
# INSTALL_DIR="/usr/local"

# # This will install the latest Boost and the Boost Python library
# echo "We're going to download and install the latest Boost"

# echo "Installing some dependencies"
# sudo apt-get install -y build-essential libicu-dev libbz2-dev autotools-dev

# echo "Now downloading Boost"
# cd ${TEMP_DIR}
# mkdir boost
# wget ${BOOST_URL} -O ${TEMP_DIR}/boost.tar.gz
# # verify sha256 sum
# if ! echo "${BOOST_SHA256_SUM} boost.tar.gz" | sha256sum -c; then
#     echo "Checksum does not match. Aborting!!"
#     exit 1
# fi

# tar -xzf boost.tar.gz -C ./boost --strip-components=1

# echo "Now installing Boost and compiled libraries"
# cd boost
# ./bootstrap.sh --prefix=${INSTALL_DIR} --with-libraries=all


# sudo ./b2 -j 4 install

# echo "Boost and are installed to $INSTALL_DIR"

#########################################################################
# INSTALL CGAL
########################################################################
CGAL_VER='CGAL-4.12.1'
TEMP_DIR="$(mktemp -d)"
CGAL_RELEASE_URL='https://github.com/CGAL/cgal/releases/download/releases%2F'${CGAL_VER}'/'${CGAL_VER}'.tar.xz'

echo "Installing some dependencies"
sudo apt-get update
sudo apt-get install build-essential libgmp-dev libmpfr-dev zlib1g-dev
sudo apt-get install libgl1 libstdc++6 libgcc1 libc6 libntl-dev 
# install CGAL for Python and build it from source for C++
echo "Downloading ${CGAL_VER}"

# download the source tarball
cd ${TEMP_DIR}
wget ${CGAL_RELEASE_URL} $TEMP_DIR
tar xf ${CGAL_VER}.tar.* 

echo "Installing ${CGAL_VER}"
cd $CGAL_VER
cmake -DWITH_examples=OFF -DWITH_demos=OFF -DWITH_CGAL_Qt5=OFF .
make -j 4
sudo make install
