#!/bin/bash

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
VERSION=3.11
BUILD="3"
TEMP_DIR="$(mktemp -d)"

echo "This will download and install the latest stable version of CMAKE"
echo "First we'll remove cmake"

sudo apt-get purge cmake

echo "Now going to download cmake v$VERSION.$BUILD"

cd ${TEMP_DIR}
wget https://cmake.org/files/v$VERSION/cmake-$VERSION.$BUILD.tar.gz
tar -xzvf cmake-$VERSION.$BUILD.tar.gz
cd cmake-$VERSION.$BUILD/

echo "Now going to configure cmake"

./bootstrap

echo "Now build cmake"
make -j 4
# TODO Look up checkinstall flags to make this automatic
sudo make install

echo "CMake $VERSION.$BUILD installed"

sudo ldconfig

##############################################################################
# INSTALL EIGEN
###############################################################################
EIGEN_VER=3.3.4 INSTALL_DIR="/usr/local/include"
EIGEN_RELEASE_URL="https://github.com/eigenteam/eigen-git-mirror/archive/${EIGEN_VER}.tar.gz"
TEMP_DIR="$(mktemp -d)"
WORK_DIR="$(pwd)"

# download  latest release of eigen

echo "We're going to download Eigen ${EIGEN_VER} and install to ${INSTALL_DIR}"
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

sudo ldconfig

##########################################################################
# INSTALL BOOST
#########################################################################

#########################################################################
# INSTALL CGAL
########################################################################
