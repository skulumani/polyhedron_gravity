#!/bin/bash

CGAL_VER='CGAL-4.12.1'
TEMP_DIR="$(mktemp -d)"
CGAL_RELEASE_URL='https://github.com/CGAL/cgal/releases/download/releases%2F'${CGAL_VER}'/'${CGAL_VER}'.tar.xz'

# check if the temp dir was created
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

echo "Installing some dependencies"
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    sudo apt-get update
    sudo apt-get install cmake build-essential libgmp-dev libmpfr-dev zlib1g-dev
    sudo apt-get install libgl1 libstdc++6 libgcc1 libc6 libntl-dev 
elif [[ "$OSTYPE" == "darwin"* ]]; then
        # Mac OSX
fi
# install CGAL for Python and build it from source for C++
echo "Downloading ${CGAL_VER}"

# download the source tarball
cd ${TEMP_DIR}
wget ${CGAL_RELEASE_URL} $TEMP_DIR
tar xf ${CGAL_VER}.tar.* 

echo "Installing ${CGAL_VER}"
cd $CGAL_VER
cmake -DWITH_examples=ON -DWITH_demos=ON -DWITH_CGAL_Qt5=OFF .
make -j 4
# make examples
# make demos
if [[ "$OSTYPE" == "linux-gnu" ]]; then
    sudo checkinstall make install
fi
