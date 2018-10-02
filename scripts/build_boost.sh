#!/bin/bash
BOOST_VER=1.67.0
BOOST_SHA256_SUM="8aa4e330c870ef50a896634c931adf468b21f8a69b77007e45c444151229f665"
BOOST_URL="https://dl.bintray.com/boostorg/release/${BOOST_VER}/source/boost_1_67_0.tar.gz"
TEMP_DIR="$(mktemp -d)"
INSTALL_DIR="/usr/local"

# This will install the latest Boost and the Boost Python library
echo "We're going to download and install the latest Boost"
read -p "Press Enter to continue"

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
sudo apt-get install -y build-essential g++ libicu-dev libbz2-dev autotools-dev

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

# copy to /usr/local/include
# echo "Now copying to ${INSTALL_DIR}/boost"
# sudo mv boost/boost ${INSTALL_DIR}

echo "Remove old boost"
sudo rm -rf /usr/include/boost
sudo rm /usr/local/lib/libboost*

echo "Now installing Boost and compiled libraries"
cd boost
./bootstrap.sh --prefix=${INSTALL_DIR} --with-libraries=all --with-python=$HOME/anaconda3/bin/python3

echo "Now you need to add the following line:"

echo "using python : 3.6 : /home/shankar/anaconda3 : /home/shankar/anaconda3/include/python3.6m ;" 
echo "Located inside ${TEMP_DIR}/boost/project-config.jam"

read -p "Press enter when done"

sudo ./b2 -j 4 install

echo "Boost and Boost-Python are installed to $INSTALL_DIR"
read -p "Press Enter to exit"
