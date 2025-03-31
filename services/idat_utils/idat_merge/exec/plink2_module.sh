#!/bin/bash

# Set the PLINK version
PLINK_VERSION="2.3"

sudo apt update
sudo apt install -y unzip

# Set the URL for the latest PLINK 1.9 release
PLINK_URL="https://s3.amazonaws.com/plink2-assets/alpha6/plink2_linux_x86_64_20241114.zip"

# Set the VM directory where PLINK will be installed
PLINK_INSTALL_DIR="/home/levineks/bin/plink2"

# Create the PLINK installation directory
mkdir -p $PLINK_INSTALL_DIR

# Download and unzip the latest PLINK release
curl -L $PLINK_URL -o plink2.zip
unzip plink2.zip -d $PLINK_INSTALL_DIR

# Create the module file
mkdir -p /etc/modulefiles/plink2
cat <<EOF > /etc/modulefiles/plink2/${PLINK_VERSION}
#%Module
set plink_root ${PLINK_INSTALL_DIR}/plink2
prepend-path PATH $plink_root
prepend-path LD_LIBRARY_PATH $plink_root
EOF