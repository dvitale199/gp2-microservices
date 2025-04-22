#!/bin/bash

# Set the PLINK version
PLINK_VERSION="1.9"

sudo apt update
sudo apt install -y unzip

# Set the URL for the latest PLINK 1.9 release
PLINK_URL="https://s3.amazonaws.com/plink1-assets/plink_linux_x86_64_20220402.zip"

# Set the VM directory where PLINK will be installed
PLINK_INSTALL_DIR="./plink1.9"

# Create the PLINK installation directory
mkdir -p $PLINK_INSTALL_DIR

# Download and unzip the latest PLINK 1.9 release
curl -L $PLINK_URL -o plink1.9.zip
unzip plink1.9.zip -d $PLINK_INSTALL_DIR

# Create the module file
mkdir -p /etc/modulefiles/plink1.9
cat <<EOF > /etc/modulefiles/plink1.9/${PLINK_VERSION}
#%Module
set plink_root ${PLINK_INSTALL_DIR}
prepend-path PATH \$plink_root
prepend-path LD_LIBRARY_PATH \$plink_root
EOF