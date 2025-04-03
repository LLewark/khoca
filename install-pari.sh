# Helper script to install PARI for github workflows

# Exit on error
set -e

if [ "$PARI_VERSION" = "" ]; then
    PARI_VERSION="pari-2.17.1"
fi

if [ "$1" = "" ]; then
    PREFIX="/usr"
else
    PREFIX=$1
fi

PARI_URL="http://pari.math.u-bordeaux.fr/pub/pari/unix"

# Download PARI sources
curl --no-verbose "$PARI_URL/$PARI_VERSION.tar.gz" -o pari.tgz

# Install
mkdir -p Pari42
tar xzf pari.tgz -C Pari42
cd Pari42/*
./Configure --prefix=$PREFIX
make gp
make install
