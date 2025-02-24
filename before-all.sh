# Script to be called as cibuildwheel before-all

set -e

CIBW_PLATFORM=$(uname)
CIBW_BUILD=$AUDITWHEEL_POLICY

echo D $CIBW_PLATFORM B $CIBW_BUILD

if [[ $CIBW_PLATFORM == "Linux" ]]; then
    echo "platform is Linux."
    if [[ $CIBW_BUILD == *"manylinux"* ]]; then
        echo "building manylinux"
        yum -q update && yum -q -y install gmp-devel
    else
        echo "building musllinux"
        apk add gmp-dev
    fi
    sh install-pari.sh
    export CC="g++"
elif [[ $CIBW_PLATFORM == "Darwin" ]]; then
    echo "platform is macOS."
    sh install-pari-msys2.sh
    brew install libomp
    export LDFLAGS="-L/opt/homebrew/opt/libomp/lib"
    export CPPFLAGS="-I/opt/homebrew/opt/libomp/include"
    export CC="g++"
elif [[ $CIBW_PLATFORM == *"MSYS_NT"* ]]; then
    echo "platform is Windows."
    # sh install-pari-msys2.sh already in workflow
else
    echo "unknown platform: $CIBW_PLATFORM"
fi
