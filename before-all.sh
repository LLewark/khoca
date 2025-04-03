# Script to be called as cibuildwheel before-all

set -e

CIBW_PLATFORM=$(uname)
CIBW_BUILD=$AUDITWHEEL_POLICY

DEBUG="W"

echo "D $CIBW_PLATFORM B $CIBW_BUILD DEBUG $DEBUG"

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
elif [[ $CIBW_PLATFORM == "Darwin" ]]; then
    echo "platform is macOS."
    if [[ $DEBUG == "Y" ]] || [[ $DEBUG == "M" ]]; then
        echo "libc++ $(find /usr -name libc++*)"
        echo "llvm $(find /usr -name llvm*)"
    fi
    brew install libomp gmp
    sh install-pari.sh "$(pwd)/Pari42"
    if [[ $DEBUG == "Y" ]] || [[ $DEBUG == "M" ]]; then
        echo "__ZlsRNSt3 $(grep -r __ZlsRNSt3 /usr)"
        echo "LDD gmp 10 $(otool -L /opt/homebrew/opt/gmp/lib/libgmp.10.dylib)"
        echo "LDD gmpxx 4 $(otool -L /opt/homebrew/opt/gmp/lib/libgmpxx.4.dylib)"
        echo "LDD gmp $(otool -L /opt/homebrew/opt/gmp/lib/libgmp.dylib)"
        echo "LDD gmpxx $(otool -L /opt/homebrew/opt/gmp/lib/libgmpxx.dylib)"
    fi
elif [[ $CIBW_PLATFORM == *"MINGW64_NT"* ]]; then
    echo "platform is Windows."
    bash install-pari.sh "$(pwd)/Pari42"
    ln -fs /ucrt64/include/gmp.h /usr/include
    if [[ $DEBUG == "Y" ]] || [[ $DEBUG == "W" ]]; then
        echo "LDD u $(ldd /c/msys64/ucrt64/lib/libgmp.a)"
        echo "LDD u dll.a $(ldd /c/msys64/ucrt64/lib/libgmp.dll.a)"
        echo "LDD u dll $(ldd /c/msys64/ucrt64/bin/libgmp-10.dll)"
        echo "LDD u xx dll $(ldd /c/msys64/ucrt64/bin/libgmpxx-4.dll)"
        echo "LDD m $(ldd /c/msys64/mingw64/lib/libgmp.a)"
        echo "LDD m dll.a $(ldd /c/msys64/mingw64/lib/libgmp.dll.a)"
        echo "LDD m dll $(ldd /c/msys64/mingw64/bin/libgmp-10.dll)"
        echo "LDD m xx dll $(ldd /c/msys64/mingw64/bin/libgmpxx-4.dll)"
    fi
else
    echo "unknown platform: $CIBW_PLATFORM"
fi
