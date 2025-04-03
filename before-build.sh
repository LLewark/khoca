# Script to be called as cibuildwheel before-build

set -e

CIBW_PLATFORM=$(uname)
CIBW_BUILD=$AUDITWHEEL_POLICY

DEBUG="W"

echo "D $CIBW_PLATFORM B $CIBW_BUILD DEBUG $DEBUG"

if [[ $CIBW_PLATFORM == "Linux" ]]; then
    echo "platform is Linux."
elif [[ $CIBW_PLATFORM == "Darwin" ]]; then
    echo "platform is macOS."
    if [[ $DEBUG == "Y" ]] || [[ $DEBUG == "M" ]]; then
        echo "GMP: $(find / 2>/dev/null -name gmp.h)"
        echo "PARI: $(find / 2>/dev/null -name pari.h)"
        echo "LIBGMP: $(find / 2>/dev/null -name libgmp.*)"
        echo "LIBPARI: $(find / 2>/dev/null -name libpari.*)"
        echo "GMPXX: $(find / 2>/dev/null -name *gmpxx.*)"
    fi
elif [[ $CIBW_PLATFORM == *"MINGW64_NT"* ]]; then
    echo "platform is Windows."
    if [[ $DEBUG == "Y" ]] || [[ $DEBUG == "W" ]]; then
        echo "GMP: $(find /c/msys64 -name gmp.h)"
        echo "PARI: $(find /d -name pari.h)"
        echo "LIBGMP: $(find /c/msys64 -name libgmp.*)"
        echo "LIBPARI: $(find /d -name libpari.*)"
        echo "GMPXX: $(find /d -name *gmpxx.*)"
    fi
else
    echo "unknown platform: $CIBW_PLATFORM"
fi
