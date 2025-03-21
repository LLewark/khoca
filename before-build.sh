# Script to be called as cibuildwheel before-all

set -e

CIBW_PLATFORM=$(uname)
CIBW_BUILD=$AUDITWHEEL_POLICY

echo D $CIBW_PLATFORM B $CIBW_BUILD

if [[ $CIBW_PLATFORM == "Linux" ]]; then
    echo "platform is Linux."
elif [[ $CIBW_PLATFORM == "Darwin" ]]; then
    echo "platform is macOS."
    echo "GMP: $(find / 2>/dev/null -name gmp.h)"
    echo "PARI: $(find / 2>/dev/null -name pari.h)"
    echo "LIBGMP: $(find / 2>/dev/null -name libgmp.*)"
    echo "LIBPARI: $(find / 2>/dev/null -name libpari.*)"
    echo "GMPXX: $(find / 2>/dev/null -name *gmpxx.*)"
elif [[ $CIBW_PLATFORM == *"MINGW64_NT"* ]]; then
    echo "platform is Windows."
    echo "GMP: $(find /c/msys64 -name gmp.h)"
    echo "PARI: $(find /d -name pari.h)"
    echo "LIBGMP: $(find /c/msys64 -name libgmp.*)"
    echo "LIBPARI: $(find /d -name libpari.*)"
    echo "GMPXX: $(find /d -name *gmpxx.*)"
else
    echo "unknown platform: $CIBW_PLATFORM"
fi
