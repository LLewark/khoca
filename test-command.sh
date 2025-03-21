# Script to be called as cibuildwheel before-all

set -e

CIBW_PLATFORM=$(uname)
CIBW_BUILD=$AUDITWHEEL_POLICY

echo D $CIBW_PLATFORM B $CIBW_BUILD

if [[ $CIBW_PLATFORM == "Linux" ]]; then
    echo "platform is Linux."
    echo "SO: $(find . 2>/dev/null -name *.so*)"
elif [[ $CIBW_PLATFORM == "Darwin" ]]; then
    echo "platform is macOS."
    echo "SO: $(find . 2>/dev/null -name *.so*)"
    echo "DYLIB: $(find . 2>/dev/null -name *.dylib*)"
elif [[ $CIBW_PLATFORM == *"MINGW64_NT"* ]]; then
    echo "platform is Windows."
    echo "SO: $(find . 2>/dev/null -name *.so*)"
    echo "DLL: $(find . 2>/dev/null -name *.dll*)"
else
    echo "unknown platform: $CIBW_PLATFORM"
fi
pytest
