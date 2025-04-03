# Script to be called as cibuildwheel test-command

set -e

CIBW_PLATFORM=$(uname)
CIBW_BUILD=$AUDITWHEEL_POLICY
PWD=$(pwd)
PROJECT=$1
VENV=$(cd ..; pwd)
KHOCA_PY=$(find $VENV -name khoca.py)
KHOCA=$(dirname $KHOCA_PY)

DEBUG="W"

echo "D $CIBW_PLATFORM B $CIBW_BUILD DEBUG $DEBUG"
echo "PWD $PWD, PROJECT $PROJECT KHOCA $KHOCA"

echo "LS $(ls -ltr $KHOCA/bin)"

if [[ $CIBW_PLATFORM == "Linux" ]]; then
    echo "platform is Linux."
    echo "SO: $(find $KHOCA 2>/dev/null -name *.so*)"
elif [[ $CIBW_PLATFORM == "Darwin" ]]; then
    echo "platform is macOS."
    echo "dylib: $(find $KHOCA 2>/dev/null -name *.dylib*)"
    if [[ $DEBUG == "Y" ]] || [[ $DEBUG == "M" ]]; then
        tar -cvf $PROJECT/../khoca.tar $KHOCA
        echo "DYLIB: $(find $PROJECT 2>/dev/null -name *.dylib*)"
        echo "otool: $(otool -L $KHOCA/bin/pui.cpython-3*-darwin.so)"
    fi
elif [[ $CIBW_PLATFORM == *"MINGW64_NT"* ]]; then
    echo "platform is Windows."
    echo "dylib: $(find $KHOCA 2>/dev/null -name *.dll*)"
    if [[ $DEBUG == "Y" ]] || [[ $DEBUG == "W" ]]; then
        echo "DLL: $(find $PROJECT 2>/dev/null -name *.dll*)"
    fi
else
    echo "unknown platform: $CIBW_PLATFORM"
fi

python $PROJECT/tests/test.py
