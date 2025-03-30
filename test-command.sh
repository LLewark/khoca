# Script to be called as cibuildwheel before-all

set -e

CIBW_PLATFORM=$(uname)
CIBW_BUILD=$AUDITWHEEL_POLICY
PWD=$(pwd)
PROJECT=$1

echo D $CIBW_PLATFORM B $CIBW_BUILD
echo "PWD $PWD, PROJECT $PROJECT"

echo "SO: $(find $PROJECT 2>/dev/null -name *.so*)"
echo "LS $(ls -ltr $PROJECT/bin)"

if [[ $CIBW_PLATFORM == "Linux" ]]; then
    echo "platform is Linux."
    cd $PROJECT && pytest && cd $PWD
elif [[ $CIBW_PLATFORM == "Darwin" ]]; then
    echo "platform is macOS."
    tar -cvf $PROJECT/../khoca.tar $PROJECT
    echo "DYLIB: $(find $PROJECT 2>/dev/null -name *.dylib*)"
    python -c "from khoca import InteractiveCalculator; InteractiveCalculator()"
elif [[ $CIBW_PLATFORM == *"MINGW64_NT"* ]]; then
    echo "platform is Windows."
    echo "DLL: $(find $PROJECT 2>/dev/null -name *.dll*)"
else
    echo "unknown platform: $CIBW_PLATFORM"
fi
