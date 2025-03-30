# Script to be called as cibuildwheel before-all

set -e

CIBW_PLATFORM=$(uname)
CIBW_BUILD=$AUDITWHEEL_POLICY

echo D $CIBW_PLATFORM B $CIBW_BUILD

if [[ $CIBW_PLATFORM == "Linux" ]]; then
    echo "platform is Linux."
elif [[ $CIBW_PLATFORM == "Darwin" ]]; then
    echo "platform is macOS."
    echo "Analyse on $(pwd)"
    echo "PUI $(find . -name '*pui*')"
    echo "PyInit_pui $(grep -r PyInit_pui .)"
elif [[ $CIBW_PLATFORM == *"MINGW64_NT"* ]]; then
    echo "platform is Windows."
    echo "Analyse on $(pwd)"
    echo "PUI $(find . -name '*pui*')"
    echo "libparicrt64 $(find /d/a/ -name 'libparicrt64*')"
    echo "PyInit_pui $(grep -r PyInit_pui .)"
    echo "__gmpz_fdiv_q_ui $(grep -r __gmpz_fdiv_q_ui .)"
    echo "__gmpz_fdiv_q_ui $(grep -r __gmpz_fdiv_q_ui /c/msys64/ucrt64/)"
    echo "pari_err_last $(grep -r pari_err_last .)"
    echo "win32ctrlc $(grep -r win32ctrlc .)"
    echo "__int64 $(grep -r __int64 .)"
else
    echo "unknown platform: $CIBW_PLATFORM"
fi
