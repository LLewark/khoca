# Script to be called after cibuildwheel even on failures to gather analysing data

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
    echo "__ZlsRNSt3 $(grep -r __ZlsRNSt3__113basic_ostreamIcNS_11char_traitsIcEEEEPK12__mpq_struct .)"
elif [[ $CIBW_PLATFORM == *"MINGW64_NT"* ]]; then
    echo "platform is Windows."
    echo "Analyse on $(pwd)"
    echo "PUI $(find . -name '*pui*')"
    echo "libparicrt64 $(find /d/a/ -name 'libparicrt64*')"
    echo "PyInit_pui $(grep -r PyInit_pui .)"
    echo "win32ctrlc $(grep -r win32ctrlc /c/msys64)"
    echo "__int64 $(grep -r __int64 /c/msys64/mingw64/lib)"
else
    echo "unknown platform: $CIBW_PLATFORM"
fi
