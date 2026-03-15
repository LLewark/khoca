#!/usr/bin/env python

# A packaging of Khoca that allows easy installation using python's pip
#
# Sebastian Oehms, 09/24/2021
# seb.oehms@gmail.com
#
# Run "python setup.py package" and it will automatically download all the
# necessary sources and create a tar ball suitable for pip.
#
# We can upload with "twine upload -r pypi khoca-...tar.gz"

import glob, os


# Without the next line, we get an error even though we never
# use distutils as a symbol
from setuptools import distutils

from distutils.core import Extension
from setuptools import setup
from Cython.Build import cythonize
from platform import system

join = os.path.join
env = os.environ

v_file = join(os.path.dirname(__file__), '_version.py')
version = '0.9'
with open(v_file) as f:
    # read the current version from file
    code = f.read()
    exec(code, locals())

khoca_dir      = '.'
src_dir        = 'src'
bin_dir        = 'bin'
data_dir       = 'data'
converters_dir = 'converters'

khoca_pkg      = 'khoca'
bin_pkg        = join(khoca_pkg, bin_dir)
data_pkg       = join(khoca_pkg, data_dir)
converters_pkg = join(khoca_pkg, converters_dir)

pui_name = 'khoca.bin.pui'

Linux = (system() == 'Linux')
MacOS = (system() == 'Darwin')
Windows = (system() == 'Windows')

include_dirs = []
library_dirs = []
extra_objects = []
extra_compile_args = ['-c', '-D__STDC_LIMIT_MACROS', '-Wall']
extra_link_args = ['-lpthread', '-lstdc++', '-t']
libraries = []

if Linux:
    extra_compile_args += ['-fopenmp', '-std=c++11', '-shared', '-fPIC', '-O3']
    extra_link_args += ['-z defs']
    libraries = ['gmp','gmpxx','pari']
elif MacOS:
    locdir = ('Pari42')
    pari_include_dir = join(locdir, 'include')
    pari_library_dir = join(locdir, 'lib')
    hombrew_lib = '/opt/homebrew/lib/'
    include_dirs += ['/opt/homebrew/opt/libomp/include', '/opt/homebrew/include/', pari_include_dir]
    library_dirs += ['/opt/homebrew/opt/libomp/lib', hombrew_lib, pari_library_dir]
    extra_compile_args += ['-std=c++11', '-shared', '-fPIC', '-O3', '-mmacosx-version-min=10.9', '-Wno-unreachable-code', '-Wno-unreachable-code-fallthrough']
    libraries = ['gmp','gmpxx', 'pari']
elif Windows:
    locdir = ('Pari42')
    pari_include_dir = join(locdir, 'include')
    pari_library_dir = join(locdir, 'bin')
    gmp_include_dir = r'C:\msys64\usr\include'
    gmp_library_dir = r'C:\msys64\mingw64\lib'
    gmp_library_dir_bin = r'C:\msys64\mingw64\bin'

    include_dirs += [pari_include_dir, gmp_include_dir]
    library_dirs += [gmp_library_dir, gmp_library_dir_bin, pari_library_dir]
    extra_compile_args += ['/DDISABLE_INLINE', '/openmp', '/std:c11', '/LD']
    extra_link_args = [join(gmp_library_dir, 'libgmp.dll.a'), join(gmp_library_dir, 'libgmpxx.dll.a'), join(pari_library_dir, 'libpari.dll.a')]

def local_scheme(version):
    return ""

def collect_source(path, pattern, depth=0):
    """
    Find all files with the given extension under path.
    """
    result = []
    for l in range(depth + 1):
        path_components = path.split('/') + l * ['*'] + [pattern]
        result += glob.glob(join(*path_components))
    return result

def create_extension(name, cpps, includes):
    return Extension(
    name,
    sources=cpps,
    language='c++',
    include_dirs=includes,
    library_dirs=library_dirs,
    extra_objects=extra_objects,
    extra_compile_args=extra_compile_args,
    libraries=libraries,
    extra_link_args=extra_link_args)

template_cpp_files = collect_source('', '*Templates.cpp', depth=2)
cpp_files_with_templ = collect_source('', '*.cpp', depth=2)
cpp_files = [cpp for cpp in cpp_files_with_templ if cpp not in template_cpp_files]

data_files = collect_source('data/', '*')

print('cpp_files', cpp_files)
print('data_files', data_files)

include_dirs += ['', 'python_interface', 'planar_algebra', 'krasner']

pui_ext = create_extension(pui_name, cpp_files, include_dirs)

setup(name = khoca_pkg,
      version=version,
      license='GPLv2+',
      zip_safe=False,
      packages=[khoca_pkg, bin_pkg, converters_pkg, data_pkg],
      package_dir={
          khoca_pkg:      khoca_dir,
          bin_pkg:        bin_dir,
          converters_pkg: converters_dir,
          data_pkg:       data_dir
      },
      ext_modules=[pui_ext]+cythonize('src/python_interface/pui.pyx'),
      package_data={khoca_pkg: data_files},
      include_package_data=True,
      install_requires=[]
)
