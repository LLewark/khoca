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

khoca_dir      = '.'
src_dir        = 'src'
bin_dir        = 'bin'
data_dir       = 'data'
converters_dir = 'converters'

khoca_pkg      = 'khoca'
bin_pkg        = join(khoca_pkg, bin_dir)
data_pkg       = join(khoca_pkg, data_dir)
converters_pkg = join(khoca_pkg, converters_dir)

pui_name = join(bin_pkg, 'pui')

paridir = join('libcache', 'pari')
gmpdir = join('libcache', 'gmp')
pari_include_dir = join(paridir, 'include')
gmp_include_dir = join(gmpdir, 'include')
pari_library_dir = join(paridir, 'lib')
pari_static_library = join(pari_library_dir, 'libpari.a')
gmp_library_dir = join(gmpdir, 'lib')
gmp_static_library = join(gmp_library_dir, 'libgmp.a')

include_dirs = []
extra_objects = []
extra_compile_args = ['-c', '-D__STDC_LIMIT_MACROS', '-Wall']
extra_link_args = ['-lpthread', '-lstdc++', '-t']
libraries = []
if system() == 'Linux':
    extra_compile_args += ['-fopenmp', '-std=c++11', '-shared', '-fPIC', '-O3']
    extra_link_args += ['-z defs']
    libraries = ['gmp','gmpxx','pari']
elif system() == 'Darwin':
    include_dirs += [gmp_include_dir, pari_include_dir, '/opt/homebrew/opt/libomp/include']
    extra_compile_args += ['-std=c++11', '-shared', '-fPIC', '-O3', '-mmacosx-version-min=10.9', '-Wno-unreachable-code', '-Wno-unreachable-code-fallthrough']
    extra_link_args += ['-L/opt/homebrew/opt/libomp/lib', '-L/opt/homebrew/lib/']
    libraries = [pari_static_library, gmp_static_library]
elif system() == 'Windows':
    include_dirs += [gmp_include_dir, pari_include_dir]
    include_dirs += [r'C:\Program Files (x86)\Windows Kits\10\Include\10.0.22000.0\um',
                     r'C:\Program Files (x86)\Windows Kits\10\Include\10.0.22000.0\ucrt',
                     r'C:\Program Files (x86)\Windows Kits\10\Include\10.0.22000.0\shared']
    extra_compile_args += ['/DDISABLE_INLINE', '/openmp', '/std:c11', '/LD']
    extra_link_args += [join('Windows', 'crt', 'libparicrt64.a'), 'advapi32.lib', 'legacy_stdio_definitions.lib', join('Windows', 'crt', 'get_output_format64.o')]
    extra_objects += [r'C:\Program Files (x86)\Windows Kits\10\Lib\10.0.22000.0\um\x64\Uuid.lib',
                     r'C:\Program Files (x86)\Windows Kits\10\Lib\10.0.22000.0\um\x64\kernel32.lib',
                     r'C:\Program Files (x86)\Windows Kits\10\Lib\10.0.22000.0\ucrt\x64\ucrt.lib',
                     r'C:\msys64\ucrt64\lib\gcc\x86_64-w64-mingw32\14.2.0\libgcc.a']
    libraries = [pari_static_library, gmp_static_library]

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


long_description = """
Khoca is a packaging of Lukas Lewark's Khovanov homology calculator
[khoca](https://github.com/LLewark/khoca/)  for easy installation in
[sage](http://www.sagemath.org/).

From the original page:

"""

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(join(this_directory, 'README.md'), encoding='utf-8') as f:
    long_description += f.read()




setup(name = khoca_pkg,
      zip_safe = False,
      description = 'Khoca as pip installable package',
      long_description = long_description,
      long_description_content_type='text/markdown',
      keywords = 'Knot theory, Khovanov homology',
      classifiers = [
           'Development Status :: 3 - Alpha',
           'Intended Audience :: Science/Research',
           'License :: OSI Approved :: GNU General Public License v2 or later (GPLv2+)',
           'Operating System :: POSIX :: Linux',
           'Programming Language :: C++',
           'Programming Language :: Python',
           'Programming Language :: Cython',
           'Topic :: Scientific/Engineering :: Mathematics',
      ],
      author = 'Sebastian Oehms',
      author_email = 'seb.oehms@gmail.com',
      url = 'https://github.com/soehms/khoca/',
      license='GPLv2+',
      packages = [khoca_pkg, bin_pkg, converters_pkg, data_pkg],
      package_dir = {
          khoca_pkg:      khoca_dir,
          bin_pkg:        bin_dir,
          converters_pkg: converters_dir,
          data_pkg:       data_dir
      },
      ext_modules = [pui_ext]+cythonize('src/python_interface/pui.pyx'),
      package_data = {khoca_pkg: data_files},
      include_package_data=True,
      use_scm_version={"local_scheme": local_scheme},
      setup_requires=['setuptools_scm'],
      install_requires=[]
)
