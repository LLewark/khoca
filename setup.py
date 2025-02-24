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

join = os.path.join

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

def local_scheme(version):
    return ""

def collect_source(path, pattern, depth=0):
    """
    Find all files with the given extension under path.
    """
    result = []
    for l in range(depth + 1):
        path_components = path.split('/') + l * ['*'] + [pattern]
        result += glob.glob(os.path.join(*path_components))
    return result

def create_extension(name, cpps, includes):
    return Extension(
    name,
    sources = cpps,
    include_dirs = includes,
    language = 'c++',
    extra_compile_args=['-c', '-std=c++11', '-D__STDC_LIMIT_MACROS', '-Wall', '-shared', '-fPIC', '-O3', '-fopenmp'],
    libraries = ['gmp','gmpxx','pari'],
    extra_link_args =['-z defs', '-lpthread', '-lstdc++', '-t'])

template_cpp_files = collect_source('', '*Templates.cpp', depth=2)
cpp_files_with_templ = collect_source('', '*.cpp', depth=2)
cpp_files = [cpp for cpp in cpp_files_with_templ if cpp not in template_cpp_files]

data_files = collect_source('data/', '*')

print('cpp_files', cpp_files)
print('data_files', data_files)

include_dirs = ['', 'python_interface', 'planar_algebra', 'krasner']

pui_ext = create_extension(pui_name, cpp_files, include_dirs)


long_description = """
Khoca is a packaging of Lukas Lewark's Khovanov homology calculator
[khoca](https://github.com/LLewark/khoca/)  for easy installation in
[sage](http://www.sagemath.org/).

From the original page:

"""

this_directory = os.path.abspath(os.path.dirname(__file__))
with open(os.path.join(this_directory, 'README.md'), encoding='utf-8') as f:
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
