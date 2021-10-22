CC=g++
PYTHONVERSION=3.8
GENERALCFLAGS=-c -std=c++11 -I/usr/include/python$(PYTHONVERSION) -D__STDC_LIMIT_MACROS -Wall -Wextra -march=native

STATIC_DIRECTIVE=-lpython$(PYTHONVERSION) -lgmpxx -lgmp -L/usr/local/lib -lpari
# To compile a version of the library that contains static copies of mpir and pari, set variable STATIC to 1.
ifeq ($(STATIC),1)
  STATIC_DIRECTIVE=-Wl,--no-as-needed -ldl -lutil -static-libgcc -lpython$(PYTHONVERSION) -L/usr/local/lib /usr/local/lib/libmpirxx.a /usr/local/lib/libmpir.a /usr/local/lib/libpari.a
endif

CFLAGSPYXO=-Wno-write-strings -Wimplicit-fallthrough=0 -fno-strict-aliasing -fwrapv -pthread
SHAREDCFLAGS=-shared -fPIC

# HOW TO DEBUG:
# Turn on Debugging here in Makefile; change the switches in shared.h (and maybe compile only with one bitsize for speed); and turn on debugging-bool-constant in khoroho.py.
# If multithreading is not disabled, there will be plenty of errors, because sanity checks are done while other threads are active.

# FAST COMPILE
# CFLAGS=$(GENERALCFLAGS) $(SHAREDCFLAGS) -O0 -ggdb -g3 -fopenmp
# COMPILERFLAGSINFO=\"Fast compilation\"

# DEBUGGING
# CFLAGS=$(GENERALCFLAGS) $(SHAREDCFLAGS) -Og -ggdb -g3 -fopenmp
# COMPILERFLAGSINFO=\"Debug\"

# FAST, BUT DEBUG INFO
# CFLAGS=$(GENERALCFLAGS) $(SHAREDCFLAGS) -O3 -g3 -fopenmp
# COMPILERFLAGSINFO=\"Fast program, but including debug info\"

# NOT DEBUGGING
CFLAGS=$(GENERALCFLAGS) $(SHAREDCFLAGS) -O3 -fopenmp
COMPILERFLAGSINFO=\"Release\"

CMPIRFLAGS=-I/usr/local/include
LDFLAGS=$(STATIC_DIRECTIVE) $(SHAREDCFLAGS) -z defs -lpthread -lstdc++ -Wl,-t

PYX=src/python_interface/pui.pyx
PYXCPP=src/python_interface/pui.cpp
PYXO=src/python_interface/pui.o
PYI=src/python_interface/pythonInterface.cpp
PYIO=src/python_interface/pythonInterface.o
SOURCES=src/planar_algebra/planar_algebra.cpp src/shared.cpp src/krasner/krasner.cpp src/planar_algebra/sparsemat.cpp src/planar_algebra/coefficient_rings.cpp src/planar_algebra/smith.cpp src/planar_algebra/coefficient_rings_explicitTemplates.cpp
OBJECTS=src/planar_algebra/planar_algebra.o src/shared.o src/krasner/krasner.o src/planar_algebra/sparsemat.o src/planar_algebra/coefficient_rings.o src/planar_algebra/smith.o
EXECUTABLE=bin/pui.so


.PHONY : clean

all:	$(EXECUTABLE)

$(EXECUTABLE): $(SOURCES) $(OBJECTS) $(PYX) $(PYXCPP) $(PYXO) $(PYI) $(PYIO) Makefile
	$(CC) $(OBJECTS) $(PYIO) $(PYXO) -o $@ $(LDFLAGS)

$(PYXCPP): $(PYX) Makefile
	cython -3 --cplus $(PYX)

$(PYXO): $(PYXCPP) src/python_interface/pythonInterface.h src/shared.h Makefile
	$(CC) $(CFLAGS) $(CFLAGSPYXO) $< -o $@

src/python_interface/pythonInterface.o: src/python_interface/pythonInterface.cpp src/planar_algebra/planar_algebra.h src/krasner/krasner.h src/shared.h src/planar_algebra/sparsemat.h src/planar_algebra/coefficient_rings.h Makefile
	echo "#define COMPILERFLAGSINFO $(COMPILERFLAGSINFO)" > src/compilerFlagsInfo.cpp
	$(CC) $(CMPIRFLAGS) $(CFLAGS) $< -o $@

src/planar_algebra/planar_algebra.o: src/planar_algebra/planar_algebra.cpp src/planar_algebra/planar_algebra.h src/shared.h src/planar_algebra/explicitTemplates.cpp src/planar_algebra/sparsemat.h src/planar_algebra/coefficient_rings.h src/krasner/krasner.h Makefile
	$(CC) $(CMPIRFLAGS) $(CFLAGS) $< -o $@

src/shared.o: src/shared.cpp src/shared.h Makefile
	$(CC) $(CFLAGS) $< -o $@

src/krasner/krasner.o: src/krasner/krasner.cpp src/krasner/krasner.h src/krasner/krasnerExplicitTemplates.cpp src/planar_algebra/planar_algebra.h src/shared.h src/planar_algebra/sparsemat.h src/planar_algebra/coefficient_rings.h Makefile
	$(CC) $(CMPIRFLAGS) $(CFLAGS) $< -o $@

src/planar_algebra/sparsemat.o: src/planar_algebra/sparsemat.cpp src/planar_algebra/planar_algebra.h src/krasner/krasner.h src/shared.h src/planar_algebra/sparsemat.h src/planar_algebra/sparsematExplicitTemplates.cpp src/planar_algebra/coefficient_rings.h Makefile
	$(CC) $(CMPIRFLAGS) $(CFLAGS) $< -o $@
	
src/planar_algebra/coefficient_rings.o: src/planar_algebra/coefficient_rings.cpp src/planar_algebra/coefficient_rings.h src/shared.h src/planar_algebra/coefficient_rings_explicitTemplates.cpp Makefile
	$(CC) $(CMPIRFLAGS) $(CFLAGS) $< -o $@

src/planar_algebra/smith.o:  src/planar_algebra/smith.cpp src/planar_algebra/planar_algebra.h src/shared.h src/planar_algebra/explicitTemplates.cpp src/planar_algebra/sparsemat.h src/planar_algebra/coefficient_rings.h Makefile
	$(CC) $(CMPIRFLAGS) $(CFLAGS) $< -o $@

clean:
	rm -f $(OBJECTS) $(PYIO) $(PYXO) $(TESTER) $(EXECUTABLE) $(PYXCPP)
