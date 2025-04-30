# cython: language_level=3
# distutils: language = c++

#
#   src/python_interface/pui.pyx --- Part of khoca, a knot homology calculator
#
# Copyright (C) 2018 Lukas Lewark <lukas@lewark.de>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

## @file pui.pyx
#  @brief Compiled into pui.cpp by cython, is called by khoroho.py, calls pythonInterface.cpp
#
#  Details 22.

import sys
from ctypes import *
from libcpp.vector cimport vector
from libcpp.deque cimport deque
from libcpp.string cimport string

cdef extern from "pythonInterface.h":
	cdef cppclass ComplexStack:
		ComplexStack(unsigned short mod, vector[int] F, int N, unsigned int girth, int verbose) except +
		int loadComplexFromFile(int idx, string fileName, int fileFormat)
		void deleteComplex(int idx)
		int tensorComplexes(int idx, int firstIdx, int secondIdx, int progress)
		void glueComplex(int idx, int gluePoint1, int gluePoint2)
		void deleteNonIsos(int idx)
		void printHomology(int idx)
		void calculateHomology(int idx, string &result)
		void calculateIntegralHomology(int idx, string &result, int progress)
		int simplifyComplexParallely(int idx, int numThreads, int progress)
		void saveComplexToFile(int idx, string fileName, int fileFormat)
		int copyComplex(int fromIdx, int toIdx)
		int dualizeComplex(int fromIdx, int toIdx)
		void setRoot(int r)
		void resetPage()
		void stepPage()
		void resetSimplificationsCounter(int idx)
		int getPage()
		int firstFreeIdx()
		void reducify(int idx)
		void printCompileInfo()
		void outputDetailed(int idx)

cdef ComplexStack* stack = NULL
cdef unsigned short globalMod

cdef class PyComplexStack(object):
	r"""
	Python wrapper for cpp class ComplexStack.
	"""
	cdef ComplexStack* former_stack

	def __enter__(self):
		r"""
		Save former global ``ComplexStack`` and reset current.
		"""
		global stack
		self.former_stack = stack
		stack = NULL

	def __exit__(self, exc_type, exc_val, exc_tb):
		r"""
		Delete current global ``ComplexStack`` and restore former stack.
		"""
		global stack
		del stack
		stack = self.former_stack


def pPrintCompileInfo():
	global stack
	stack.printCompileInfo()

def pPrintHomology(idx):
	global stack
	stack.printHomology(idx)

def pSetRoot(root):
	global stack
	stack.setRoot(root)

def pReducify(idx):
	global stack
	stack.reducify(idx)

def pDual(target, idx):
	global stack
	stack.dualizeComplex(idx, target)

def pDeleteComplex(idx):
	global stack
	stack.deleteComplex(idx)

def pSaveComplexToFile(idx, fileName):
	global stack
	stack.saveComplexToFile(idx, fileName, 0)

def pFirstFreeIdx():
	global stack
	return stack.firstFreeIdx()

def pDeleteNonIsos(idx):
	global stack
	stack.deleteNonIsos(idx)

def pGlueReduced(idx, gluePoint, num_threads, progress):
	global stack
	stack.glueComplex(idx, gluePoint, gluePoint - 1)
	stack.simplifyComplexParallely(idx, num_threads, progress)

def pIniStack(mod, frobenius, N, girth, verbose):
	global stack
	global globalMod
	if (stack != NULL):
		return
	stack = new ComplexStack(mod, frobenius, N, girth, verbose)
	globalMod = mod

def pLoadComplexFromFile(idx, s):
	global stack
	stack.loadComplexFromFile(idx, s, 0);

def pResetSimplificationsCounter(idx):
	global stack
	stack.resetSimplificationsCounter(idx)

def pCopyHomology(fromIdx, toIdx):
	global stack
	stack.copyComplex(fromIdx, toIdx)

def printFormattedHomology(s, printCommand, sortDegrees = True, printWidth = False):
	l = eval(s)
	l = [[i[2], i[0], i[1], i[3]] for i in l]
	# i[0]: torsion, i[1]: t-degree, i[1]: q-degree, i[3]: coefficient
	if sortDegrees:
		l.sort()
		m = [l[0]]
		for i in l[1:]:
			if (i[:-1] == m[-1][:-1]):
				m[-1][-1] += i[-1]
			else:
				m.append(i)
	else:
		m = l
	out = ""
	diagonalIndices = []
	for i in m:
		if (i[3] == 0):
			continue
		diagonalIndices.append(2 * i[1] + i[2])
		if (i[3] != 1):
			out += str(i[3])
		out += "t^" + str(i[1]) + "q^" + str(i[2])
		if (i[0] != 0):
			out += "[" + str(i[0]) + "]"
		out += " + "
	printCommand(out[:-3])
	if printWidth:
		printCommand("Width: " + str(1 + (max(diagonalIndices) - min(diagonalIndices)) // 2))

def pCalculateEquiHomology(idx, nice, num_threads, printCommand, printEquivariant, printWidth=False):
	global stack
	cdef string s
	# Over the integers, spectral sequences are not supported yet
	printCommand("Equivariant homology:")
	stack.calculateHomology(idx, s)
	printFormattedHomology(s, printCommand, False)
	stack.outputDetailed(idx)


def pCalculateHomology(idx, nice, num_threads, printCommand, printEquivariant, progress, printWidth=False):
	global stack
	cdef string s
	# Over the integers, spectral sequences are not supported yet
	if printEquivariant:
		printCommand("Equivariant homology:")
		stack.calculateHomology(idx, s)
		printFormattedHomology(s, printCommand, False)
		stack.outputDetailed(idx)
	else:
		printCommand("Non-equivariant homology:")
		if (globalMod == 0):
			stack.calculateIntegralHomology(idx, s, progress)
			printFormattedHomology(s, printCommand)
		else:
			hasBeenSimplified = 1
			printCommand("Page 1:")
			stack.calculateHomology(idx, s)
			printFormattedHomology(s, printCommand, printWidth)
			higherPageExists = False
			# 0: not finished
			# 1: only higher-degree isos left
			# 2: all matrices are zero
			# 4: only higher-degree isos left && nothing changed
			# 5: all matrices are zero && nothing changed
			while (hasBeenSimplified % 3 != 2):
				hasBeenSimplified = 0
				changes = False
				stack.stepPage()
				stack.resetSimplificationsCounter(idx)
				hasBeenSimplified = stack.simplifyComplexParallely(idx, num_threads, progress)
				if hasBeenSimplified < 3:
					higherPageExists = True
					printCommand("Page " + str(stack.getPage() + 1) + (" = infinity" if hasBeenSimplified % 3 == 2 else "") + ":")
					stack.calculateHomology(idx, s)
					printFormattedHomology(s, printCommand, printWidth)
			if (not higherPageExists):
				printCommand("The spectral sequence collapses on the first page.\n")
			stack.resetPage()
	return eval(s)

def pSum(target, idx1, idx2, num_threads, progress):
	global stack
	stack.tensorComplexes(target, idx1, idx2, progress)
	stack.glueComplex(target, 0, 2)
	hasBeenSimplified = 0
	stack.simplifyComplexParallely(target, num_threads, progress)
	
def pCalcSubTangleTree(s, num_threads, progress):
	global stack
	UNCOMPUTED = 0
	COMPUTED = 1
	COMPUTING = 2
        
	def calcSubTangleTree(idx):
		if (status[idx] == COMPUTED):
			return 0;
		status[idx] = COMPUTING
		if (s.types[idx] >= 0):
			loadType(idx + s.shift, s.types[idx])
			for i in s.selfGlue[idx]:
				stack.glueComplex(idx + s.shift, i[0], i[1])
				stack.simplifyComplexParallely(idx + s.shift, num_threads, progress)
		else:
			calcSubTangleTree(s.daughters[idx])
			calcSubTangleTree(s.sons[idx]) 
			stack.tensorComplexes(idx + s.shift, s.daughters[idx] + s.shift, s.sons[idx] + s.shift, progress)
			stack.deleteComplex(s.daughters[idx] + s.shift) # at the moment, we have a tree, so this is ok
			stack.deleteComplex(s.sons[idx] + s.shift)
			stack.resetSimplificationsCounter(idx + s.shift)
			for i in range(len(s.daughterGluePoints[idx])):
				stack.glueComplex(idx + s.shift, s.daughterGluePoints[idx][i], s.sonGluePoints[idx][i])
				stack.simplifyComplexParallely(idx + s.shift, num_threads, progress)

			# Delete after tensoring if this is "the last time" the subtangle was needed
		status[idx] = COMPUTED

	status = [UNCOMPUTED for i in range(len(s.daughters))]
	calcSubTangleTree(0)

def loadType(idx, typeIdx):
	global stack
	cdef string s
	data_dir = 'data'
	from os import path
	try:
		import khoca
		data_dir = path.join(path.dirname(khoca.__file__), data_dir)
	except ImportError:
		pass

	encode = 'utf-8'
	if (typeIdx == 0):
		s = bytes(path.join(data_dir, 'KrasnerPlus.bin'), encode)
	elif (typeIdx == 1):
		s = bytes(path.join(data_dir, 'KrasnerMinus.bin'), encode)
	elif (typeIdx == 2):
		s = bytes(path.join(data_dir, 'KhovanovPlus.bin'), encode)
	elif (typeIdx == 3):
		s = bytes(path.join(data_dir, 'KhovanovMinus.bin'), encode)
	else:
		print("Unknown type: " + str(typeIdx))
		sys.exit()
	stack.loadComplexFromFile(idx, s, 0)
