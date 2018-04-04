#!/usr/bin/python3

#
#    binToHuman.py --- Part of khoca, a knot homology calculator
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

# Reads from stdin, converts from the binary format to the human readable.

import sys

qShift_t = -16
tShift_t = qShift_t
dot_t = 16
vidx_t = 32
vidx2_t = 64
boundary_t = -8
uint32_t = 32

def myRead(size):
	assert(size % 8 == 0)
	return int.from_bytes(sys.stdin.buffer.read(int(abs(size / 8))), byteorder='little', signed=(size < 0))


def krasnerTangle(boundarySize, indent):
	print("  " * indent + "Tangle: Q-shift = " + str(myRead(qShift_t)) + " + " + str(myRead(qShift_t)) + "*N, " + str(myRead(boundary_t)) + " circles.")
	print("  " * (indent + 1) + "Data: " + " ".join([str(myRead(boundary_t)) for i in range(boundarySize)]))

def vecTangles(boundarySize, indent):
	tanglesSize = myRead(vidx2_t)
	print("  " * indent + "Module: Sum of " + str(tanglesSize) + " tangles.")
	for i in range(tanglesSize):
		krasnerTangle(boundarySize, indent + 2)
	deloopStackSize = myRead(vidx2_t)
	print("  " * (indent + 1) + "The deloop-stack has size " + str(deloopStackSize) + " and consists of the elements:" + " ".join([str(myRead(vidx2_t)) for i in range(deloopStackSize)]))

def mpirUInt():
	size = myRead(8)
	if (size == 0):
		return "0"
	else:
		return str(myRead(size * 8))

def mpirRat(coefficientRing):
	if (coefficientRing == 0):
		sgn = myRead(8)
		if (sgn == 0):
			return str(myRead(32))
		else:
			return str(-myRead(32))
	elif (coefficientRing == 1):
		sgn = myRead(8)
		num = mpirUInt()
		den = mpirUInt()
		if (sgn != 0):
			num = "-" + num
		return num + "/" + den
	elif (coefficientRing < 256):
		return str(myRead(8))
	else:
		return str(myRead(16))

def krasnerCobo(indent, coefficientRing):
	s = "  " * indent + mpirRat(coefficientRing) + " * cobordism with "
	size = myRead(boundary_t)
	s += str(size) + " facets and the following dot distribution: " + " ".join([str(myRead(dot_t)) for i in range(size)])
	print(s)

def LCCobos(indent, coefficientRing):
	size = myRead(vidx2_t)
	print("  " * indent + "LC of " + str(size) + " cobordisms:")
	for i in range(size):
		krasnerCobo(indent + 1, coefficientRing)

def matLCCobos(indent, coefficientRing):
	colCount = myRead(uint32_t)
	entryCount = myRead(uint32_t)
	print("  "*indent + "Map (as a matrix in yale format): " + str(colCount) + " columns, " + str(entryCount) + " entries.")
	print("  "*(indent + 1) + "Listing entries:")
	for i in range(entryCount):
		LCCobos(indent + 2, coefficientRing)
	print("  "*(indent + 1) + "Listing column indices:" + " ".join([str(myRead(uint32_t)) for i in range(entryCount)]))
	rowPtrSize = myRead(uint32_t)
	print("  "*(indent + 1) + "Number of rows plus one = " + str(rowPtrSize) + ", and the row pointers: " + " ".join([str(myRead(uint32_t)) for i in range(rowPtrSize)]))
	invertiblesSize = myRead(uint32_t)
	print("  "*(indent + 1) + "Number of invertible entries = " + str(invertiblesSize) + ", and their indices: " + " ".join([str(myRead(uint32_t)) for i in range(invertiblesSize)]))

def complex(indent):
	coefficientRing = myRead(tShift_t)
	print("Coefficient ring is " + str(coefficientRing))
	globalTShift = myRead(qShift_t)
	boundarySize = myRead(boundary_t)
	size = myRead(vidx_t)
	print("  "*indent + "Chain complex: global t-shift = " + str(globalTShift) + ", " + str(boundarySize) + " boundary points, " + str(size) + " modules.")
	print("  "*(indent + 1) + "Listing modules:")
	for i in range(size):
		vecTangles(boundarySize, indent + 2)
	print("  "*(indent + 1) + "Listing maps:")
	for i in range(size - 1):
		matLCCobos(indent + 2, coefficientRing)

complex(0)
