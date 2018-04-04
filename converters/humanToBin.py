#!/usr/bin/python3

#
#    humanToBin.py --- Part of khoca, a knot homology calculator
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



import sys, re, math

qShift_t = -16
tShift_t = qShift_t
dot_t = 16
vidx_t = 32
vidx2_t = 64
boundary_t = -8
uint32_t = 32

def myWrite(val, size):
	assert(size % 8 == 0)
	sys.stdout.buffer.write((val).to_bytes(int(abs(size) / 8), byteorder='little', signed=(size < 0)))

def myPop(size):
	res = nlist.pop(0)
	myWrite(res, size)
	return res

def mpirUInt(x):
	assert(x >= 0)
	if (x == 0):
		myWrite(0, 8)
	else:
		size = int(math.log(x, 256)) + 1
		myWrite(size, 8)
		myWrite(x, size * 8)
	

def mpirRat(coefficientRing):
	if (coefficientRing == 0):
		x = nlist.pop(0)
		if (x >= 0):
			myWrite(0, 8)
		else:
			myWrite(1, 8)
		myWrite(abs(x), 32)
	elif (coefficientRing == 1):
		num = nlist.pop(0)
		den = nlist.pop(0)
		if (num < 0):
			myWrite(1, 8)
		else:
			myWrite(0, 8)
		mpirUInt(abs(num))
		mpirUInt(den)
	elif (coefficientRing < 256):
		x = nlist.pop(0)
		myWrite(x, 8)
	else:
		x = nlist.pop(0)
		myWrite(x, 16)


def krasnerTangle(dataSize):
	myPop(qShift_t) # get qShift
	myPop(qShift_t) # get qShift times N
	myPop(boundary_t) # get Circle coun
	for i in range(dataSize):
		myPop(boundary_t)


def vecTangles(boundarySize):
	tanglesSize = myPop(vidx2_t)
	for i in range(tanglesSize):
		krasnerTangle(boundarySize)
	deloopStackSize = myPop(vidx2_t)
	for i in range(deloopStackSize):
		myPop(vidx2_t)

def krasnerCobo(coefficientRing):
	mpirRat(coefficientRing)
	
	size = myPop(boundary_t)
	for i in range(size):
		myPop(dot_t)

def LCCobos(coefficientRing):
	size = myPop(vidx2_t)
	for i in range(size):
		krasnerCobo(coefficientRing)

def matLCCobos(coefficientRing):
	colCount = myPop(uint32_t)
	entryCount = myPop(uint32_t)
	for i in range(entryCount):
		LCCobos(coefficientRing)
	for i in range(entryCount):
		myPop(uint32_t)
	rowPtrSize = myPop(uint32_t)
	for i in range(rowPtrSize):
		myPop(uint32_t)
	invertiblesSize = myPop(uint32_t)
	for i in range(invertiblesSize):
		myPop(uint32_t)

def complex():
	coefficientRing = myPop(tShift_t) # global T-shift
	myPop(qShift_t) # global T-shift
	boundarySize = myPop(boundary_t) # boundary size
	complexSize = myPop(vidx_t) # complex size ( = number of modules)
	for i in range(complexSize):
		vecTangles(boundarySize)
	for i in range(complexSize - 1):
		matLCCobos(coefficientRing)

nlist = [int(i) for i in re.findall(r'-?\d+', sys.stdin.read())];
complex()
