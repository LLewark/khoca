#!/usr/bin/python3

#
#    KrasnerGaussToMyPDLib.py --- Part of khoca, a knot homology calculator
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

import sys, re

def inSublist(x, l):
	for i in range(len(l)):
		if x in l[i]:
			return i
	raise

def splitSublist(x, y, l):
	i = inSublist(x, l)
	idx1 = l[i].index(x)
	idx2 = l[i].index(y)
	minIdx = min(idx1, idx2)
	maxIdx = max(idx1, idx2)
	new1 = l[i][:minIdx] + l[i][maxIdx + 1:]
	new2 = l[i][minIdx + 1:maxIdx]
	l.pop(i)
	l.append(new1)
	l.append(new2)
	
def KrasnerGaussToMyPDmain(gauss):
	tmp = list(range(1, 2 * len(gauss) + 1))
	decompressed_gauss = [0] * 2 * len(gauss)
	signs = []
	for i in range(2 * len(gauss)):
		signs.append([True])
	mypd = []
	for i in gauss:
		x = tmp[0]
		decompressed_gauss[x - 1] = abs(i) - 1
		decompressed_gauss[abs(i) - 1] = x - 1
		if i < 0:
			signs[abs(i) - 1] = False
			signs[x - 1] = False
		tmp = tmp[1:]
		tmp.remove(abs(i))
	del tmp
	
	isInside = []
	for i in range(2 * len(gauss)):
		isInside.append([True, False])
	insidePart = [list(range(2 * len(gauss)))]
	outsidePart = [list(range(2 * len(gauss)))]
	idle = False
	
	while not idle:
		idle = True
		for i in range(2 * len(gauss)):
			if (len(isInside[i]) == 2):
				if inSublist(i, insidePart) != inSublist(decompressed_gauss[i], insidePart):
					isInside[i].remove(True)
					isInside[decompressed_gauss[i]].remove(True)
				if inSublist(i, outsidePart) != inSublist(decompressed_gauss[i], outsidePart):
					isInside[i].remove(False)
					isInside[decompressed_gauss[i]].remove(False)
				if isInside[i] == []:
					raise
				elif isInside[i] == [True]:
					splitSublist(i, decompressed_gauss[i], insidePart)
					idle = False
				elif isInside[i] == [False]:
					splitSublist(i, decompressed_gauss[i], outsidePart)
					idle = False
		if idle:
			for i in range(2 * len(gauss)):
				if (len(isInside[i]) == 2):
					isInside[i].remove(False)
					isInside[decompressed_gauss[i]].remove(False)
					splitSublist(i, decompressed_gauss[i], insidePart)
					idle = False
					break
	
	pd = []
	
	for i in range(len(decompressed_gauss)):
		x = decompressed_gauss[i]
		if i < x:
			toadd = [x, (x + 1) % (2 * len(gauss)), abs(i), (abs(i) + 1) % (2 * len(gauss))]
			toaddPd1 = [2 * x, 2 * abs(i) + 1, 2 * x + 1, (2 * abs(i) + 2) % (4 * len(gauss))]
			toaddPd2 = [2 * abs(i), 2 * x + 1, 2 * abs(i) + 1, (2 * x + 2) % (4 * len(gauss))]
	
			if isInside[i][0] and signs[i]:
				toaddPd1 = [toaddPd1[0], toaddPd1[3], toaddPd1[2], toaddPd1[1]]
				toaddPd2 = [toaddPd2[0], toaddPd2[3], toaddPd2[2], toaddPd2[1]]
			elif isInside[i][0] or signs[i]:
				toaddPd1 = [toaddPd1[1], toaddPd1[0], toaddPd1[3], toaddPd1[2]]
				toaddPd2 = [toaddPd2[1], toaddPd2[0], toaddPd2[3], toaddPd2[2]]
	
			if isInside[i][0]:
				toadd.reverse()
			if signs[i]:
				toadd = toadd[1:] + toadd[:1]
				toadd = [0] + toadd
			else:
				toadd = toadd[1:] + toadd[:1]
				toadd = [1] + toadd
			mypd.append(toadd)
			pd.append(toaddPd1)
			pd.append(toaddPd2)
	
	return mypd
