#!/usr/bin/python3

#
#    pseudoBraidToKrasnerGaussLib.py --- Part of khoca, a knot homology calculator
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


import sys

def pseudoBraidToKrasnerGaussMain(inputList):
	strands = max(inputList) + 1
	if (strands % 2 == 1):
		strands += 1

	direction = -1
	enum = [0 for i in range(len(inputList))]
	result = [0 for i in range(2*len(inputList))]
	num = 1
	for i in range(strands):
		idx = [j for j in range(len(inputList)) if (abs(inputList[j]) in [i, i+1])]
		if (direction == -1):
			idx.reverse()
		direction *= -1
		for j in idx:
			if (abs(inputList[j]) == i+1):
				enum[j] = num
			else:
				result[enum[j] - 1] = num * inputList[j] / abs(inputList[j])
			num += 1

	return [int(i) for i in result if (i != 0)]
