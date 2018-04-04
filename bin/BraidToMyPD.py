#
#    braidToMyPD.py --- Part of khoca, a knot homology calculator
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

def chrToInt(i):
	t = ord(i)
	if (t in range(97, 123)):
		return t - 96
	elif (t in range(65, 91)):
		return 64 - t
	else:
		print("Unallowed char in braid: \"" + i + "\"")
		return 0

def BraidToMyPD(s):
	l = [chrToInt(i) for i in s]
	assert(not 0 in l)
	stumps = list(range(max(max(l), -min(l)) + 1))
	m = len(stumps)
	r = []
	for i in l:
		c = [3, m, m + 1, stumps[abs(i) - 1], stumps[abs(i)]]
		if i < 0:
			c = [2, c[3], c[4], c[1], c[2]]
		stumps[abs(i) - 1] = m + 1
		stumps[abs(i)] = m
		m += 2
		r.append(c)
	x = [item for subl in r for item in subl[1:]]
	for i in range(len(stumps)):
		if (not stumps[i] in x):
			print("At the moment, all strands of a braid must have at least one crossing.")
			assert(False)
		j = x.index(stumps[i])
		r[int(j / 4)][1 + j % 4] = i 
	return r	
