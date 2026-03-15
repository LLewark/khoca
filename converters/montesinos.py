#!/usr/bin/python3

#
#    montesinos.py --- Part of khoca, a knot homology calculator
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


from math import gcd
import sys, re

def sgn(x):
	assert(x != 0)
	return 1 if x > 0 else -1

# Increase last entry of a list x by d, return the result
def changeLast(x, d):
	assert(len(x) > 0)
	return [i for i in x[:-1]] + [x[len(x)-1] + d]

# Change sign of all entries of a list, return result
def neg(x):
	return [-i for i in x]

# Accepts as input two integers p, q, with q non-zero.
# Returns a continued fraction [a_1, ... a_n] = a_n + 1/(a_(n-1) + ...)
# for the rational p / q, where a_1, ... a_(n-1) are even numbers,
# and a_n is even if p or q is, and odd if p and q are.
def cf(p,q):
	if (p == 0): return [0]
	if (p < 0): return neg(cf(-p,q))
	if (p > q): return changeLast( cf( ((p + q) % (2*q)) - q, q), 2 * int((p + q) / (2*q)))
	if (((p * q) % 2) == 1): return changeLast( neg(cf(q - p, q)), 1)
	return changeLast(neg(cf(2*p - q, p)), 2) + [0]

def parse(s):                                 
	reg = re.compile("(-?[0-9]+)(/[0-9]+)?")
	x = reg.findall(sys.argv[1])
	return [[int(i[0]), 1 if (i[1] == "") else int(i[1][1:])] for i in x]

# Read from input
l = parse(sys.argv[1])
evenDen = None
j = 0
for i in l:
	assert(i[1] != 0)
	g = gcd(i[0],i[1])
	i[0] = int(i[0] / g)
	i[1] = int(i[1] / g)
	if (j == 0) and (len(l) == 1) and (l[0][0] % 2 == 1) and (l[0][1] % 2 == 1):
		l[0][1] += l[0][0]
	if (i[1] % 2 == 0) and (evenDen == None):
		evenDen = j
	j += 1
del j

L = [cf(l[i][0], l[i][1]) for i in range(len(l)) if i != evenDen]
s = sum([i[len(i) - 1] for i in L])
for i in L:
	i[len(i) - 1] = 0
if (evenDen != None):
	L.insert(evenDen, cf(l[evenDen][0] + l[evenDen][1] * s, l[evenDen][1]))
	s = 0
del l, evenDen
evens = [i for i in range(len(L)) if (len(L[i]) % 2 == 0)]
if not (((len(evens) == 1) and (s % 2 == 0)) or ((len(evens) == 0) and (s % 2 == 1))):
	print("This is not a knot.")
	sys.exit()
if (s % 2 != 0):
	print("Cannot produce a matched diagram.")
	sys.exit()
# cycle even to the front
L = L[evens[0]:] + L[:evens[0]]
# adapt sign (now: handedness, after: sign of crossing)
for i in L:
	for j in range(len(i)):
		if ((j % 2) == (len(i) % 2)):
			i[j] *= -1
                       
# part coming from s

T = [[[] for j in range(int(sum([abs(k) for k in i]) / 2))] for i in L]
T.insert(0, [[] for i in range(int(abs(s) / 2))])
    
#part coming from the even tangle

shift = 0
for i in range(int(abs(s) / 2)):
	T[0][i].append(sgn(s) * (shift + i + 1))
shift += int(abs(s) / 2)

minishift = 0
for i in range(int(len(L[0]) / 2)):
	for j in range(int(abs(L[0][2 * i]) / 2)):
		T[1][minishift + j].append(sgn(L[0][2 * i])*(j + shift + 1))
	shift += int(abs(L[0][2*i]) / 2)
	minishift += int(abs(L[0][2 * i + 1]) / 2) + int(abs(L[0][2*i]) / 2)

for i in range(int(abs(s) / 2)):
	T[0][int(abs(s) / 2) - i - 1].append(sgn(s) * (shift + i + 1))
shift += int(abs(s) / 2)

for i in range(len(L)):
	minishift = 0
	for j in range(len(L[len(L) - i - 1])):
		for k in range(int(abs(L[len(L) - i - 1][-j-1]) / 2)):
			T[len(L) - i][-k - minishift-1].append(sgn(L[len(L) - i - 1][-1-j])*(k + shift + 1))
		minishift += int(abs(L[len(L) - i - 1][-1-j]) / 2)
		shift += int(abs(L[len(L) - i - 1][-1-j]) / 2)
	minishift = 0
	for j in range(int(len(L[len(L) - i - 1]) / 2)):
		minishift += int(abs(L[len(L) - i - 1][2 * j]) / 2)
		for k in range(int(abs(L[len(L) - i - 1][2 * j + 1]) / 2)):
			T[len(L) - i][k + minishift].append(sgn(L[len(L) - i - 1][2 * j + 1])*(k + shift + 1))
		minishift += int(abs(L[len(L) - i - 1][2 * j + 1]) / 2)
		shift += int(abs(L[len(L) - i - 1][2 * j + 1]) / 2)

for i in range(len(L) - 1):
	minishift = 0
	for j in range(int(len(L[i + 1]) / 2)):
		for k in range(int(abs(L[i + 1][2 * j]) / 2)):
			T[i + 2][k + minishift].append(sgn(L[i + 1][2 * j])*(k + shift + 1))
		minishift += int(abs(L[i + 1][2 * j + 1]) / 2) + int(abs(L[i + 1][2*j]) / 2)
		shift += int(abs(L[i + 1][2 * j]) / 2)

T2 = [[min(abs(i[0]),abs(i[1])),max(i) if (i[0] > 0) else min(i)] for j in T for i in j]
T2.sort()
T3 = [i[1] for i in T2]
print(str(T3).replace(" ",""))
