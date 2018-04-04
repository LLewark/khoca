#!/usr/bin/python3

#
#    KrasnerGaussToMyPD.py --- Part of khoca, a knot homology calculator
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

# 2n points on a circle are enumerated from 1 to 2n.
# Each of the points is paired with another one,
# and the connection may be positive or negative.
# This is encoded as a list of integers of even length:
# The element at index i specifies to which point the
# point i is connected; its sign is the sign of the connection.
# To compress, drop half of the list items, namely for each
# connection the higher number. The encoding remains injective.

import sys, re
sys.path.append('.')
import KrasnerGaussToMyPDLib

if len(sys.argv) != 2:
	print("Expecting one argument.")
	sys.exit()
if sys.argv[1] == "-":
	x = sys.stdin.read()
else:
	x = sys.argv[1]
gauss = [int(j) for j in re.findall("-?[0-9]+", x)]
print(str(KrasnerGaussToMyPDLib.KrasnerGaussToMyPDmain(gauss)).replace(" ",""))
