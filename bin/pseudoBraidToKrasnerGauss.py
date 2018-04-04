#!/usr/bin/python3

#
#    pseudoBraidToKrasnerGauss.py --- Part of khoca, a knot homology calculator
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

# Accepted input: A list of signed integers, where i stands for
# (sigma_i)^2 and -i for (sigma_i)^(-2). The braid is taken to
# have an even number of strands, which are then connected as follows:
# On the left, 1 to 2, 3 to 4 and so on; on the right, 1 to the last, 2 to 3, 4 to 5 and so on.
# Output is a list of signed integers representing the Krasner Gauss
# encoding of knot (so the output may be fed to KrasnerGaussToMyPD.py).

import sys, re
sys.path.append('.')
import pseudoBraidToKrasnerGaussLib

if len(sys.argv) != 2:
	print("Expecting one argument.")
	sys.exit()
if sys.argv[1] == "-":
	x = sys.stdin.read()
else:
	x = sys.argv[1]
braid = [int(j) for j in re.findall("-?[0-9]+", x)]
print(str(pseudoBraidToKrasnerGaussLib.pseudoBraidToKrasnerGaussMain(braid)).replace(" ",""))
