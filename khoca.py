#!/usr/bin/python3
## @file

#
#    khoca.py --- Part of khoca, a knot homology calculator
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

# The file to be run by the user.
# The main work is done by pui.pyx, which is called by this file.

import sys
sys.path.insert(0, './bin/')

import os

def interactive():
    r"""
    Checks if the module is used as part of :class:`InteractiveCalculator`.
    """
    head, tail = os.path.split(__file__)
    if not head:
        return False
    head, tail = os.path.split(head)
    if tail == 'khoca':
        return True
    else:
        return False

def run_command_exit(print_command, err_mess):
    r"""
    Exit :func:`run_command_line` on error.
    """
    if interactive():
        raise ValueError(err_mess)
    else:
        print_command(err_mess)
        sys.exit()

def set_options(verbose_mode, progress_mode):
    r"""
    Set the global options ``verbose`` and ``progress``.
    """
    global verbose, progress
    verbose = verbose_mode
    progress = progress_mode

import random
import re
import math

sys.path.append('.')
if interactive():
    from khoca.bin.pui import pCalcSubTangleTree, pCalculateHomology, pGlueReduced, pCopyHomology, pFirstFreeIdx, pSaveComplexToFile, pDeleteComplex, pLoadComplexFromFile, pSum, pDual, pDeleteNonIsos, pResetSimplificationsCounter, pReducify, pIniStack, pSetRoot, pPrintHomology, pPrintCompileInfo, pCalculateEquiHomology
    from khoca.bin import KrasnerGaussToMyPDLib, pseudoBraidToKrasnerGaussLib, BraidToMyPD
else:
    from pui import pCalcSubTangleTree, pCalculateHomology, pGlueReduced, pCopyHomology, pFirstFreeIdx, pSaveComplexToFile, pDeleteComplex, pLoadComplexFromFile, pSum, pDual, pDeleteNonIsos, pResetSimplificationsCounter, pReducify, pIniStack, pSetRoot, pPrintHomology, pPrintCompileInfo, pCalculateEquiHomology
    import KrasnerGaussToMyPDLib, pseudoBraidToKrasnerGaussLib, BraidToMyPD


debugging = False
NUM_THREADS = 1 if debugging else 12


## Converts a link diagram given in the Planar-diagram notation to the mypd format.
# @param pd A list of lists of four integers each, such as [[0,1,2,3],[1,4,5,2],[4,0,3,5] for the positive trefoil.
# @return That link diagram in mypd-notation, type "0" for positive, type "1" for negative crossings.
def pd_to_mypd(pd):
    result = list()
    for i in pd:
        if (i[3] == ((i[1] + 1) % (2 * len(pd)))):
            result.append([2,i[0] - 1,i[1] - 1,i[2] - 1,i[3] - 1])
        else:
            result.append([3,i[1] - 1,i[2] - 1,i[3] - 1,i[0] - 1])
    return result

def generate_binary_tree(length, s = 0):
    if (length == 1):
        return s
    return [generate_binary_tree(length / 2, s), generate_binary_tree((length / 2) + (length % 2), s + (length / 2))]

def generate_standard_tree(l):
    if (len(l) == 1):
        return l[0]
    return [generate_standard_tree(l[:-1]), l[len(l) - 1]]

def isBetterSignature(s1, s2):
    s1Copy = list(s1)
    s2Copy = list(s2)
    s1Copy.sort(reverse = True)
    s2Copy.sort(reverse = True)
    if (s1Copy < s2Copy): return True
    if (s1Copy > s2Copy): return False
    return (s1 > s2)

# Using a randomised algorithm, produces a linear tree such that intermediate tangles have small girth
def generate_girth_optimised_tree(pd):
    pd = [[i] + pd[i][1:] for i in range(len(pd))]
    iterations = 1000

    bestOrder = []
    bestSignature = []
    girth = 0
    for k in range(iterations):
        copyPd = list(pd)
        openEnds = set()
        order = []
        signature = []
        while len(copyPd) > 0:
            girthChange = [4 - 2 * len(openEnds & set(i[1:])) for i in copyPd]
            bestGirthChange = min(girthChange)
            bestCrossings = [i for i in range(len(girthChange)) if girthChange[i] == bestGirthChange]
            newCrossing = bestCrossings[0] if debugging else random.choice(bestCrossings)
            order.append(copyPd[newCrossing][0])
            openEnds = openEnds ^ set(copyPd.pop(newCrossing)[1:])
            signature.append(len(openEnds))
        if (bestOrder == []) or (isBetterSignature(signature, bestSignature)):
            bestOrder = order
            bestSignature = signature
            girth = max(bestSignature)
    return (girth, generate_standard_tree(bestOrder))

def eliminatedoubles(l):
    i = 0
    while (i < len(l)):
        j = i + 1
        found = False
        while (j < len(l)):
            if l[i] == l[j]:
                l.pop(j)
                l.pop(i)
                found = True
                break
            j = j + 1
        if not found:
            i = i + 1

## Stores the information which elementary tangles to load and how and in which order to glue them.
# In particular, all index shifts have been done, so that calcSubTangleTree may
# call functions in pythonInterface.cpp without further ado.
# Details: This represents a directed tree with vertices numbered 0,1,....
# Each vertex of number i has at most two children, whose number is saved in
# sons[i] and daughers[i] (-1 meaning no child). types saves whether the tangle at the vertex
# is elementary and should be loaded from a file (then types[i] is 0,1,2,3, see function loadType in pui.pyx),
# or the tangle is obtained by gluing son and daughter (then types[i] is -1).
# boundaryPoints[i] is an orderd list of the tangles' boundary points (in global numbering).
# when gluing e.g. a daughter and son with boundaryPoints 2,3,0,1 and 4,6,3,2, respectively, one
# first takes disjoint union (giving a tangle with boundaryPoints 2,3,0,1,4,6,3,2), and then
# glues the matching boundaryPoints; which match which is saved in daughterGluePoints and sonGluePoints,
# which are 0,0 and 7,5 here (because index shifts are already contained).
class CSchedule:
    def initialise(self, mypd):
        nodes = 2 * len(mypd) - 1
        self.shift = 0
        self.daughters = [-1] * nodes
        self.sons = [-1] * nodes
        self.types = [-1] * nodes
        self.selfGlue = [list() for i in range(nodes)]
        self.boundaryPoints = [list() for i in range(nodes)]
        self.sonGluePoints = [list() for i in range(nodes)]
        self.daughterGluePoints = [list() for i in range(nodes)]
        return nodes + 1

    def output(self):
        sys.stderr.write(str(self.daughters) + str(self.sons) + str(self.types) + str(self.boundaryPoints) + str(self.daughterGluePoints) + str(self.sonGluePoints) + str(self.selfGlue) + "\n")


    def getFromMypdAndTree(self, mypd, tree, idx):
        if (type(tree) == type(int())):
            self.types[idx] = mypd[tree][0]
            self.boundaryPoints[idx] = mypd[tree][1:]
            unitedBoundaryPoints = list(self.boundaryPoints[idx])
            i = 0
            while (i < len(unitedBoundaryPoints)):
                j = i + 1
                found = False
                while (j < len(unitedBoundaryPoints)):
                    if unitedBoundaryPoints[i] == unitedBoundaryPoints[j]:
                        self.selfGlue[idx].append([i,j])
                        unitedBoundaryPoints.pop(j)
                        unitedBoundaryPoints.pop(i)
                        found = True
                        break
                    j = j + 1
                if not found:
                    i = i + 1
            return idx + 1
        else:
            self.daughters[idx] = idx + 1
            self.sons[idx] = self.getFromMypdAndTree(mypd, tree[0], idx + 1)
            result = self.getFromMypdAndTree(mypd, tree[1], self.sons[idx])
            beforeDeletion = (self.boundaryPoints[self.daughters[idx]] + self.boundaryPoints[self.sons[idx]])
            self.boundaryPoints[idx] = [ x for x in beforeDeletion if beforeDeletion.count(x) == 1]
            shift = len(self.boundaryPoints[self.daughters[idx]])
            daughterBP = list(self.boundaryPoints[self.daughters[idx]])
            eliminatedoubles(daughterBP)
            sonBP = list(self.boundaryPoints[self.sons[idx]])
            eliminatedoubles(sonBP)
            unitedBoundaryPoints = list(daughterBP)
            unitedBoundaryPoints.extend(sonBP)
            i = 0
            while (i < len(unitedBoundaryPoints)):
                j = i + 1
                found = False
                while (j < len(unitedBoundaryPoints)):
                    if unitedBoundaryPoints[i] == unitedBoundaryPoints[j]:
                        self.daughterGluePoints[idx].append(i)
                        self.sonGluePoints[idx].append(j)
                        unitedBoundaryPoints.pop(j)
                        unitedBoundaryPoints.pop(i)
                        found = True
                        break
                    j = j + 1
                if not found:
                    i = i + 1
            return result

def reducePD(mypd):
    x = [item for subl in mypd for item in subl[1:]]
    j = x.index(max(x))
    mypd[int(j / 4)][1 + j % 4] +=1
    return mypd

def calcIt(mypd, shift):
    global verbose, progress
    (girth, tree) = generate_girth_optimised_tree(mypd)
    pIniStack(stackMod, stackFrobenius, N, girth, verbose)
    pSetRoot(stackRoot)

    s = CSchedule()
    result = s.initialise(mypd)
    s.getFromMypdAndTree(mypd, tree, 0)
    s.shift = shift

    pCalcSubTangleTree(s, NUM_THREADS, progress)
    return result

def getInts(s, signed):
    return [int(j) for j in re.findall( ("-?" if signed else "") + "[0-9]+", s[3:])]

## @var mypd
# The internal format for a closed object (as a link), glued together from smaller tangles (as tangles).
#
# mypd is a list of tangles; an object is a list of integers.
# The first of that integers is the type of the object (e.g. a positive or a
# negative crossing), all following integers are the numbers of the boundary
# points of that object, the order of which matters.
# %Boundary point numbers need not be consecutive; each boundary point number
# must appear exactly twice in the whole mypd.

def isAllowedStackMod(a):
    if not a.isdigit():
        return False
    x = int(a)
    if (x <= 3):
        return True
    for i in range(2, int(math.sqrt(x) + 2)):
        if (x % i == 0):
            return False
    return True

def run_commandline(argv, printCommand, progress):
        NUM_ARGUMENTS = 3
        global stackMod, stackFrobenius, N, stackRoot

        # Parsing the arguments
        if (len(argv) <= NUM_ARGUMENTS):
            printCommand(HELP_TEXT)
            return 1

        if isAllowedStackMod(argv[1]):
            stackMod = int(argv[1])
        else:
            printCommand("First argument must be 0, 1 or a prime.")
            return 1

        N = 0
        if argv[2][0] == "e":
            N = int(argv[2][1:])
            l = []
            if (N < 2):
                printCommand("e must be followed by an integer that is at least 2.")
                return 1
        else:
            l = [int(i) for i in re.findall("-?[0-9]+", argv[2])]
            if (len(l) == 0):
                printCommand("Second argument must be a list of at least one integer, or e (for equivariant).")
                return 1
        stackFrobenius = l
        equivariant = (len(stackFrobenius) == 0)

        if argv[3].isdigit() or ((argv[3][0] in ["+","-"]) and (argv[3][1:].isdigit())):
            stackRoot = int(argv[3])
        else:
            printCommand("Third argument must be a signed integer.")
            return 1

        # doing the work
        shift = 0
        idxTranslator = [0]
        res = []
        for i in argv[(NUM_ARGUMENTS + 1):]:
            if i[:3].capitalize() == "Nat":
                l = getInts(i, False)
                if (len(l) % 5 != 0) or (len(l) == 0):
                    run_command_exit(printCommand, "The native code (\"" + i[3:] +  "\") following the command \"nat\" must be a non-empty list of integers whose length is divisible by five.")

                mypd = [list(j) for j in zip(*[iter(l)]*5)]
                mypd = reducePD(mypd)
                shift += calcIt(mypd, shift)
                idxTranslator += [shift]
            if i[:7].capitalize() == "Special":
                l = eval(i[7:])
                shift += calcIt(l, shift)
                idxTranslator += [shift]
            elif i[:5].capitalize() == "Gauss":
                l = getInts(i, True)
                mypd = reducePD(KrasnerGaussToMyPDLib.KrasnerGaussToMyPDmain(l))
                shift += calcIt(mypd, shift)
                idxTranslator += [shift]
            elif i[:11].capitalize() == "PseudoBraid":
                m = getInts(i, True)
                l = pseudoBraidToKrasnerGaussLib.pseudoBraidToKrasnerGaussMain(m)
                mypd = reducePD(KrasnerGaussToMyPDLib.KrasnerGaussToMyPDmain(l))
                shift += calcIt(mypd, shift)
                idxTranslator += [shift]
            elif i[:5].capitalize() == "Braid":
                mypd = reducePD(BraidToMyPD.BraidToMyPD(i[5:]))
                shift += calcIt(mypd, shift)
                idxTranslator += [shift]
            elif i[:2].capitalize() == "Pd":
                l = getInts(i, False)
                l = [x - 1 for x in l]
                mypd = reducePD(pd_to_mypd([list(j) for j in zip(*[iter(l)]*4)]))
                print(mypd)
                shift += calcIt(mypd, shift)
                idxTranslator += [shift]
            elif i[:10].capitalize() == "Calcnonred":
                l = getInts(i, False)
                if (len(l) != 1):
                    run_command_exit(printCommand, "\"Calcnonred\" must be followed by exactly one unsigned integer.")
                idx = idxTranslator[l[0]]
                firstFreeIdx = shift
                shift += 1
                idxTranslator[-1] += 1
                pResetSimplificationsCounter(idx)
                pCopyHomology(idx, firstFreeIdx)
                printCommand("Result:")
                printCommand("Unreduced Homology:")
                pGlueReduced(firstFreeIdx, 1, NUM_THREADS, progress)
                s = pCalculateHomology(firstFreeIdx, False, NUM_THREADS, printCommand, equivariant, progress)
                res.append(s)
                pDeleteComplex(firstFreeIdx)
            elif i[:4].capitalize() == "Calc":
                nice = (i[4:8].capitalize() == "Nice")
                l = getInts(i, False)
                if (len(l) != 1):
                    run_command_exit(printCommand, "\"Calc\" must be followed by exactly one unsigned integer.")
                idx = idxTranslator[l[0]]
                firstFreeIdx = shift
                secondFreeIdx = shift + 1
                shift += 2
                idxTranslator[-1] += 2
                pResetSimplificationsCounter(idx)
                pCopyHomology(idx, firstFreeIdx)
                pCopyHomology(idx, secondFreeIdx)
                printCommand("Result:")
                printCommand("Reduced Homology:")
                pReducify(secondFreeIdx)
                s = pCalculateHomology(secondFreeIdx, nice, NUM_THREADS, printCommand, equivariant, progress, True)
                res.append(s)
                printCommand("Unreduced Homology:")
                pGlueReduced(firstFreeIdx, 1, NUM_THREADS, progress)
                s = pCalculateHomology(firstFreeIdx, nice, NUM_THREADS, printCommand, equivariant, progress, True)
                res.append(s)
                pDeleteComplex(firstFreeIdx)
                pDeleteComplex(secondFreeIdx)
            elif i[:4].capitalize() == "Save":
                l = re.match("[0-9]+", i[4:])
                if (l == None):
                    run_command_exit(printCommand, "\"Save\" must be followed by an unsigned integer and a filename.")
                pSaveComplexToFile(idxTranslator[int(l.group(0))], i[(4 + len(l.group(0))):])
            elif i[:4].capitalize() == "Load":
                pLoadComplexFromFile(shift, i[4:])
                shift += 1
                idxTranslator += [shift]
            elif i[:4].capitalize() == "Dual":
                l = re.match("[0-9]+", i[4:])
                if (l == None):
                    run_command_exit(printCommand, "\"Dual\" must be followed by an unsigned integer.")
                pDual(shift, int(l.group(0)))
                shift += 1
                idxTranslator += [shift]
            elif i[:3].capitalize() == "Sum":
                l = getInts(i, False)
                if (len(l) != 2):
                    run_command_exit(printCommand, "\"Sum\" must be followed by exactly two unsigned integer.")
                pSum(shift, idxTranslator[l[0]], idxTranslator[l[1]], NUM_THREADS, progress)
                shift += 1
                idxTranslator += [shift]
            else:
                run_command_exit(printCommand, str(i) + " is not a valid command.")
        return res


def parseOptions(argv):
    global verbose, progress
    o = [x for x in argv[1:] if (x[0] == '-')]
    for i in o:
        if (i == '-v'):
            verbose = 1
        elif (i == '-p'):
            progress = 1
        elif (i == '-h'):
            print(HELP_TEXT)
            sys.exit()
        else:
            print("Unknown option " + i + ".")
            sys.exit()
    return [argv[0]] + [x for x in argv[1:] if (x[0] != '-')]

HELP_TEXT = "Expecting at least three arguments:\n(1) The coefficient ring; 0 for integers, 1 for rationals, a prime p for the finite field with p elements,\n(2) the vector [a_0, ... a_{N-1}] defining the Frobenius algebra F[X]/(X^N+a_{N-1}X^{N-1}+...+a_0) or e followed by the number N for equivariant,\n(3) a root of the polynomial given in (2).\nAny following argument is interpreted as command. For example,\n./khoca.py 0 0.0 0 braidaBaB calc0\ncomputes integral Khovanov sl(2)-homology of the figure-eight knot.\nRefer to the README file or to http://lewark.de/lukas/khoca.html for more detailed help."
verbose = 0
progress = 0

if (verbose):
    pPrintCompileInfo()
    sys.stderr.write(("Multithreading with " + str(NUM_THREADS) + " threads.\n") if (NUM_THREADS > 1) else "No multithreading.\n")

if not interactive():
    a = parseOptions(sys.argv)
    run_commandline(a, print, progress)
