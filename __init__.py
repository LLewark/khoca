#!/usr/bin/python3
## @file

#
#    khoca/__init__.py --- Part of khoca, a knot homology calculator
#
# Copyright (C) 2021 Sebastian Oehms <seb.oehms@gmail.com>
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

# This file contains a Python wrapper class of the command line version
# of Khoca

# The main work is done by pui.pyx via run_commandline, which is called
# by the instance call of this class.

from khoca.khoca import run_commandline
from khoca.bin.pui import PyComplexStack


class InteractiveCalculator:
    r"""
    Class to allow the usage of `Khoca` interactively in a Python session.

    EXAMPLES::

        >>> from khoca import InteractiveCalculator
        >>> KH = InteractiveCalculator()
        >>> KH
        Khovanov homology calculator for Frobenius algebra: Z[X] / (1*X^2).
        >>> KH('braidaBaB')    # doctest: +NORMALIZE_WHITESPACE
        [[[-2, 4, 0, 1], [-1, 2, 0, 1], [0, 0, 0, 1], [1, -2, 0, 1], [2, -4, 0, 1],
          [-2, 4, 0, 0], [-1, 4, 0, 0], [-1, 2, 0, 0], [0, 2, 0, 0], [0, 0, 0, 0],
          [1, 0, 0, 0], [1, -2, 0, 0], [2, -2, 0, 0]], [[-2, 3, 0, 1], [-2, 5, 0, 1],
          [-1, 1, 0, 1], [-1, 3, 0, 1], [0, -1, 0, 1], [0, 1, 0, 1], [1, -3, 0, 1],
          [1, -1, 0, 1], [2, -5, 0, 1], [2, -3, 0, 1], [-1, 3, 2, 1], [-2, 3, 0, -1],
          [-1, 3, 0, -1], [-2, 5, 0, 0], [-1, 5, 0, 0], [-1, 1, 0, 0], [0, 1, 0, 0],
          [-1, 3, 0, 0], [0, 3, 0, 0], [0, -1, 0, 0], [1, -1, 0, 0], [0, 1, 0, 0],
          [1, 1, 0, 0], [2, -3, 2, 1], [1, -3, 0, -1], [2, -3, 0, -1], [1, -1, 0, 0],
          [2, -1, 0, 0]]]
        >>> res, mess = KH('braidaaa', print_messages=True); mess # doctest: +NORMALIZE_WHITESPACE
        ['Result:', 'Reduced Homology:', 'Non-equivariant homology:',
         't^-3q^8 + t^-2q^6 + t^0q^2', 'Unreduced Homology:',
         'Non-equivariant homology:', 't^-3q^9 + t^-2q^5 + t^0q^1 + t^0q^3 + t^-2q^7[2]']
        >>> KH((1, 1, 1)) == KH('braidaaa')
        True
        >>> KH([[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]) == KH('braidaBaB')
        True
        >>> KH((1, -2, 1, -2)) == KH('braidaBaB')
        True
        >>>
    """
    def __init__(self, coefficient_ring=0, frobenius_algebra=(0, 0), root=0, equivariant=None):
        r"""
        Constructor.

        INPUT:

            - ``coefficient_ring`` -- coefficient ring ``F``of the homology. It
              can be given by an integer (or string convertible to an integer):
              -  ``0`` (default) the ring of integers
              -  ``1``  the rational field
              -  a prime for the corresponding finite field
            - ``frobenius_algebra`` -- the Frobenius algebra ``F[x]/p`` given by
              the coefficients of the normed polynomial ``p`` (default is
              ``p = X^2``) as a tuple of length ``deg(p)`` where the constant
              term is the first entry.
            - ``root`` -- a root of the polynomial ``p`` (default is zero).
            - ``equivariant`` optional integer ``n > 1`` giving the degree of
              ``sl(n)`` for equivariant homology. If this is given then the input
              to the keyword arguments  ``frobenius_algebra`` and ``root`` will
              be ignored. If it is not given then non equivariant homology will
              be calculated.

        EXAMPLES::

            >>> from khoca import InteractiveCalculator
            >>> InteractiveCalculator(1, (0, 1), 0)
            Khovanov homology calculator for Frobenius algebra: Q[X] / (1*X^2 + 1*X).
            >>> from khoca import InteractiveCalculator
            >>> KH = InteractiveCalculator(equivariant=3); KH
            Khovanov homology calculator for Frobenius algebra: Z[a, b, c][X] / (1*X^3 + c*X^2 + b*X + a).
        """
        if type(coefficient_ring) == str:
            self._coefficient_ring = coefficient_ring
        else:
            try:
                self._coefficient_ring = str(int(coefficient_ring))
            except TypeError:
                raise TypeError('Coefficient ring must be declared by an integer or string')
        if type(frobenius_algebra) == str:
            self._frobenius_algebra = frobenius_algebra
        else:
            try:
                self._frobenius_algebra = '%s' %list(frobenius_algebra)
            except TypeError:
                raise TypeError('Frobenius algebra must be declared by a tuple, list or string')
        if type(root) == str:
            self._root = root
        else:
            try:
                self._root = str(int(root))
            except TypeError:
                raise TypeError('Root must be declared by an integer or string')

        self._equivariant = None
        if equivariant:
            try:
                 equivariant = int(equivariant)
            except TypeError:
                raise TypeError('equivariant must be given as integer')
            if equivariant <= 1:
                raise ValueError('equivariant must be larger than 1')
            self._equivariant = 'e%s' %equivariant

        self._description = None
        self._stack = PyComplexStack()
        self('braidaA', ' ') # sets self._description

    def __repr__(self):
        r"""
        String representation of ``self``.

        EXAMPLES::

            >>> from khoca import InteractiveCalculator
            >>> InteractiveCalculator(5, (1, 0), 0)
            Khovanov homology calculator for Frobenius algebra: F_5[X] / (1*X^2 + 1).
        """
        return 'Khovanov homology calculator for %s' %self._description


    def __call__(self, link, command=None, print_messages=False, verbose=False, progress=False):
        r"""
        Instance call to apply the calculator to a link with respect to a
        certain command.

        INPUT:

            - ``link`` -- the link to which the calculation should be done.
              It can be given as a Tuple which will be interpreted as a braid
              in Tietze form which closure is the link in question. Further
              a list of lists is accepted which will be interpreted as a list
              of crossings in pd-notation. Alternatively you can declare the
              link by a string. The following formats are accepted:

              - ``BraidX`` for a braid formatted as in knotscape (``a`` = first
                Artin generator, ``A`` = its inverse, ``b`` = second Artin
                generator, etc.). This works only for ``sl(2)`` homology,
                otherwise output is nonsensical
              - ``PdX`` for input in PD notation (as e.g. given on !KnotInfo).
                Again, this works only for ``sl(2)`` homology, otherwise output
                is nonsensical
              - ``GaussX`` for bipartite knot given as a matched diagram,
                following the convention explained in the section below. This
                works for ``sl(N)`` homology for all ``N``

            - ``command`` -- a command given as string as explained in the
              command line version of the functionality. This can be:

              - ``MirrorX`` takes the dual of the result at spot ``X``
              - ``SumXY`` computes the homology of the connected sum of the
                results saved at spots ``X`` and ``Y`` (numbers separated by
                a non-digit character)
              - ``CalcX` outputs the result saved at spot X. If you forget
                this command, the program will have no output

              If this keyword argument is not given it will be interpreted
              as ``Calc0`` by default.

              - ``print_messages`` boolean (default is ``False``). If set to
                ``True`` the print output of the command line version is
                returned as a list on the second output position. By default
                all print messages to ``stdout`` of the command line version
                are suppressed
              - ``verbose`` boolean (default is ``False``). If it is set to
                ``True`` all print messages to ``stdout`` together with special
                verbose messages of the command line version are printed
              - ``progress`` boolean (default is ``False``). If it is set to
                ``True`` progress bars will be printed

        OUTPUT:

            Two lists of quadruples the first for reduced and the second for
            unreduced homology. In such a quadruple the first item stands for
            the ``t``-degree, the second for the ``q``-degree, the third for
            the torsion and the last item stands for the coefficient of the
            corresponding summand of the homology.

            If the option ``print_messages`` is given, than the output contains
            a second result which is a list of all suppressed print messages.


        ..NOTE:

            The program keeps a stack of computed homologies, enumerated 0,1,2... .
            Each link declaration and each command except ``CalcX`` puts a new
            homology on that stack, whereas ``CalcX`` prints the homology at a
            certain spot. This is mainly useful to compute homology of sums of
            knots. As an example of such a stack operations see the first example
            below.

        EXAMPLES::

            >>> from khoca import InteractiveCalculator
            >>> KH = InteractiveCalculator(1, (1, 0), 0); KH
            Khovanov homology calculator for Frobenius algebra: Q[X] / (1*X^2 + 1).
            >>> KH('braidaaa', 'dual0 sum0+1 braidaBaB sum2+3 calc4', print_messages=True) # doctest: +NORMALIZE_WHITESPACE
            ([[[-5, 10, 0, 1], [-4, 8, 0, 1], [-4, 8, 0, 1], [-3, 6, 0, 1],
               [-3, 6, 0, 1], [-3, 6, 0, 1], [-2, 4, 0, 1], [-2, 4, 0, 1],
               [-2, 4, 0, 1], [-2, 4, 0, 1], [-2, 4, 0, 1], [-2, 4, 0, 1],
               [-1, 2, 0, 1], [-1, 2, 0, 1], [-1, 2, 0, 1], [-1, 2, 0, 1],
               [-1, 2, 0, 1], [-1, 2, 0, 1], [-1, 2, 0, 1], [0, 0, 0, 1],
               [0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 0, 1], [0, 0, 0, 1],
               [0, 0, 0, 1], [0, 0, 0, 1], [1, -2, 0, 1], [1, -2, 0, 1],
               [1, -2, 0, 1], [1, -2, 0, 1], [1, -2, 0, 1], [1, -2, 0, 1],
               [1, -2, 0, 1], [2, -4, 0, 1], [2, -4, 0, 1], [2, -4, 0, 1],
               [2, -4, 0, 1], [2, -4, 0, 1], [2, -4, 0, 1], [3, -6, 0, 1],
               [3, -6, 0, 1], [3, -6, 0, 1], [4, -8, 0, 1], [4, -8, 0, 1],
               [5, -10, 0, 1]], [[0, -1, 0, 1], [0, 1, 0, 1]]],
             ['Result:', 'Reduced Homology:', 'Non-equivariant homology:',
              'Page 1:', 't^-5q^10 + 2t^-4q^8 + 3t^-3q^6 + 6t^-2q^4 + 7t^-1q^2 + 7t^0q^0
                       + 7t^1q^-2 + 6t^2q^-4 + 3t^3q^-6 + 2t^4q^-8 + t^5q^-10',
              'The spectral sequence collapses on the first page.\n',
              'Unreduced Homology:', 'Non-equivariant homology:',
              'Page 1:', 't^-5q^11 + t^-4q^7 + t^-4q^9 + t^-3q^5 + 2t^-3q^7 + 2t^-2q^3
                       + 4t^-2q^5 + 4t^-1q^1 + 3t^-1q^3 + 4t^0q^-1 + 4t^0q^1 + 3t^1q^-3
                       + 4t^1q^-1 + 4t^2q^-5 + 2t^2q^-3 + 2t^3q^-7 + t^3q^-5 + t^4q^-9
                       + t^4q^-7 + t^5q^-11',
              'Page 3 = infinity:', 't^0q^-1 + t^0q^1'])
        """
        # Todo: allow Python type inputs like tuples ..
        from .khoca import set_options
        verbose = int(verbose); progress = int(progress)
        set_options(verbose, progress)
        arg_list = [None]
        arg_list.append(self._coefficient_ring)
        if self._equivariant:
            arg_list.append(self._equivariant)
        else:
            arg_list.append(self._frobenius_algebra)
        arg_list.append(self._root)


        if type(link) in (tuple, list):
            # Allow more Python like ways to declare the link
            try:
                link_list = [int(j) for j in link]
                # interpreted as braid in Tietze-form
                offset = {-1:64, 1:96}
                from math import copysign
                link_knotscape = ''
                for j in link_list:
                    link_knotscape += chr(abs(j) + offset[copysign(1,j)])
                link = 'braid' + link_knotscape
            except TypeError:
                try:
                    link_list = [list(j) for j in link]
                    # interpreted as pd-Code
                    link = 'pd' + str(link_list)
                except TypeError:
                    pass

        arg_list.append(link)

        if not command:
            command = 'calc0'
        cmd_list = command.split()
        arg_list += cmd_list

        # redirect print outputs to a list to
        # return it optionally to the user
        print_output = []
        def no_print(s):
            print_output.append(s)

        if verbose:
            if print_messages:
                raise ValueError('All messages are printed in verbose mode')
            with self._stack:
                res = run_commandline(arg_list, print, progress)
            return res

        # capture print output from the shared library
        # and use the to set the description resp.
        # return it optionally to the user
        import py
        capture = py.io.StdCaptureFD(err=False, in_=False)

        try:
            with self._stack:
                res = run_commandline(arg_list, no_print, progress)
        except:
            out, err = capture.reset()
            raise

        out, err = capture.reset()
        if not self._description and out:
            self._description = out.strip('\n')

        if print_messages:
            if  out:
                print_output += [out][:-1]
            return res, print_output
        return res

    def khoca_version(self):
        r"""
        Return the version if ``Khoca``.

        EXAMPLES::

            >>> from khoca import InteractiveCalculator
            >>> KH = InteractiveCalculator()
            >>> KH.khoca_version() >= '1.1'
            True
        """
        from khoca._version import version
        return version
