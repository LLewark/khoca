#!/usr/bin/python3
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#
# This file contains test functions to be perfomed on cibuildwheel
# It performs the example of :class:`InteractiveCalculator` and
# does some cross-checks against results calculated with KnotJob
# and stored in the  databases at the web-pages
# `KnotInfo <https://knotinfo.math.indiana.edu/>`__ and
# `LinkInfo <https://linkinfo.sitehost.iu.edu/>`__


from database_knotinfo import link_list
from khoca import InteractiveCalculator
from sympy import poly
from sympy.parsing.sympy_parser import (parse_expr, standard_transformations, implicit_multiplication)

transformations = standard_transformations + (implicit_multiplication,)
varnames = ['q', 't', 'T', '1/q', '1/t', '1/T']

def khoca2dict(l, characteristic=0, reduced=False):
    """
    According to  src/python_interface/pui.pyx  printFormattedHomology:
    l = eval(s)
    l = [[i[2], i[0], i[1], i[3]] for i in l]
    # i[0]: torsion, i[1]: t-degree, i[1]: q-degree, i[3]: coefficient

    here:
    # i[2]: torsion, i[0]: t-degree, i[1]: q-degree, i[3]: coefficient
    """
    if reduced:
        data = l[0]
    else:
        data = l[1]

    data_as_dict = {}
    for i in data:
        # i[0]: degree, i[1]: height, i[2]: torsion i[3]: rank
        d, h, t, r = i
        d = int(t/2) - d # for compatibility with KnotInfo
        if (h, d, t) in data_as_dict:
            data_as_dict[(h, d, t)] += r
        else:
            data_as_dict[(h, d, t)] = r
    return {k: v for k, v in data_as_dict.items() if v != 0}

def convert_key(k, gens):
    r"""
    Convert a key ``k`` from the SymPy dictionary of the Khovanov polynomial from KnotInfo
    to a key of the dictionary returned by ``khoca2dict``.
    """
    res = [0, 0, 0]
    for i in range(len(k)):
        j = varnames.index(gens[i])
        if j < 3:
            res[j] += k[i]
        else:
            res[j - 3] -= k[i]
    return tuple(res)

def knotInfo2dict(string):
    r"""
    Return a dictionary representing the Khovanov polynomial from KnotInfo given in ``string``.
    """
    new_string = parse_expr(string.replace('^', '**'), transformations=transformations)
    p = poly(new_string)
    gens = [str(g) for g in p.gens]
    p_dict = p.as_dict()
    p_dict_red = {convert_key(k, gens): v for k, v in p_dict.items()}
    return p_dict_red

def cmp_knotinfo(max_cross_number=9, multi_component_links=False):
    r"""
    Perform a cross-check between the results obtained from ``Khoca`` and the
    results from ``KnotJob`` listed in the KnotInfo database.
    """
    if multi_component_links:
       print('Start KnotInfo cross-check of multicomponent links up to %s crossings' % max_cross_number)
       phrase='links'
       # for multicomponent link KnotInfo just provides rational homology
       KH = InteractiveCalculator(coefficient_ring=1)
    else:
       print('Start KnotInfo cross-check of knots up to %s crossings' % max_cross_number)
       phrase='knots'
       KH = InteractiveCalculator()
    kl = link_list(proper_links=multi_component_links)
    klred = [k for k in kl if kl.index(k) > 1 and int(k['crossing_number']) <= max_cross_number]
    count = 0
    for k in klred:
        name = k['name']
        braid = eval(k['braid_notation'].replace('{', '(').replace('}', ')'))
        if multi_component_links:
            braid = braid[1]
            kp_raw = k['khovanov_polynomial']
        else:
            kp_raw = k['khovanov_unreduced_integral_polynomial']
        kh_raw = KH(braid)
        kp_dict = knotInfo2dict(kp_raw)
        kh_dict = khoca2dict(kh_raw)
        count += 1
        if kp_dict != kh_dict:
            print('different results for %s: kh %s <> ki %s' % (name, kh_dict, kp_dict))
            assert(kp_dict == kh_dict)
        if count % 20 == 0:
            print('comparison of %s %s passed' % (count, phrase))

    print('comparison of %s %s completely passed' % (count, phrase))

def test_khoca():
    r"""
    Perform the tests.
    """
    # Test according to the examples of __init__.py
    KH = InteractiveCalculator()
    print(KH)
    print(KH('braidaBaB'))
    res_trefoil, messages = KH('braidaaa', print_messages=True)
    print(messages)
    assert(KH((1, 1, 1)) == res_trefoil)
    res_figure_eight = KH('braidaBaB')
    assert(KH([[4,2,5,1],[8,6,1,5],[6,3,7,4],[2,7,3,8]]) == res_figure_eight)
    assert(KH((1, -2, 1, -2)) == res_figure_eight)
    print(InteractiveCalculator(1, (0, 1), 0))
    KH = InteractiveCalculator(equivariant=3)
    print(KH)
    KH = InteractiveCalculator(1, (1, 0), 0)
    print(KH)
    print(KH('braidaaa', 'dual0 sum0+1 braidaBaB sum2+3 calc4', print_messages=True))

    # Cross-checks against the KnotInfo database
    cmp_knotinfo()
    cmp_knotinfo(multi_component_links=True)


test_khoca()
