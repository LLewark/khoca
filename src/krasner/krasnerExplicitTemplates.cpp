/*
 *
 *    src/krasner/krasnerExplicitTemplates.cpp --- Part of khoca, a knot homology calculator
 *
 * Copyright (C) 2018 Lukas Lewark <lukas@lewark.de>
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef getsize
#define INSTANTIATE_PRINTSIZE(COEFF, BITSIZE) \
template void KrasnerCobo<COEFF, BITSIZE>::printSize( \
        std::vector<word64> &) const;
#else
#define INSTANTIATE_PRINTSIZE(COEFF, BITSIZE)
#endif

#ifndef notests
#define SANITY(COEFF, BITSIZE)
#else
#define SANITY(COEFF, BITSIZE)
#endif

#define INSTANTIATE_ALL(COEFF, BITSIZE) \
template          KrasnerCobo<COEFF, BITSIZE>::KrasnerCobo(); \
template          KrasnerCobo<COEFF, BITSIZE>::KrasnerCobo(const KrasnerTangle &); \
template          KrasnerCobo<COEFF, BITSIZE>::KrasnerCobo(std::ifstream &f, bool); \
template qShift_t KrasnerCobo<COEFF, BITSIZE>::degree(boundary_t boundarySize) const; \
template qShift_t KrasnerCobo<COEFF, BITSIZE>::relativeDegree() const; \
template void     KrasnerCobo<COEFF, BITSIZE>::compose(const KrasnerCobo<COEFF, BITSIZE> &, std::vector<this_t> &, const tangle_t &, const tangle_t &, const tangle_t &) const; \
template void     KrasnerCobo<COEFF, BITSIZE>::setToUnion(const KrasnerTangle&, const KrasnerTangle&, const KrasnerTangle&, const KrasnerTangle&, const this_t &o1, const this_t &o2); \
template void     KrasnerCobo<COEFF, BITSIZE>::glue(const KrasnerTangle&, const KrasnerTangle&, const boundary_t[2], boundary_t, LCCobos<this_t> &); \
template void     KrasnerCobo<COEFF, BITSIZE>::print() const; \
template void     KrasnerCobo<COEFF, BITSIZE>::setToRandom(); \
template void     KrasnerCobo<COEFF, BITSIZE>::modifyDeloopCopy(int, bool, std::vector<this_t> &, const KrasnerTangle&, const KrasnerTangle&); \
template bool     KrasnerCobo<COEFF, BITSIZE>::operator<(const this_t &) const; \
template bool     KrasnerCobo<COEFF, BITSIZE>::operator==(const this_t &) const; \
template bool     KrasnerCobo<COEFF, BITSIZE>::operator!=(const this_t &) const; \
template int      KrasnerCobo<COEFF, BITSIZE>::reducify(); \
template bool     KrasnerCobo<COEFF, BITSIZE>::isInvertible(const tangle_t &, const tangle_t &) const; \
template bool     KrasnerCobo<COEFF, BITSIZE>::isEmpty() const; \
template void     KrasnerCobo<COEFF, BITSIZE>::writeToBin(std::ofstream &f) const; \
SANITY(COEFF, BITSIZE) \

#define INSTANTIATE_META(BITSIZE) \
INSTANTIATE_ALL(MInteger, BITSIZE) \
INSTANTIATE_ALL(MRational, BITSIZE) \
INSTANTIATE_ALL(FF<uint8_t>, BITSIZE) \
INSTANTIATE_ALL(FF<uint16_t>, BITSIZE) \
INSTANTIATE_ALL(Polynomial<MInteger>, BITSIZE) \
INSTANTIATE_ALL(Polynomial<MRational>, BITSIZE) \
INSTANTIATE_ALL(Polynomial<FF<uint8_t> >, BITSIZE) \
INSTANTIATE_ALL(Polynomial<FF<uint16_t> >, BITSIZE)

MULTICALLER

#undef INSTANTIATE_ALL
#undef INSTANTIATE_META
#undef INSTANTIATE_PRINTSIZE
#undef SANITY
