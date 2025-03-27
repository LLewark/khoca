/*
 *
 *    src/planar_algebra/explicitTemplates.cpp --- Part of khoca, a knot homology calculator
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
#define INSTANTIATE_PRINTSIZE(COBO) \
template void Complex<COBO>::printSize(std::vector<word64> &) const; \
template void LCCobos<COBO>::printSize(std::vector<word64> &) const;
#else
#define INSTANTIATE_PRINTSIZE(COBO)
#endif
#define INSTANTIATE_ALL(COBO) \
INSTANTIATE_PRINTSIZE(COBO) \
template Complex<COBO>::Complex(const Complex<COBO> &, const Complex<COBO> &);\
template void Complex<COBO>::glue(const boundary_t[2]); \
template int Complex<COBO>::simplifyOnce(int, int, int); \
template int Complex<COBO>::simplifyOnceAtTDegree(int, int, const std::string*, int, int); \
template Complex<COBO>::Complex(std::ifstream &); \
template void Complex<COBO>::writeToBin(std::ofstream &) const; \
template void Complex<COBO>::deleteNonIsos(); \
template void Complex<COBO>::reducify(int); \
template int Complex<COBO>::isSimplified(int qDiff, int tDegree); \
template int Complex<COBO>::isSimplified(int qDiff); \
template void LCCobos<COBO>::compose(const LCCobos<COBO> &, const tangle_t&, const tangle_t&, const tangle_t&); \
template LCCobos<COBO>::LCCobos(std::ifstream &, bool); \
template void LCCobos<COBO>::writeToBin(std::ofstream &) const; \
template void LCCobos<COBO>::add(LCCobos<COBO>); \
template void LCCobos<COBO>::setToNegInv(LCCobos<COBO>); \
template bool LCCobos<COBO>::isInvertible(const tangle_t &, const tangle_t &) const; \
template void LCCobos<COBO>::add(COBO x); \
template std::ostream& LCCobos<COBO>::detailedOutput(std::ostream &) const;

/** Macro parameters cannot contain commas, so we have to use these typedefs.
 * The following workaround requires C++11.
 */
template <unsigned bitSize> using KC0 = KrasnerCobo<MInteger, bitSize>;
template <unsigned bitSize> using KC1 = KrasnerCobo<MRational, bitSize>;
template <unsigned bitSize> using KC8 = KrasnerCobo<FF<uint8_t>, bitSize>;
template <unsigned bitSize> using KC16 = KrasnerCobo<FF<uint16_t>, bitSize>;
template <unsigned bitSize> using PKC0  = KrasnerCobo<Polynomial<MInteger>, bitSize>;
template <unsigned bitSize> using PKC1  = KrasnerCobo<Polynomial<MRational>, bitSize>;
template <unsigned bitSize> using PKC8  = KrasnerCobo<Polynomial<FF<uint8_t> >, bitSize>;
template <unsigned bitSize> using PKC16 = KrasnerCobo<Polynomial<FF<uint16_t> >, bitSize>;


#define INSTANTIATE_META(BITSIZE) \
INSTANTIATE_ALL(KC0<BITSIZE>) \
INSTANTIATE_ALL(KC1<BITSIZE>) \
INSTANTIATE_ALL(KC8<BITSIZE>) \
INSTANTIATE_ALL(KC16<BITSIZE>) \
INSTANTIATE_ALL(PKC0<BITSIZE>) \
INSTANTIATE_ALL(PKC1<BITSIZE>) \
INSTANTIATE_ALL(PKC8<BITSIZE>) \
INSTANTIATE_ALL(PKC16<BITSIZE>)

MULTICALLER

#undef INSTANTIATE_ALL
#undef INSTANTIATE_META
#undef INSTANTIATE_PRINTSIZE

template class FF<uint8_t>;
template class FF<uint16_t>;
template class Polynomial<MInteger>;
template class Polynomial<MRational>;
template class Polynomial<FF<uint8_t> >;
template class Polynomial<FF<uint16_t> >;

template class VecTangles<KrasnerTangle>;
template void VecTangles<KrasnerTangle>::getQs(qShift_t &qMax, qShift_t &qMin, std::vector<int> &idxTranslator, std::vector<int> &qDims) const;
