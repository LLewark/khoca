/*
 *
 *    src/planar_algebra/sparsematExplicitTemplates.cpp --- Part of khoca, a knot homology calculator
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
template void SparseMat<LCCobos<COBO> >::printSize(std::vector<word64> &) const;
#else
#define INSTANTIATE_PRINTSIZE(COBO)
#endif

#ifndef notests
#define INSTANTIATE_ISSANE(COBO) \
template bool SparseMat<LCCobos<COBO> >::isSane() const;
#else
#define INSTANTIATE_ISSANE(COBO)
#endif

#define INSTANTIATE_ALL(COBO) \
INSTANTIATE_PRINTSIZE(COBO) \
INSTANTIATE_ISSANE(COBO) \
template bool SparseMat<LCCobos<COBO> >::hasInvertibles() const; \
template SparseMat<LCCobos<COBO> >::idx_t SparseMat<LCCobos<COBO> >::numberOfInvertibles() const; \
template void SparseMat<LCCobos<COBO> >::noLongerInvertible(idx_t idx); \
template void SparseMat<LCCobos<COBO> >::isNowInvertible(idx_t idx); \
template SparseMat<LCCobos<COBO> >::SparseMat(idx_t, idx_t); \
template SparseMat<LCCobos<COBO> >::SparseMat(std::ifstream&, bool); \
template void SparseMat<LCCobos<COBO> >::writeToBin(std::ofstream&) const; \
template void SparseMat<LCCobos<COBO> >::eraseRow(idx_t); \
template void SparseMat<LCCobos<COBO> >::eraseCol(idx_t); \
template void SparseMat<LCCobos<COBO> >::addCols(idx_t n); \
template void SparseMat<LCCobos<COBO> >::addRows(idx_t n); \
template void SparseMat<LCCobos<COBO> >::setLastEntry(idx_t, LCCobos<COBO>, bool isInv); \
template void SparseMat<LCCobos<COBO> >::setLastEntry(idx_t, idx_t, LCCobos<COBO>&&, bool); \
template LCCobos<COBO>* SparseMat<LCCobos<COBO> >::getEntry(idx_t, idx_t, idx_t*); \
template const LCCobos<COBO>* SparseMat<LCCobos<COBO> >::getEntry(idx_t, idx_t, idx_t*) const; \
template int SparseMat<LCCobos<COBO> >::setEntry(idx_t, idx_t, LCCobos<COBO>, bool, idx_t*); \
template void SparseMat<LCCobos<COBO> >::copyRow(idx_t, idx_t); \
template void SparseMat<LCCobos<COBO> >::copyCol(idx_t, idx_t); \
template void SparseMat<LCCobos<COBO> >::extractRow(idx_t, idx_t, valCont_t &, uintCont_t &) const; \
template void SparseMat<LCCobos<COBO> >::extractCol(idx_t, idx_t, valCont_t &, uintCont_t &) const; \
template SparseMat<LCCobos<COBO> >::idx_t SparseMat<LCCobos<COBO> >::getColCount() const; \
template SparseMat<LCCobos<COBO> >::idx_t SparseMat<LCCobos<COBO> >::getRowCount() const; \
template SparseMat<LCCobos<COBO> >::idx_t SparseMat<LCCobos<COBO> >::getEntryCount() const; \
template void SparseMat<LCCobos<COBO> >::reserve(int, int, int); \
template SparseMat<LCCobos<COBO> > SparseMat<LCCobos<COBO> >::setToDual(SparseMat<LCCobos<COBO> > const&); \
template bool SparseMat<LCCobos<COBO> >::stepToNextInv(idx_t&, idx_t&, idx_t&, LCCobos<COBO> *&); \
template void SparseMat<LCCobos<COBO> >::swap(SparseMat<LCCobos<COBO> >&); \
template void SparseMat<LCCobos<COBO> >::setRowNumber(idx_t); \
template class GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>; \
template void GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::setToMatBegin(SparseMat<LCCobos<COBO> > const&); \
template bool GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::isOn() const; \
template void GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::stepAlongMat(); \
template const LCCobos<COBO>* SMconstIterator<LCCobos<COBO> >::getVal() const; \
template void SMIterator<LCCobos<COBO> >::stepAlongRow(bool); \
template void SMIterator<LCCobos<COBO> >::stepAlongCol(bool); \
template void SMIterator<LCCobos<COBO> >::stepAlongMat(bool); \
template LCCobos<COBO>* SMIterator<LCCobos<COBO> >::getVal(); \
template void GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::setToRowBegin(SparseMat<LCCobos<COBO> >&, GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::idx_t); \
template void GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::stepAlongRow(); \
template bool GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::isOn() const; \
template GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::idx_t GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::getCol() const; \
template void GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::stepAlongMat(); \
template void GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::setToMatBegin(SparseMat<LCCobos<COBO> >&); \
template GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::idx_t GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::getRow() const; \
template GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::idx_t GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::getIdx() const; \
template void GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::setToColBegin(SparseMat<LCCobos<COBO> >&, GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::idx_t); \
template void GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::stepAlongCol(); \
template GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::idx_t GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::getRow() const; \
template GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::idx_t GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::getIdx() const; \
template void GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::stepAlongCol(); \
template void GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::stepAlongRow(); \
template void GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::setToColBegin(const SparseMat<LCCobos<COBO> >&, GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::idx_t); \
template void GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::setToRowBegin(const SparseMat<LCCobos<COBO> >&, GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::idx_t); \
template GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::idx_t GeneralIterator<SparseMat<LCCobos<COBO> > const, LCCobos<COBO> const>::getCol() const; \
template void GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::correctIdx(GeneralIterator<SparseMat<LCCobos<COBO> >, LCCobos<COBO> >::idx_t); \
template std::ostream& SparseMat<LCCobos<COBO> >::detailedOutput(std::ostream &) const; \



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
#undef INSTANTIATE_ISSANE
