/*
 *
 *    src/planar_algebra/planar_algebra.cpp --- Part of khoca, a knot homology calculator
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

#include <mutex>
#include <condition_variable>
#include <thread>

#include <vector>
#include <deque>
#include <bitset>
#include <stdlib.h>
#include <inttypes.h>
#include <sstream>
#include <algorithm>
// Next line is necessary due to a bug of mpir
#include <stddef.h>
#include <gmp.h>
#include <fstream>
#include <iostream>
#include <iomanip>

#include "../shared.h"

#include <assert.h>

#include "coefficient_rings.h"
#include "sparsemat.h"
#include "planar_algebra.h"
#include "../krasner/krasner.h"

#include "explicitTemplates.cpp"

/** Sizes.
 * 0. Complex
 * 1. Vector
 * 2. Tangles
 * 3. Matrix
 * 4. LCCobos
 * 5. Cobos
 * 6. Coefficients
 * 7. New Cobos
 */
#ifndef getsize
template <class cobordism_tpl>
void Complex<cobordism_tpl>::printSize(std::vector<word64> &s) const {
    s.at(0) += sizeof(*this);
    s.at(1) += sizeof(vecTangles_t) * vecTangles.capacity();
    s.at(3) += sizeof(matLCCobos_t) * matLCCobos.capacity();
    for (auto i = vecTangles.cbegin(); i != vecTangles.cend(); ++i)
        i->printSize(s);
    for (auto i = matLCCobos.cbegin(); i != matLCCobos.cend(); ++i)
        i->printSize(s);
}

template <class tangle_tpl>
void VecTangles<tangle_tpl>::printSize(std::vector<word64> &s) const {
    s.at(1) += sizeof(word64) * deloopStack.capacity();
    for (typename tangleCont_t::const_iterator i = tangles.begin();
            i != tangles.end(); ++i)
        i->printSize(s);
}

template <class cobordism_tpl>
void MatLCCobos<cobordism_tpl>::printSize(std::vector<word64> &s) const {
    morphisms.printSize(s);
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::printSize(std::vector<word64> &s) const {
    for (auto i = cobordisms.cbegin(); i != cobordisms.cend(); ++i)
        i->printSize(s);
}
#endif

template <class cobordism_tpl>
Complex<cobordism_tpl>* Complex<cobordism_tpl>::setToDualConcrete(
        const Complex<cobordism_tpl> &other) {
    Complex<cobordism_tpl> *dual = new Complex<cobordism_tpl> (
            other.boundary.size(),
            - other.matLCCobos.size() - other.globalTShift);
    dual->vecTangles.reserve(other.vecTangles.size());
    dual->matLCCobos.reserve(other.matLCCobos.size());
    for (auto i = other.vecTangles.rbegin(); i != other.vecTangles.rend(); ++i)
        dual->vecTangles.push_back(std::move(vecTangles_t::setToDual(*i)));
    for (auto i = other.matLCCobos.rbegin(); i != other.matLCCobos.rend(); ++i)
        dual->matLCCobos.push_back(std::move(matLCCobos_t::setToDual(*i)));
    return dual;
}

template <class tangle_tpl>
VecTangles<tangle_tpl> VecTangles<tangle_tpl>::setToDual(
        const VecTangles<tangle_tpl> &other) {
    VecTangles<tangle_tpl> dual;
    dual.deloopStack = other.deloopStack;
    dual.tangles.reserve(other.tangles.size());
    for (auto i = other.tangles.begin(); i != other.tangles.end(); ++i)
        dual.tangles.push_back(std::move(Tangle<tangle_tpl>::setToDual(*i)));
    return dual;
}

template <class tangle_tpl>
tangle_tpl Tangle<tangle_tpl>::setToDual(const tangle_tpl &other) {
    tangle_tpl dual = other;
    dual.qShift *= -1;
    return dual;
}

template <class cobordism_tpl>
MatLCCobos<cobordism_tpl> MatLCCobos<cobordism_tpl>::setToDual(
        const MatLCCobos<cobordism_tpl> &other) {
    MatLCCobos<cobordism_tpl> dual;
    
    dual.morphisms = std::move(SparseMat<LCCobos<cobordism_tpl> >::setToDual(
                other.morphisms));
    return dual;
}

template <class cobordism_tpl>
LCCobos<cobordism_tpl>::LCCobos(std::ifstream &f, bool intCoefficients) {
    uint64_t size;
    readFromBinTpl(f, size);
    cobordisms.reserve(size);
    for (word64 i = W64LIT(0); i < (word64)size; ++i)
        cobordisms.emplace_back(f, intCoefficients);
}

template <class tangle_tpl>
VecTangles<tangle_tpl>::VecTangles(std::ifstream &f, boundary_t boundarySize) {
    uint64_t size, deloopStackSize;
    readFromBinTpl(f, size);
    tangles.reserve(size);
    for (word64 i = W64LIT(0); i < (word64)size; ++i)
        tangles.emplace_back(f, boundarySize);
    readFromBinTpl(f, deloopStackSize);
    for (word64 i = W64LIT(0); i < (word64)deloopStackSize; ++i) {
        uint64_t v;
        readFromBinTpl(f, v);
        deloopStack.emplace_back(std::move(v));
    }
}

template <class cobordism_tpl>
Complex<cobordism_tpl>::Complex(std::ifstream &f) {
    boundary_t boundarySize;
    uint32_t vecTanglesSize;
    uint16_t coefficientRing;
    readFromBinTpl(f, coefficientRing);
    if ((coefficientRing != 0) && (coefficientRing !=
                cobordism_tpl::coeff_t::coefficientTypeToUint())) { 
        std::cerr << "File to be loaded is over ground ring " << coefficientRing
            << ", but complex is over ground ring " <<
            cobordism_tpl::coeff_t::coefficientTypeToUint() << ".";
        throw;
    }
    readFromBinTpl(f, globalTShift);
    readFromBinTpl(f, boundarySize);
    readFromBinTpl(f, vecTanglesSize);
    boundary.setToSize(boundarySize);
    vecTangles.reserve(vecTanglesSize);
    word64 matLCCobosSize = vecTanglesSize ? W64LIT(vecTanglesSize - 1) : W64LIT(0);
    matLCCobos.reserve(matLCCobosSize);
    for (word64 i = W64LIT(0); i < (word64)vecTanglesSize; ++i)
        vecTangles.emplace_back(f, boundarySize);
    for (word64 i = W64LIT(0); i < matLCCobosSize; ++i)
        matLCCobos.emplace_back(f, coefficientRing == 0);
    assert(isSane());
}


template <class cobordism_tpl>
std::ostream& LCCobos<cobordism_tpl>::detailedOutput(std::ostream &os) const {
    for (auto i = cobordisms.begin(); i != cobordisms.end(); ++i) {
        i->detailedOutput(os);
        if ((i + 1) != cobordisms.end())
            os << " + ";
    }
    return os;
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::writeToBin(std::ofstream &f) const {
    writeToBinTpl(f, (uint64_t)cobordisms.size());
    for (typename std::vector<cobordism_tpl>::const_iterator i =
            cobordisms.begin(); i != cobordisms.end(); ++i)
        i->writeToBin(f);
}

template <class cobordism_tpl>
void MatLCCobos<cobordism_tpl>::writeToBin(std::ofstream &f) const {
    morphisms.writeToBin(f);
}

template <class tangle_tpl>
void VecTangles<tangle_tpl>::writeToBin(std::ofstream &f) const {
    writeToBinTpl(f, (uint64_t)tangles.size());
    for (typename tangleCont_t::const_iterator i = tangles.begin();
            i != tangles.end(); ++i)
        i->writeToBin(f);
    writeToBinTpl(f, (uint64_t)deloopStack.size());
    for (word64 i = W64LIT(0); i < (word64)deloopStack.size(); ++i)
        writeToBinTpl(f, deloopStack.at(i));
}

template <class cobordism_tpl>
void Complex<cobordism_tpl>::writeToBin(std::ofstream &f) const {
    writeToBinTpl(f, cobordism_tpl::coeff_t::coefficientTypeToUint());
    writeToBinTpl(f, globalTShift);
    writeToBinTpl(f, (boundary_t)boundary.size());
    writeToBinTpl(f, (uint64_t)vecTangles.size());
    for (typename vecTanglesCont_t::const_iterator i = vecTangles.begin();
            i != vecTangles.end(); ++i)
        i->writeToBin(f);
    for (typename matLCCobosCont_t::const_iterator i = matLCCobos.begin();
            i != matLCCobos.end(); ++i)
        i->writeToBin(f);
}

#ifndef notests
template <class cobordism_tpl>
bool Complex<cobordism_tpl>::isSane() const {
    if ((vecTangles.size() != 0) &&
            (vecTangles.size() != matLCCobos.size() + 1)) {
	insane();
        return false;
    }
    for (int i = 0; i < (int)vecTangles.size(); ++i) {
        if (! vecTangles.at(i).isSane(&boundary))
            return false;
        if ((i) && (! matLCCobos.at(i - 1).isSane(&boundary,
                        &vecTangles.at(i - 1).getTangles(),
                        &vecTangles.at(i).getTangles())))
            return false;
        if (i && ((i + 1) != (int)vecTangles.size()) &&
                matLCCobos.at(i-1).multIsNonZero(matLCCobos.at(i), 
                    vecTangles.at(i - 1).getTangles(),
                    vecTangles.at(i).getTangles(),
                    vecTangles.at(i + 1).getTangles())) {
	    insane();
            return false;
	}
    }
    return true;
}

template <class cobordism_tpl>
bool MatLCCobos<cobordism_tpl>::multIsNonZero(const MatLCCobos<cobordism_tpl>
        &other, const tangleCont_t &lower, const tangleCont_t &middle,
        const tangleCont_t &upper) const {

    assert(lower.size() == morphisms.getColCount());
    assert(middle.size() == other.morphisms.getColCount());
    assert(middle.size() == morphisms.getRowCount());
    assert(upper.size() == other.morphisms.getRowCount());

    for (long i = 0; i < (long)lower.size(); ++i)
        for (long j = 0; j < (long)upper.size(); ++j) {
            constIterator_t l, r;
            l.setToColBegin(this->morphisms, i);
            r.setToRowBegin(other.morphisms, j);
            LCCobos_t result;
            while (l.isOn() && r.isOn()) {
                const long lIdx = l.getRow();
                const long rIdx = r.getCol();
                if (lIdx == rIdx) {
                    LCCobos_t toAdd = *(l.getVal());
                    toAdd.compose(*(r.getVal()), lower.at(i), middle.at(lIdx),
                            upper.at(j));
                    result.add(std::move(toAdd));
                }
                if (lIdx <= rIdx)
                    l.stepAlongCol();
                else
                    r.stepAlongRow();
            }
            if (! result.isZero()) {
                std::cout << "Mult of matrices does not give zero.\n";
                return true;
            }
        }
    return false;
}

template <class tangle_tpl>
bool VecTangles<tangle_tpl>::isSane(const Boundary *b) const {
    word64 countDeloopables = W64LIT(0);
    for (typename tangleCont_t::const_iterator i = tangles.begin();
            i != tangles.end(); ++i) {
        if (! i->isSane(b))
            return false;
        if (i->hasLoop())
            countDeloopables += 1;
    }
    for (std::vector<word64>::const_iterator i = deloopStack.begin();
            i != deloopStack.end(); ++i)
        if (! tangles.at(*i).hasLoop()) {
	    insane();
            return false;
	}
    return true;
}

template <class cobordism_tpl>
bool MatLCCobos<cobordism_tpl>::isSane(const Boundary *b,
        const tangleCont_t *domain, const tangleCont_t *codomain) const {
    if (! morphisms.isSane())
        return false;

    bool detail = b || domain || codomain;
    if (detail) {
        assert(b && domain && codomain);

        typename SparseMat<LCCobos<cobordism_tpl> >::uintCont_t inv2;
        for (constIterator_t i = getConstIterator(); i.isOn();
                i.stepAlongMat()) {
            if (i.getVal()->isInvertible(domain->at(i.getCol()),
                        codomain->at(i.getRow())))
                inv2.push_back(i.getIdx());
            if (i.getVal()->isZero()) {
		insane();
                return false;
	    }
        }
        if (inv2 != morphisms.getInvertibles()) {
	    insane();
            return false;
	}
    }
    for (constIterator_t i = getConstIterator(); i.isOn(); i.stepAlongMat()) {
        if (detail) {
            if (! i.getVal()->isSane(b, &domain->at(i.getCol()),
                        &codomain->at(i.getRow())))
                return false;
        } else if (! i.getVal()->isSane())
            return false;
    }
    return true;
}

template <class cobordism_tpl>
bool LCCobos<cobordism_tpl>::isSane(const Boundary *b, const tangle_t *domain,
        const tangle_t *codomain) const {
    for (typename std::vector<cobordism_tpl>::const_iterator i =
            cobordisms.begin(); i != cobordisms.end(); ++i) {
        if (! i->isSane(b, domain, codomain))
            return false;
        if (! i->coefficient.isNonZero()) {
	    insane();
            return false;
	}
        if ((i != cobordisms.begin()) && (! (*(i - 1) < *i))) {
	    insane();
            return false;
	}
    }
    return true;
}
#endif

template <class cobordism_tpl>
void MatLCCobos<cobordism_tpl>::deloop(word64 idx, int copies,
        const tangleCont_t &lowerTangles, const tangleCont_t &upperTangles,
        bool left) {

    assert(isSane());
    if (left)
        morphisms.copyRow(idx, copies);
    else
        morphisms.copyCol(idx, copies);
    SMIterator<LCCobos<cobordism_tpl> > i;
    for (int j = 0; j <= copies; ++j) {
        const int k = j ? (left ? (morphisms.getRowCount() - copies + j - 1) :
                (morphisms.getColCount() - copies + j - 1)) : idx;
        if (left)
            i.setToRowBegin(morphisms, k);
        else
            i.setToColBegin(morphisms, k);
        while (i.isOn()) {
            const tangle_t &lowerTangle = lowerTangles.at(
                    left ? i.getCol() : k);
            const tangle_t &upperTangle = upperTangles.at(
                    left ? k : i.getRow());
            i.getVal()->modifyDeloopCopy(j, left, lowerTangle, upperTangle);
            const bool becameZero = i.getVal()->isZero();
            const bool becameInvertible = i.getVal()->isInvertible(
                    lowerTangle, upperTangle);
            if (! becameInvertible)
                morphisms.noLongerInvertible(i.getIdx());
            else 
                morphisms.isNowInvertible(i.getIdx());

            if (left)
                i.stepAlongRow(becameZero);
            else
                i.stepAlongCol(becameZero);
        }
    }
    assert(isSane());
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::modifyDeloopCopy(int kind, bool left,
        const typename cobordism_tpl::tangle_t &lower,
        const typename cobordism_tpl::tangle_t &upper) {
    std::vector<cobordism_tpl> oldCobos;
    oldCobos.swap(cobordisms);
    for (int i = 0; i < (int)oldCobos.size(); ++i)
        oldCobos.at(i).modifyDeloopCopy(kind, left, cobordisms, lower, upper);
    sortAndFactor();
}

template <class cobordism_tpl>
bool LCCobos<cobordism_tpl>::isInvertible(const tangle_t &lowerTangle,
        const tangle_t &upperTangle) const {
    return ((cobordisms.size() == 1) &&
            (cobordisms.front().coefficient.isInvertible()) &&
            (cobordisms.front().isInvertible(lowerTangle, upperTangle)));
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::add(cobordism_tpl x) {
    assert(isSane());
    const typename cobosCont_t::iterator i = std::lower_bound(
            cobordisms.begin(), cobordisms.end(), x);
    if ((i != cobordisms.end()) && (*i == x)) {
        i->coefficient += x.coefficient;
        if (! i->coefficient.isNonZero())
            cobordisms.erase(i);
    } else
        cobordisms.insert(i, std::move(x));
    assert(isSane());
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::factor() {
    int i = 0;
    while (i < (int)cobordisms.size()) {
        int j = i + 1;
        while ((j < (int)cobordisms.size()) &&
                (cobordisms.at(i) == cobordisms.at(j)))
            ++j;
        if (j > i + 1) {
            for (int k = i + 1; k < j; ++k)
                cobordisms.at(i).coefficient += cobordisms.at(k).coefficient;
            const bool notZero = cobordisms.at(i).coefficient.isNonZero();
            cobordisms.erase(cobordisms.begin() + i + (notZero ? 1 : 0),
                    cobordisms.begin() + j);
            if (notZero)
                ++i;
        } else
            ++i;
    }
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::sortAndFactor() {
    std::sort(cobordisms.begin(), cobordisms.end());
    factor();
}

template <class cobordism_tpl>
int Complex<cobordism_tpl>::isSimplified(int qDiff, int tDegree) {
    assert(isSane());
    assert((tDegree >= 0) && (tDegree < (int)vecTangles.size()));
    idx_t row, col;
    if (! vecTangles.at(tDegree).deloopingDone())
        return 0;
    if (tDegree < (int)matLCCobos.size()) {
        if (matLCCobos.at(tDegree).gaussianEliminationDone(row, col,
                    vecTangles.at(tDegree), vecTangles.at(tDegree + 1), qDiff))
            return 0;
        if (matLCCobos.at(tDegree).hasInvertibles())
            return 1;
        else return 2;
    } else
        return 2;
}

template <class cobordism_tpl>
int Complex<cobordism_tpl>::isSimplified(int qDiff) {
    int result = 2;
    for (int i = 0; i < (int)vecTangles.size(); ++i)
        switch (isSimplified(qDiff, i)) {
            case 0: return 0;
            case 1: result = 1;
        }
    return result;
}

template <class cobordism_tpl>
int Complex<cobordism_tpl>::simplifyOnceAtTDegree(
        int qDiff, int tDegree, const std::string* s, int numThreads, int progress) {
    assert(isSane());
    if (((tDegree < (int)matLCCobos.size()) &&
                tryToGauss(tDegree, qDiff, numThreads))
            || (tryToDeloop(tDegree))) {
        assert(isSane());
        if (progress)
            showProgressBar(s);
        return 0;
    }
    return ((tDegree < (int)matLCCobos.size()) &&
            matLCCobos.at(tDegree).hasInvertibles()) ? 1 : 2;
}

template <class cobordism_tpl>
bool Complex<cobordism_tpl>::tryToGauss(int i, int qDiff, int numThreads) {
    idx_t row, col;
    if (matLCCobos.at(i).gaussianElimination(row, col, vecTangles.at(i),
                vecTangles.at(i + 1), qDiff, numThreads)) {
        vecTangles.at(i).erase(col);
        vecTangles.at(i + 1).erase(row);
        if (i > 0)
            matLCCobos.at(i - 1).eraseRow(col);
        if ((i + 1) < (int)matLCCobos.size())
            matLCCobos.at(i + 1).eraseCol(row);
        assert(isSane());
        return true;
    }
    return false;
}

template <class cobordism_tpl>
bool Complex<cobordism_tpl>::tryToDeloop(int i) {
    word64 idx;
    const word64 oldSize = W64LIT(vecTangles.at(i).size());
    if ((idx = vecTangles.at(i).simplifyOnce()) != -1) {
        const int copies = vecTangles.at(i).size() - oldSize;
        if (i > 0)
            matLCCobos.at(i - 1).deloop(
                    idx, copies, vecTangles.at(i - 1).getTangles(),
                    vecTangles.at(i).getTangles(), true);
        if ((i + 1) < (int)vecTangles.size())
            matLCCobos.at(i).deloop(idx, copies, vecTangles.at(i).getTangles(),
                    vecTangles.at(i + 1).getTangles(), false);
        assert(isSane());
        return true;
    }
    return false;
}

template <class cobordism_tpl>
void Complex<cobordism_tpl>::showProgressBar(const std::string* s) {
    static std::mutex cprogressMutex;
    int simplificationsToBeDone = 0;
    for (typename matLCCobosCont_t::const_iterator i = matLCCobos.begin();
            i != matLCCobos.end(); ++i)
        simplificationsToBeDone += i->numberOfInvertibles();
    for (typename vecTanglesCont_t::const_iterator i = vecTangles.begin();
            i != vecTangles.end(); ++i)
        simplificationsToBeDone += i->deloopsToBeDone();
    simplificationsCount += 1;
    int percent = 100 * simplificationsCount /
        (simplificationsToBeDone + simplificationsCount);
    {
        std::unique_lock<std::mutex> lck (cprogressMutex);
        io::cprogress() << " " << std::setw(3) << percent << "% ["
            << std::string(percent/2, '=') << std::string(50 - percent/2,'.')
            << "] " << simplificationsCount << "/"
            << simplificationsToBeDone + simplificationsCount;
        if (s) io::cprogress() << *s;
        io::cprogress() << "\033[K\r" << std::flush;
    }
}

template <class cobordism_tpl>
int Complex<cobordism_tpl>::simplifyOnce(int qDiff, int numThreads, int progress) {
    assert(isSane());
    bool didSimplify = false;
    bool thereAreInvertibles = false;
    for (int i = 0; i < (int)matLCCobos.size(); ++i) {
        if (tryToGauss(i, qDiff, numThreads)) {
            didSimplify = true;
            break;
        }
        thereAreInvertibles = thereAreInvertibles ||
            matLCCobos.at(i).hasInvertibles();
    }
    if (! didSimplify)
        for (int i = 0; i < (int)vecTangles.size(); ++i)
            if (tryToDeloop(i)) {
                didSimplify = true;
                break;
            }
    assert(isSane());
    if (didSimplify && progress) {
        showProgressBar(nullptr);
        return 0;
    } else return thereAreInvertibles ? 1 : 2;
}

template <class cobordism_tpl>
typename MatLCCobos<cobordism_tpl>::LCCobos_t*
        MatLCCobos<cobordism_tpl>::gaussianEliminationDone(idx_t &row,
                idx_t &col, const vecTangles_t &domains,
                const vecTangles_t &codomains, int qDiff) {
    if (domains.size()) {
        Boundary b(domains.at(0).boundarySize());
        assert(isSane(&b, &(domains.getTangles()), &(codomains.getTangles())));
    }
    idx_t invIdx = -1;
    row = 0;
    LCCobos_t* iso;
    while (morphisms.stepToNextInv(invIdx, row, col, iso)) {
        const tangle_t &domain = domains.at(col);
        const tangle_t &codomain = codomains.at(row);
        assert(iso->isInvertible(domain, codomain));
        if (codomain.getQShift() != domain.getQShift() + qDiff)
            continue;
        return iso;
    }
    return nullptr;
}

template <class cobordism_tpl>
void MatLCCobos<cobordism_tpl>::gaussThread(const std::vector<LCCobos_t>
        &sameRowLCs, const std::vector<LCCobos_t> &sameColLCs,
        const std::vector<idx_t> sameRowIdxs, const std::vector<idx_t>
        sameColIdxs, bool compoIsTrivial, const LCCobos_t &isoInv, coeff_t y,
        const vecTangles_t &domains, const vecTangles_t &codomains,
        const tangle_t &domain, const tangle_t &codomain, int threadId,
        int numThreads, std::vector<LCCobos_t*> &result, std::mutex &resultMtx,
        std::condition_variable &compFinishedCv) {
    if ((int)sameRowLCs.size() <= threadId)
        return;
    for (auto i = sameRowLCs.begin() + threadId; i != sameRowLCs.end();) {
        idx_t c = sameRowIdxs.at(i - sameRowLCs.begin());
        for (auto j = sameColLCs.begin(); j != sameColLCs.end(); ++j) {
            idx_t r = sameColIdxs.at(j - sameColLCs.begin());
            LCCobos_t *toAdd = new LCCobos_t(*i);
            if (compoIsTrivial) {
                toAdd->switchSign();
                *toAdd *= y;
            } else
                toAdd->compose(isoInv, domains.at(c), codomain, domain);
            toAdd->compose(*j, domains.at(c), domain, codomains.at(r));
            {
                std::unique_lock<std::mutex> lck(resultMtx);
                result.at((j - sameColLCs.begin()) * sameRowLCs.size() +
                        (i - sameRowLCs.begin())) = toAdd;
            }
            compFinishedCv.notify_one();
        }
        if (sameRowLCs.end() - i > numThreads)
            i += numThreads;
        else
            i = sameRowLCs.end();
    }
}

/* In Pseudocode:
 * 1) Check that there really is something invertible where indicated, otherwise exit.
 * 2) Compute the inverse and save it.
 * 3) Extract the column and the row vector.
 * 4) Compute all products in the background.
 * 5) Recontrusct the matrix by iterating through the old matrix entries.
 * a) A double-iteration: One iterator goes through all old matrix entries; he always
 * knows at which index, row and col he is. The second iterator goes through the entries
 * in the row-vector, and for each, through the entries in the col vector. He also always
 * knows in which row and which col he is (but does not have an index). Double-iteration
 * works like this: at each step, work on the iterator that comes first, or on both if
 * they are at the same spot. "Working" means putting the value into newVal, the respective
 * column into newCol, updating newRowPtr if necessary, and writing newInvertibles as well;
 * and then stepping the iterator.
 */
template <class cobordism_tpl>
bool MatLCCobos<cobordism_tpl>::gaussianElimination(idx_t &row, idx_t &col,
        const vecTangles_t &domains, const vecTangles_t &codomains, int qDiff,
        int numThreads) {

    assert(numThreads >= 0);

    LCCobos_t* iso = gaussianEliminationDone(row, col, domains, codomains, qDiff);
    if (! iso)
        return false;
    const tangle_t &domain = domains.at(col);
    const tangle_t &codomain = codomains.at(row);
    LCCobos_t isoInv;
    coeff_t y;
    const bool compoIsTrivial = iso->compositionIsTrivial(domain, codomain);
    if (compoIsTrivial) {
        iso->getCoeff(y);
        y.inv();
    } else 
        isoInv.setToNegInv(*iso);

    std::vector<idx_t> sameRowIdxs, sameColIdxs;
    std::vector<LCCobos_t*> products;
    std::vector<std::thread> threads(numThreads);
    std::mutex resultMtx;
    std::condition_variable compFinishedCv;
        std::vector<LCCobos_t> sameRowLCs, sameColLCs;
        morphisms.extractRow(row, col, sameRowLCs, sameRowIdxs);
        morphisms.extractCol(col, row, sameColLCs, sameColIdxs);
        products.resize(sameRowIdxs.size() * sameColIdxs.size(), nullptr);
        for (auto i = threads.begin(); i != threads.end(); ++i)
            *i = std::thread(&MatLCCobos<cobordism_tpl>::gaussThread, this,
                    std::ref(sameRowLCs), std::ref(sameColLCs),
                    std::ref(sameRowIdxs), std::ref(sameColIdxs),
                    compoIsTrivial, std::ref(isoInv), y, std::ref(domains),
                    std::ref(codomains), std::ref(domain), std::ref(codomain),
                    i - threads.begin(), threads.size(), std::ref(products),
                    std::ref(resultMtx), std::ref(compFinishedCv));

    SparseMat<LCCobos<cobordism_tpl> > oldMorphisms(0, morphisms.getColCount() - 1);
    oldMorphisms.swap(morphisms);
    morphisms.reserve(oldMorphisms.getEntryCount() + (sameColIdxs.size() - 1) * (sameRowIdxs.size() - 1), oldMorphisms.getRowCount() - 1, oldMorphisms.numberOfInvertibles() + sameColIdxs.size() * sameRowIdxs.size());
    auto colIt = sameColIdxs.begin();
    auto rowIt = sameRowIdxs.begin();
    int prodIt = 0;
    bool addOn = products.size();

    SMIterator<LCCobos<cobordism_tpl> > i;
    i.setToMatBegin(oldMorphisms);
    while (i.isOn() && ((i.getCol() == col) || (i.getRow() == row)))
        i.stepAlongMat(false);
    bool oldOn = i.isOn();
    while (addOn || oldOn) {
        // Determine which is first
        bool doAdd = (addOn && ((!oldOn) || (*colIt < i.getRow()) || ((*colIt == i.getRow()) && (*rowIt <= i.getCol()))));
        bool doOld = (oldOn && ((!addOn) || (*colIt > i.getRow()) || ((*colIt == i.getRow()) && (*rowIt >= i.getCol()))));
        // Write entry of the new matrix accordingly
        if (doAdd) {
            for (;;) {
                std::unique_lock<std::mutex> lck(resultMtx);
                if (products.at(prodIt))
                    break;
                else
                    compFinishedCv.wait(lck);
            }
            if (doOld)
                products.at(prodIt)->add(std::move(*i.getVal()));
        }
        {
            idx_t toAddRow, toAddCol;
            LCCobos_t* toAddVal;
            bool isNonZero = true;
            if (doAdd) {
                toAddRow = *colIt;
                toAddCol = *rowIt;
                toAddVal = products.at(prodIt);
                isNonZero = ! products.at(prodIt)->isZero();
            }
            else if (doOld) {
                toAddRow = i.getRow();
                toAddCol = i.getCol();
                toAddVal = i.getVal();
            }
            else throw;
            if (isNonZero) {
                const bool isInv = toAddVal->isInvertible( domains.at(toAddCol), codomains.at(toAddRow));
                morphisms.setLastEntry(toAddCol - ((toAddCol > col) ? 1 : 0), toAddRow - ((toAddRow > row) ? 1 : 0), std::move(*toAddVal), isInv);
            }
        }
        // Advance iterators
        if (doAdd) {
            delete products.at(prodIt);
            ++prodIt;
            if (++rowIt == sameRowIdxs.end()) {
                rowIt = sameRowIdxs.begin();
                ++colIt;
                addOn = (colIt != sameColIdxs.end());
            }
        }
        if (doOld) {
            i.stepAlongMat(false);
            // Avoid the row and the column
            while (i.isOn() && ((i.getCol() == col) || (i.getRow() == row)))
                i.stepAlongMat(false);
            oldOn = i.isOn();
        }
    }
    morphisms.setRowNumber(oldMorphisms.getRowCount() - 1);
    for (auto i = threads.begin(); i != threads.end(); ++i)
        i->join();
    return true;
}

template <class tangle_tpl>
int VecTangles<tangle_tpl>::deloopsToBeDone() const {
    return deloopStack.size();
}

template <class tangle_tpl>
bool VecTangles<tangle_tpl>::deloopingDone() const {
    assert(isSane());
    return (deloopStack.empty());
}

template <class tangle_tpl>
word64 VecTangles<tangle_tpl>::simplifyOnce() {
    assert(isSane());
    if (deloopStack.empty())
        return -1;

    word64 criticalIndex = W64LIT(deloopStack.back());
    deloopStack.resize(deloopStack.size() - 1);
    tangle_tpl &criticalTangle = tangles.at(criticalIndex);

    tangleCont_t newCopies;
    criticalTangle.deloop(newCopies);

    if (criticalTangle.hasLoop())
        deloopStack.push_back(criticalIndex);
    for (word64 i = W64LIT(0); i < (int)newCopies.size(); ++i)
        if (newCopies.at(i).hasLoop())
            deloopStack.push_back(tangles.size() - newCopies.size() + i);

    tangles.reserve(tangles.size() + newCopies.size());
    std::move(newCopies.begin(), newCopies.end(), std::back_inserter(tangles));

    assert(isSane());
    return criticalIndex;
}

template <class tangle_tpl>
void VecTangles<tangle_tpl>::getQs(qShift_t &qMax, qShift_t &qMin,
        std::vector<int> &idxTranslator, std::vector<int> &qDims) const {
    idxTranslator.reserve(tangles.size());
    for (auto i = tangles.cbegin(); i != tangles.cend(); ++i) {
        int q = i->getQShift();
        if (i == tangles.begin()) {
            qMax = qMin = q;
            qDims.resize(1);
        }
        else {
            if (q > qMax) {
                qDims.resize((q - qMin) / 2 + 1);
                qMax = q;
            }
            if (q < qMin) {
                qDims.insert(qDims.begin(), (qMin - q) / 2, 0);
                qMin = q;
            }
        }
        idxTranslator.push_back(qDims.at((q - qMin) / 2)++);
    }
}

template <class cobordism_tpl>
void Complex<cobordism_tpl>::glue(const boundary_t gluePoints[2]) {
    assert(isSane());
    for (int i = 0; i < (int)matLCCobos.size(); ++i)
        matLCCobos.at(i).glue(gluePoints, vecTangles.at(i),
                vecTangles.at(i + 1), boundary.size());

    for (int i = 0; i < (int)vecTangles.size(); ++i)
        vecTangles.at(i).glue(gluePoints);

    boundary.glue();
#ifndef PRINTCOMPLEX
    std::cout << "Glued " << (int)gluePoints[0] << " to " << (int)gluePoints[1] << ". Output: ";
    detailedOutput(std::cout);
#endif
    assert(isSane());
}

template <class cobordism_tpl>
void Complex<cobordism_tpl>::reducify(int root) {
    for (typename matLCCobosCont_t::iterator i = matLCCobos.begin();
            i != matLCCobos.end(); ++i) {
        bool deleted = false;
        for (typename MatLCCobos<cobordism_tpl>::iterator_t j =
                i->getIterator(); j.isOn(); j.stepAlongMat(deleted)) {
            deleted = j.getVal()->reducify(root);
	    if (! deleted)
                i->isNowInvertible(j.getIdx());
	}
    }
}

template <class cobordism_tpl>
void Complex<cobordism_tpl>::deleteNonIsos() {
    for (typename matLCCobosCont_t::iterator i = matLCCobos.begin();
            i != matLCCobos.end(); ++i) {
        typename vecTanglesCont_t::iterator i2 =
            vecTangles.begin() + (i - matLCCobos.begin());
        for (typename MatLCCobos<cobordism_tpl>::iterator_t j =
                i->getIterator(); j.isOn(); j.stepAlongMat(
                    ! j.getVal()->isInvertible(i2->at(j.getCol()),
                        (i2 + 1)->at(j.getRow()))));
    }
}

template <class cobordism_tpl>
void MatLCCobos<cobordism_tpl>::deleteNonIsos() {
}

template <class tangle_tpl>
void VecTangles<tangle_tpl>::glue(const boundary_t gluePoints[2]) {
    assert(isSane());
    for (int i = 0; i < (int)tangles.size(); ++i) {
        tangles.at(i).glue(gluePoints);
        if (tangles.at(i).hasLoop())
            deloopStack.push_back(i);
    }
    assert(isSane());
}

template <class tangle_tpl>
void VecTangles<tangle_tpl>::appendTensorProduct(
        const VecTangles<tangle_tpl> &s1, const VecTangles<tangle_tpl> &s2) {
    assert(deloopStack.empty());
    assert(isSane());
    for (int i = 0; i < (int)s1.tangles.size(); ++i)
        for (int j = 0; j < (int)s2.tangles.size(); ++j) {
            tangle_tpl newTangle;
            newTangle.setToUnion(s1.tangles.at(i), s2.tangles.at(j));
            tangles.push_back(std::move(newTangle));
        }
    assert(isSane());
}

void Boundary::setToSum(const Boundary &b1, const Boundary &b2) {
    size_ = b1.size_ + b2.size_;
}

template <class cobordism_tpl>
Complex<cobordism_tpl>::Complex(const Complex<cobordism_tpl> &cc1,
        const Complex<cobordism_tpl> &cc2) :
    globalTShift(cc1.globalTShift + cc2.globalTShift) {

    const word64 size1 = W64LIT(cc1.vecTangles.size());
    const word64 size2 = W64LIT(cc2.vecTangles.size());

    // s is the sum, i is the first index.
    matLCCobos.resize(size1 + size2 - 2);
    for (word64 s = W64LIT(1); s + 1 < (int)(size1 + size2); ++s) {
        MatLCCobos<cobordism_tpl> &m = matLCCobos.at(s - 1);
        word64 colShift = W64LIT(0);

        const word64 forBegin = std::max(W64LIT(0ll), W64LIT(s - size2 + 1));
        const word64 forEnd = std::min(W64LIT(size1 - 1), W64LIT(s));
        for (word64 i = forBegin; i <= forEnd; ++i) {
            const word64 j = W64LIT(s - i);
            const VecTangles<typename cobordism_tpl::tangle_t> *domainVecs[2] =
            { i ? &(cc1.vecTangles.at(i - 1)) : nullptr,
                j ? &(cc2.vecTangles.at(j - 1)) : nullptr };
            const VecTangles<typename cobordism_tpl::tangle_t> *vecs[2] =
            { &(cc1.vecTangles.at(i)), &(cc2.vecTangles.at(j)) };
            const MatLCCobos<cobordism_tpl> *mats[2] =
            {(i > 0) ?  &cc1.matLCCobos.at(i - 1) : nullptr,
                (j > 0) ?  &cc2.matLCCobos.at(j - 1) : nullptr };
            if ((i == forBegin) && (mats[0]))
                m.addCols(mats[0]->getColCount() * vecs[1]->size());
            colShift += m.hossa(colShift, domainVecs, vecs, mats, j % 2);
        }
    }

    boundary.setToSum(cc1.boundary, cc2.boundary);
    vecTangles.resize(size1 + size2 - 1);
    for (word64 i = W64LIT(0); i < size1; ++i)
        for (word64 j = W64LIT(0); j < size2; ++j)
            vecTangles.at(i + j).appendTensorProduct(
                    cc1.vecTangles.at(i), cc2.vecTangles.at(j));

    assert(isSane());
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::compose(const LCCobos<cobordism_tpl> &other,
        const tangle_t &lower, const tangle_t &middle, const tangle_t &upper) {
    {
        Boundary b(lower.boundarySize());
        assert(isSane(&b, &lower, &middle));
        assert(other.isSane(&b, &middle, &upper));
    }
    cobosCont_t oldCobosLow;
    oldCobosLow.swap(cobordisms);
    cobordisms.reserve(oldCobosLow.size() * other.cobordisms.size());
    for (typename cobosCont_t::const_iterator i = oldCobosLow.begin();
            i != oldCobosLow.end(); ++i)
        for (typename cobosCont_t::const_iterator j = other.cobordisms.begin();
                j != other.cobordisms.end(); ++j)
            i->compose(*j, cobordisms, lower, middle, upper);
    sortAndFactor();
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::setToUnion(
        const typename cobordism_tpl::tangle_t& lowerLeft,
        const typename cobordism_tpl::tangle_t& lowerRight,
        const typename cobordism_tpl::tangle_t& upperLeft,
        const typename cobordism_tpl::tangle_t& upperRight,
        const LCCobos<cobordism_tpl> &xLeft,
        const LCCobos<cobordism_tpl> &xRight) {

    cobordisms.resize(xLeft.cobordisms.size() * xRight.cobordisms.size());
    for (int i = 0; i < (int)xLeft.cobordisms.size(); ++i)
        for (int j = 0; j < (int)xRight.cobordisms.size(); ++j)
            cobordisms.at(i * xRight.cobordisms.size() + j).setToUnion(
                    lowerLeft, lowerRight, upperLeft, upperRight,
                    xLeft.cobordisms.at(i), xRight.cobordisms.at(j));
    sortAndFactor();
    assert(! isZero());
    assert(isSane());
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::glue(const typename cobordism_tpl::tangle_t& lower,
        const typename cobordism_tpl::tangle_t& upper,
        const boundary_t gluePoint[2], boundary_t boundarySize) {

    assert(isSane());
    std::vector<cobordism_tpl> oldCobos;
    cobordisms.swap(oldCobos);
    for (typename std::vector<cobordism_tpl>::iterator i =
            oldCobos.begin(); i != oldCobos.end(); ++i)
        i->glue(lower, upper, gluePoint, boundarySize, *this);
    sortAndFactor();
    assert(isSane());
}


template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::add(LCCobos<cobordism_tpl> other) {
    assert(isSane());
    cobosCont_t oldCobos;
    oldCobos.swap(cobordisms);
    cobordisms.resize(oldCobos.size() + other.cobordisms.size());

    std::merge(std::make_move_iterator(oldCobos.begin()),
            std::make_move_iterator(oldCobos.end()), 
            std::make_move_iterator(other.cobordisms.begin()),
            std::make_move_iterator(other.cobordisms.end()),
            cobordisms.begin());
    factor();
    assert(isSane());
}

template <class cobordism_tpl>
void LCCobos<cobordism_tpl>::setToNegInv(LCCobos<cobordism_tpl> other) {
    assert(other.cobordisms.size() == 1);
    *this = std::move(other);
    cobordisms.front().setToNegInv();
    assert(isSane());
}

template <class cobordism_tpl>
bool LCCobos<cobordism_tpl>::reducify(int root) {
    for (typename std::vector<cobordism_tpl>::iterator i = cobordisms.begin();
	    i != cobordisms.end();) {
	const int numDots = i->reducify();
        if ((! root) && numDots)
            i = cobordisms.erase(i);
        else {
            for (int j = 0; j < numDots; ++j)
                i->coefficient *= root;
            ++i;
        }
    }
    sortAndFactor();
    return isZero();
}

template <class cobordism_tpl>
void MatLCCobos<cobordism_tpl>::glue(const boundary_t gluePoints[2],
        const VecTangles<typename cobordism_tpl::tangle_t> &lowerVec,
        const VecTangles<typename cobordism_tpl::tangle_t> &upperVec,
        boundary_t boundarySize) {
    assert(isSane());
    for (SMIterator<LCCobos<cobordism_tpl> > i = getIterator(); i.isOn();
            i.stepAlongMat(i.getVal()->isZero())) {
        const typename cobordism_tpl::tangle_t& lower = lowerVec.at(i.getCol());
        const typename cobordism_tpl::tangle_t& upper = upperVec.at(i.getRow());
        const bool oldInv = i.getVal()->isInvertible(lower, upper);
        i.getVal()->glue(lower, upper, gluePoints, boundarySize);
        const bool newInv = i.getVal()->isInvertible(lower, upper);
        if (oldInv && !newInv)
            morphisms.noLongerInvertible(i.getIdx());
        else if (newInv && !oldInv)
            morphisms.isNowInvertible(i.getIdx());
    }
    assert(isSane());
}

template <class cobordism_tpl>
int MatLCCobos<cobordism_tpl>::hossa(int colShift,
        const VecTangles<typename cobordism_tpl::tangle_t> *domainVecs[2],
        const VecTangles<typename cobordism_tpl::tangle_t> *vecs[2],
        const MatLCCobos<cobordism_tpl> *mats[2], bool signSwitch) {

    assert(isSane());
    assert(vecs[0] && vecs[1]);
    assert(mats[0] || mats[1]);
    assert(morphisms.isFullyGaussed());
    const int newCols = mats[1] ? mats[1]->getColCount() * vecs[0]->size() : 0;
    const int oldCols = mats[0] ? mats[0]->getColCount() * vecs[1]->size() : 0;
    const int newRows = vecs[0]->size() * vecs[1]->size();

    std::vector<LCCobos<cobordism_tpl> > idCobo[2];
    for (int i = 0; i < 2; ++i)
        if (mats[i]) {
            idCobo[1 - i].reserve(vecs[1 - i]->size());
            for (int j = 0; j < (int)vecs[1 - i]->size(); ++j)
                idCobo[1 - i].emplace_back(vecs[1 - i]->at(j));
        }

    morphisms.addCols(newCols);
    for (int i = 0; i < newRows; ++i) {
        morphisms.addRows();
        if (mats[0]) {
            constIterator_t j;
            for (j.setToRowBegin(mats[0]->morphisms, i / vecs[1]->size());
                    j.isOn(); j.stepAlongRow()) {
                LCCobos<cobordism_tpl> toAdd; 
                toAdd.setToUnion(domainVecs[0]->at(j.getCol()),
                        vecs[1]->at(i % vecs[1]->size()),
                        vecs[0]->at(i / vecs[1]->size()),
                        vecs[1]->at(i % vecs[1]->size()),
                        *(j.getVal()), idCobo[1].at(i % vecs[1]->size()));
                if (signSwitch)
                    toAdd.switchSign();
                const bool isInv = j.getVal()->isInvertible(
                        domainVecs[0]->at(j.getCol()),
                        vecs[0]->at(i / vecs[1]->size()));
                morphisms.setLastEntry(colShift + (j.getCol()) * vecs[1]->size()
                        + i % vecs[1]->size(), std::move(toAdd), isInv);
            }
        }
        if (mats[1]) {
            constIterator_t j;
            for (j.setToRowBegin(mats[1]->morphisms, i % vecs[1]->size());
                    j.isOn(); j.stepAlongRow()) {
                LCCobos<cobordism_tpl> toAdd; 
                toAdd.setToUnion(vecs[0]->at(i / vecs[1]->size()),
                        domainVecs[1]->at(j.getCol()),
                        vecs[0]->at(i / vecs[1]->size()),
                        vecs[1]->at(i % vecs[1]->size()),
                        idCobo[0].at(i / vecs[1]->size()), *(j.getVal()));
                const bool isInv =
                    j.getVal()->isInvertible(domainVecs[1]->at(j.getCol()),
                            vecs[1]->at(i % vecs[1]->size()));
                morphisms.setLastEntry(colShift + oldCols + j.getCol() +
                        mats[1]->getColCount() * (i / vecs[1]->size()),
                        std::move(toAdd), isInv);
            }
        }
    }
    assert(isSane());
    return oldCols;
}
