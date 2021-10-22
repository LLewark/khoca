/*
 *
 *    src/planar_algebra/sparsemat.cpp --- Part of khoca, a knot homology calculator
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
#include <deque>
#include <vector>
#include <bitset>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <iomanip>
#include <sstream>
#include <algorithm>
// Next line is necessary due to a bug of mpir
#include <stddef.h>
#include <gmp.h>
#include <fstream>

#include "../shared.h"

#include <assert.h>

#include "sparsemat.h"
#include "coefficient_rings.h"
#include "../planar_algebra/planar_algebra.h"
#include "../krasner/krasner.h"
#include "sparsematExplicitTemplates.cpp"

#ifndef getsize
template <class entry_tpl>
void SparseMat<entry_tpl>::printSize(std::vector<long long> &s) const {
    s.at(4) += sizeof(entry_tpl) * val.capacity();
    s.at(3) += sizeof(idx_t) *
        (colInd.capacity() + rowPtr.capacity() + invertibles.capacity());
    for (typename valCont_t::const_iterator i = val.begin();
            i != val.end(); ++i)
        i->printSize(s);
}
#endif

template <class entry_tpl>
std::ostream& SparseMat<entry_tpl>::detailedOutput(std::ostream &o) const {
    SMconstIterator<entry_tpl> i;
    i.setToMatBegin(*this);
    o << "[";
    if (getRowCount())
        for (idx_t i = 0; i < getRowCount(); ++i) {
            for (idx_t j = 0; j < getColCount(); ++j) {
                const entry_tpl* const x = getEntry(i, j);
                if (x)
                    x->detailedOutput(o);
                else
                    o << "0";
                if (j + 1 < getColCount())
                    o << ", ";
            }
            if (i + 1 < getRowCount())
                o << "; ";
        }
    return o << "]";
}

template <class ctr_tpl, class idx_tpl>
void makeSure(ctr_tpl &d, idx_tpl a, bool add) {
    const typename ctr_tpl::iterator it = std::lower_bound(
            d.begin(), d.end(), a);
    if (add) {
        if ((it == d.end()) || (*it != a))
            d.insert(it, a);
    } else {
        if ((it != d.end()) && (*it == a))
            d.erase(it, it + 1);
    }
}

template <class entry_tpl>
SparseMat<entry_tpl> SparseMat<entry_tpl>::setToDual(
        const SparseMat<entry_tpl> &other) {
    SparseMat<entry_tpl> dual;
    dual.rowPtr.resize(other.colCount + 1, 0);
    dual.colCount = other.rowPtr.size() - 1;
    for (auto i = other.colInd.begin(); i != other.colInd.end(); ++i)
        dual.rowPtr.at((*i) + 1) += 1;
    for (auto i = dual.rowPtr.begin() + 1; i != dual.rowPtr.end(); ++i)
        *i += *(i - 1);
    std::vector<idx_t> newIdxToOld(other.val.size());
    dual.colInd.resize(other.colInd.size());
    {
        std::vector<idx_t> countPerColumn(other.colCount, 0);
        auto row = other.rowPtr.begin();
        for (auto i = other.colInd.begin(); i != other.colInd.end(); ++i) {
            const idx_t newIdx = dual.rowPtr.at(*i) + countPerColumn.at(*i);
            const idx_t oldIdx = i - other.colInd.begin();
            newIdxToOld.at(newIdx) = oldIdx;
            row = std::upper_bound(row, other.rowPtr.end(), oldIdx);
            dual.colInd.at(newIdx) = row - other.rowPtr.begin() - 1;
            countPerColumn.at(*i) += 1;
        }
    }
    dual.val.reserve(other.val.size());
    std::vector<bool> isInvertible(other.val.size(), false);
    for (auto i = other.invertibles.begin(); i != other.invertibles.end(); ++i)
        isInvertible.at(*i) = true;
    dual.invertibles.reserve(other.invertibles.size());
    for (auto i = newIdxToOld.begin(); i != newIdxToOld.end(); ++i) {
        dual.val.push_back(other.val.at(*i));
        if (isInvertible.at(*i))
            dual.invertibles.push_back(i - newIdxToOld.begin());
    }
    return std::move(dual);
}

template <class entry_tpl>
bool SparseMat<entry_tpl>::hasInvertibles() const {
    return ! invertibles.empty();
}

template <class entry_tpl>
typename SparseMat<entry_tpl>::idx_t SparseMat<entry_tpl>::numberOfInvertibles()
        const {
    return invertibles.size();
}

template <class entry_tpl>
bool SparseMat<entry_tpl>::isZero() const {
    return val.empty();
}

template <class entry_tpl>
void SparseMat<entry_tpl>::noLongerInvertible(idx_t idx) {
    makeSure(invertibles, idx, false);
}

template <class entry_tpl>
void SparseMat<entry_tpl>::isNowInvertible(idx_t idx) {
    makeSure(invertibles, idx, true);
}

template <class ctr_tpl, class idx_tpl>
void modSortedContainer(
        ctr_tpl &d, idx_tpl a, idx_tpl b, int addRemove, int shift) {
    const typename ctr_tpl::iterator it = std::lower_bound(
            d.begin(), d.end(), a);
    const typename ctr_tpl::iterator it2 = std::lower_bound(it, d.end(), b);
    if (shift)
        for (typename ctr_tpl::iterator it3 = it; it3 != d.end(); ++it3)
            *it3 += (b - a) * shift;
    if (addRemove < 0)
        d.erase(it, it2);
    else if (addRemove > 0)
        d.insert(it, a);
}

template <class entry_tpl>
SparseMat<entry_tpl>::SparseMat(std::ifstream &f, bool intCoefficients) {
    uint32_t colCount_, colIndSize_, rowPtrSize_, invertiblesSize_;
    readFromBinTpl(f, colCount_);
    readFromBinTpl(f, colIndSize_);
    colCount = colCount_;
    for (uint32_t i = 0; i < colIndSize_; ++i)
        val.emplace_back(f, intCoefficients);
    for (uint32_t i = 0; i < colIndSize_; ++i) {
        uint32_t v;
        readFromBinTpl(f, v);
        colInd.push_back(std::move(v));
    }
    readFromBinTpl(f, rowPtrSize_);
    for (uint32_t i = 0; i < rowPtrSize_; ++i) {
        uint32_t v;
        readFromBinTpl(f, v);
        rowPtr.push_back(std::move(v));
    }
    readFromBinTpl(f, invertiblesSize_);
    for (uint32_t i = 0; i < invertiblesSize_; ++i) {
        uint32_t v;
        readFromBinTpl(f, v);
        invertibles.push_back(std::move(v));
    }
}

template <class entry_tpl>
void SparseMat<entry_tpl>::writeToBin(std::ofstream &f) const {
    writeToBinTpl(f, (uint32_t)colCount);
    writeToBinTpl(f, (uint32_t)colInd.size());
    for (typename valCont_t::const_iterator i = val.begin();
            i != val.end(); ++i)
        i->writeToBin(f);
    for (typename uintCont_t::const_iterator i = colInd.begin();
            i != colInd.end(); ++i)
        writeToBinTpl(f, (uint32_t)(*i));
    writeToBinTpl(f, (uint32_t)rowPtr.size());
    for (typename uintCont_t::const_iterator i = rowPtr.begin();
            i != rowPtr.end(); ++i)
        writeToBinTpl(f, (uint32_t)(*i));
    writeToBinTpl(f, (uint32_t)invertibles.size());
    for (typename uintCont_t::const_iterator i = invertibles.begin();
            i != invertibles.end(); ++i)
        writeToBinTpl(f, (uint32_t)(*i));
}

template <class entry_tpl>
SparseMat<entry_tpl>::SparseMat(idx_t rows, idx_t cols) : colCount(cols) {
    rowPtr.resize(rows + 1, 0);
}

template <class entry_tpl>
void SparseMat<entry_tpl>::setLastEntry(idx_t col, entry_tpl v, bool isInv) {
    assert(col < getColCount());
    assert(! v.isZero());
    assert((rowPtr.back() == rowPtr.at(rowPtr.size() - 2)) ||
            (colInd.empty()) || (colInd.back() < col));

    colInd.push_back(col);
    rowPtr.back() += 1;
    val.push_back(std::move(v));
    if (isInv)
        invertibles.push_back(val.size() - 1);
}

template <class entry_tpl>
void SparseMat<entry_tpl>::setLastEntry(idx_t col, idx_t row,
        entry_tpl &&v, bool isInv) {
    assert(col < getColCount());
    assert(! v.isZero());
    // The matrix is empty, or the last row is not empty
    assert((rowPtr.size() == 1) || (rowPtr.back() > rowPtr.at(rowPtr.size() - 2)));
    // The new entry is in a new row; or in the last row, to the right of the
    // last entry
    assert((row >= rowPtr.size() - 1) ||
        ((row == rowPtr.size() - 2) && (col > colInd.back())));
            
    colInd.push_back(col);
    val.push_back(v);    
    if (row >= rowPtr.size() - 1) {
        rowPtr.resize(row + 1, val.size() - 1);
        rowPtr.push_back(val.size());
    } else
        rowPtr.back() += 1;
    if (isInv)
        invertibles.push_back(val.size() - 1);
}

template <class entry_tpl>
int SparseMat<entry_tpl>::setEntry(
        idx_t row, idx_t col, entry_tpl v, bool isInv, idx_t *minIdx) {
    assert(row < getRowCount());
    assert(col < getColCount());
    if (v.isZero())
        return eraseEntry(row, col, minIdx) ? -1 : 0;

    const idx_t i = goToCol(row, col, minIdx);

    if ((i < rowPtr[row + 1]) && (colInd.at(i) == col)) {
        // Old entry is destroyed, copy of new entry is made
        val.at(i) = std::move(v);
        makeSure(invertibles, i, isInv);
        return 0;
    } else {
        // Copy of new entry is made
        val.insert(val.begin() + i, std::move(v));
        colInd.insert(colInd.begin() + i, col);
        for (idx_t j = row; j < getRowCount(); ++j)
            ++rowPtr[j + 1];
        if (isInv)
            modSortedContainer(invertibles, i, i + 1, 1, 1);
        else
            modSortedContainer(invertibles, i, i + 1, 0, 1);
        return 1;
    }
}

template <class entry_tpl>
typename SparseMat<entry_tpl>::idx_t SparseMat<entry_tpl>::goToCol(
        idx_t row, idx_t col, idx_t *minIdx) const {
    idx_t result = std::lower_bound(colInd.begin() +
            (minIdx ? std::max(rowPtr.at(row), *minIdx) : rowPtr.at(row)),
            colInd.begin() + rowPtr.at(row + 1), col) - colInd.begin();
    if (minIdx)
        *minIdx = result;
    return result;
}

template <class entry_tpl>
const entry_tpl* SparseMat<entry_tpl>::getEntry(idx_t row, idx_t col, idx_t *minIdx) const {
    return getEntry_(row, col, minIdx);
}

template <class entry_tpl>
entry_tpl* SparseMat<entry_tpl>::getEntry(idx_t row, idx_t col, idx_t *minIdx) {
    return const_cast<entry_tpl*>(getEntry_(row, col, minIdx));
}

template <class entry_tpl>
const entry_tpl* SparseMat<entry_tpl>::getEntry_(idx_t row, idx_t col, idx_t *minIdx) const {
    assert(row < getRowCount());
    assert(col < getColCount());

    const idx_t i = goToCol(row, col, minIdx);
        
    if ((i < rowPtr[row + 1]) && (colInd[i] == col))
        return &(val.at(i));
    else return nullptr;
}

template <class entry_tpl>
void SparseMat<entry_tpl>::eraseEntryByIdx(idx_t idx, idx_t row) {
    val.erase(val.begin() + idx);
    colInd.erase(colInd.begin() + idx);
    for (idx_t j = row + 1; j < rowPtr.size(); ++j)
        rowPtr.at(j) -= 1;
    modSortedContainer(invertibles, idx, idx + 1, -1, -1);
}

template <class entry_tpl>
bool SparseMat<entry_tpl>::eraseEntry(idx_t row, idx_t col, idx_t *minIdx) {
    assert(row < getRowCount());
    assert(col < getColCount());
    idx_t i = goToCol(row, col, minIdx);
    if ((i < rowPtr[row + 1]) && (colInd[i] == col)) {
        eraseEntryByIdx(i, row);
        return true;
    }
    return false;
}

template <class entry_tpl>
bool SparseMat<entry_tpl>::isSane() const {
    if ((val.size() != colInd.size()) || (rowPtr.size() == 0) ||
            (rowPtr.back() != val.size())) return false;
    for (idx_t r = 0; r < rowPtr.size() - 1; ++r) {
        if (rowPtr.at(r) > rowPtr.at(r + 1)) return false;
        idx_t oldCol = 0;
        for (idx_t i = rowPtr.at(r); i < rowPtr.at(r + 1); ++i) {
            if ((oldCol > colInd.at(i)) || (colInd.at(i) >= colCount)) {
		insane();
                return false;
	    }
            oldCol = colInd.at(i) + 1;
        }
    }

    for (idx_t i = 0; i < invertibles.size(); ++i)
        if (val.at(invertibles.at(i)).obviouslyNotInvertible()) {
	    insane();
            return false;
	}
    return true;
}

template <class entry_tpl>
void SparseMat<entry_tpl>::eraseRow(idx_t row) {
    assert(row < getRowCount());
    const idx_t s = rowPtr.at(row);
    const idx_t n = rowPtr.at(row + 1) - rowPtr.at(row);
    for (idx_t i = row + 1; i < rowPtr.size(); ++i) 
        rowPtr.at(i) -= n;
    rowPtr.erase(rowPtr.begin() + row);
    val.erase(val.begin() + s, val.begin() + s + n);
    colInd.erase(colInd.begin() + s, colInd.begin() + s + n);
    modSortedContainer(invertibles, s, s + n, -1, -1);
}

template <class entry_tpl>
void SparseMat<entry_tpl>::eraseCol(idx_t col) {
    assert(col < getColCount());
    idx_t soFarDeleted = 0;
    for (idx_t r = 0; r < rowPtr.size() - 1; ++r) {
        rowPtr.at(r + 1) -= soFarDeleted;
        idx_t i = goToCol(r, col);
        if ((i < rowPtr.at(r + 1)) && (colInd.at(i) == col)) {
            modSortedContainer(invertibles, i, i + 1, -1, -1);
            colInd.erase(colInd.begin() + i);
            val.erase(val.begin() + i);
            soFarDeleted += 1;
            rowPtr.at(r + 1) -= 1;
        }
        while (i < rowPtr.at(r + 1))
            colInd.at(i++) -= 1;
    }
    colCount -= 1;
}

template <class entry_tpl>
bool SparseMat<entry_tpl>::stepToNextInv(
        idx_t &invIdx, idx_t &row, idx_t &col, entry_tpl *&iso) {
    invIdx += 1;
    if (invIdx >= invertibles.size())
        return false;
    col = colInd.at(invertibles.at(invIdx));

    row = std::upper_bound(rowPtr.begin() + row, rowPtr.end(), 
            invertibles.at(invIdx)) - rowPtr.begin() - 1;

    iso = &(val.at(invertibles.at(invIdx)));
    return true;
}

template <class entry_tpl>
void SparseMat<entry_tpl>::extractRow(idx_t row, idx_t avoidCol,
        valCont_t &rVal, uintCont_t &rColInd) const {
    assert((row >= 0) && (row <= rowPtr.size() - 2));

    auto i = std::lower_bound(
            colInd.begin() + rowPtr.at(row),
            colInd.begin() + rowPtr.at(row + 1), avoidCol);
    if (*i == avoidCol) {
        rVal.insert(rVal.end(), val.begin() + rowPtr.at(row),
                val.begin() + (i - colInd.begin()));
        rVal.insert(rVal.end(), val.begin() + (i - colInd.begin()) + 1,
                val.begin() + rowPtr.at(row + 1));
        rColInd.insert(rColInd.end(), colInd.begin() + rowPtr.at(row), i);
        rColInd.insert(rColInd.end(), i + 1,
                colInd.begin() + rowPtr.at(row + 1));
    } else {
        rVal.insert(rVal.end(), val.begin() + rowPtr.at(row),
                val.begin() + rowPtr.at(row + 1));
        rColInd.insert(rColInd.end(), colInd.begin() + rowPtr.at(row),
                colInd.begin() + rowPtr.at(row + 1));
    }
}

template <class entry_tpl>
void SparseMat<entry_tpl>::extractCol(idx_t col, idx_t avoidRow,
        std::vector<entry_tpl> &cVal, uintCont_t &cRowInd) const {
    assert((col >= 0) && (col < colCount));
    SMconstIterator<entry_tpl> sameCol;
    for (sameCol.setToColBegin(*this, col); sameCol.isOn();
            sameCol.stepAlongCol()) {
        idx_t r = sameCol.getRow();
        if (r != avoidRow) {
            cVal.push_back(*sameCol.getVal());
            cRowInd.push_back(r);
        }
    }
}

template <class entry_tpl>
void SparseMat<entry_tpl>::setRowNumber(idx_t r) {
    rowPtr.resize(r + 1, rowPtr.back());
}

template <class entry_tpl>
typename SparseMat<entry_tpl>::idx_t SparseMat<entry_tpl>::getColCount() const {
    return colCount;
} 

template <class entry_tpl>
typename SparseMat<entry_tpl>::idx_t SparseMat<entry_tpl>::getRowCount() const {
    return rowPtr.size() - 1;
}

template <class entry_tpl>
typename SparseMat<entry_tpl>::idx_t SparseMat<entry_tpl>::getEntryCount()
    const {
    return val.size();
}

template <class entry_tpl>
void SparseMat<entry_tpl>::swap(SparseMat<entry_tpl> &other) {
    std::swap(colCount, other.colCount);
    val.swap(other.val);
    colInd.swap(other.colInd);
    rowPtr.swap(other.rowPtr);
    invertibles.swap(other.invertibles);
}

template <class entry_tpl>
void SparseMat<entry_tpl>::reserve(int sVal, int sRow, int sInv) {
    if (sVal) {
        val.reserve(sVal);
        colInd.reserve(sVal);
    }
    if (sRow)
        rowPtr.reserve(sRow + 1);
    if (sInv)
        invertibles.reserve(sInv);
}

template <class entry_tpl>
void SparseMat<entry_tpl>::addRows(idx_t n) {
    rowPtr.resize(rowPtr.size() + n, rowPtr.back());
}

template <class entry_tpl>
void SparseMat<entry_tpl>::addCols(idx_t n) {
    colCount += n;
}

template <class entry_tpl>
void SparseMat<entry_tpl>::copyRow(idx_t row, idx_t copies) {
    assert(row < getRowCount());
    assert(isSane());
    const idx_t nbEntries = rowPtr.at(row + 1) - rowPtr.at(row);
    const idx_t totalNbEntries = rowPtr.back();

    const typename uintCont_t::iterator invBegin = std::lower_bound(
            invertibles.begin(), invertibles.end(), rowPtr.at(row));
    const typename uintCont_t::iterator invEnd = std::lower_bound(
            invBegin, invertibles.end(), rowPtr.at(row + 1));

    val.reserve(val.size() + nbEntries * copies);
    colInd.reserve(colInd.size() + nbEntries * copies);
    rowPtr.reserve(rowPtr.size() + copies);
    invertibles.reserve(invertibles.size() + (invEnd - invBegin) * copies);

    for (idx_t i = 1; i <= copies; ++i) {
        std::copy_n(val.begin() + rowPtr.at(row), nbEntries,
                std::back_inserter(val));
        std::copy_n(colInd.begin() + rowPtr.at(row), nbEntries,
                std::back_inserter(colInd));

        for (typename uintCont_t::const_iterator j = invBegin; j != invEnd; ++j)
            invertibles.push_back(*j - rowPtr.at(row) + rowPtr.back());

        rowPtr.push_back(totalNbEntries + i * nbEntries);
    }
    assert(isSane());
}

template <class entry_tpl>
void SparseMat<entry_tpl>::copyCol(idx_t col, idx_t copies) {
    assert(col < getColCount());
    idx_t nbEntries = 0;
    idx_t invIt = 0;
    for (idx_t r = 0; r < getRowCount(); ++r) {
        bool originalIsInvertible = false;
        while ((invIt < invertibles.size()) &&
                (invertibles.at(invIt) < rowPtr.at(r + 1))) {
            invertibles.at(invIt) += nbEntries * copies;
            if (colInd.at(invertibles.at(invIt)) == col)
                originalIsInvertible = true;
            ++invIt;
        }
        rowPtr.at(r + 1) += nbEntries * copies;

        idx_t i = goToCol(r, col);
        if ((i < rowPtr.at(r + 1)) && (colInd.at(i) == col)) {
            ++nbEntries;
            val.insert(val.begin() + rowPtr.at(r + 1), copies, val.at(i));
            for (idx_t k = copies; k > 0; --k)
                colInd.insert(
                        colInd.begin() + rowPtr.at(r + 1), k + colCount - 1);
            if (originalIsInvertible) {
                invertibles.reserve(invertibles.size() + copies);
                for (idx_t k = copies; k > 0; --k)
                    invertibles.insert(invertibles.begin() + invIt,
                            rowPtr.at(r + 1) + k - 1);
                invIt += copies;
            }
            rowPtr.at(r + 1) += copies;
        }

    }
    colCount += copies;
}

template <class entry_tpl>
void SMIterator<entry_tpl>::stepAlongRow(bool deleteEntry) {
    assert(m);
    if (deleteEntry)
        m->eraseEntryByIdx(i, r);
    else
        ++(i);
    if (i >= m->rowPtr.at(r + 1))
        m = 0;
}

template <class entry_tpl>
void SMIterator<entry_tpl>::stepAlongCol(bool deleteEntry) {
    assert(m);
    const idx_t c = m->colInd.at(i);
    if (deleteEntry)
        m->eraseEntryByIdx(i, r);
    stepAlongCol_(c);
}

template <class entry_tpl>
void SMIterator<entry_tpl>::stepAlongMat(bool deleteEntry) {
    assert(m);
    if (deleteEntry)
        m->eraseEntryByIdx(i, r);
    else
        ++(i);
    if (i >= m->rowPtr.back())
        m = 0;
    else
        r = std::upper_bound(m->rowPtr.begin() + r + 1, m->rowPtr.end(), i)
            - m->rowPtr.begin() - 1;
}

template <class entry_tpl>
entry_tpl* SMIterator<entry_tpl>::getVal() {
    return &(m->val.at(i));
}

template <class matrix_tpl, class entry_tpl>
void GeneralIterator<matrix_tpl, entry_tpl>::setToRowBegin(
        matrix_tpl &sm, idx_t row) {
    if (sm.rowPtr.at(row + 1) > sm.rowPtr.at(row)) {
        m = &sm;
        i = m->rowPtr.at(row);
        r = row;
    } else m = 0;
}

template <class matrix_tpl, class entry_tpl>
void GeneralIterator<matrix_tpl, entry_tpl>::setToColBegin(
        matrix_tpl &sm, idx_t col) {
    m = &sm;
    r = -1;
    stepAlongCol_(col);
}

template <class matrix_tpl, class entry_tpl>
void GeneralIterator<matrix_tpl, entry_tpl>::setToMatBegin(matrix_tpl &sm) {
    m = &sm;
    i = -1;
    r = 0;
    stepAlongMat();
}

template <class matrix_tpl, class entry_tpl>
void GeneralIterator<matrix_tpl, entry_tpl>::stepAlongRow() {
    if (!m) {
        std::cerr << "Row-stepping of non-active iterator.\n";
        throw;
    }
    ++i;
    if (i >= m->rowPtr.at(r + 1)) m = 0;
}

template <class matrix_tpl, class entry_tpl>
void GeneralIterator<matrix_tpl, entry_tpl>::stepAlongCol() {
    if (!m) {
        std::cerr << "Col-stepping of non-active iterator.\n";
        throw;
    }
    stepAlongCol_(getCol());
}

template <class matrix_tpl, class entry_tpl>
void GeneralIterator<matrix_tpl, entry_tpl>::stepAlongMat() {
    if (!m) {
        std::cerr << "Mat-stepping of non-active iterator.\n";
        throw;
    }
    ++i;
    if (i >= m->rowPtr.back())
        m = 0;
    else while (i >= m->rowPtr.at(r + 1))
        ++r;
}

template <class matrix_tpl, class entry_tpl>
bool GeneralIterator<matrix_tpl, entry_tpl>::isOn() const {
    return m;
}

template <class matrix_tpl, class entry_tpl>
void GeneralIterator<matrix_tpl, entry_tpl>::correctIdx(idx_t x) {
    i = x;
}

template <class matrix_tpl, class entry_tpl>
typename GeneralIterator<matrix_tpl, entry_tpl>::idx_t
        GeneralIterator<matrix_tpl, entry_tpl>::getIdx() const {
    return i;
}

template <class matrix_tpl, class entry_tpl>
typename GeneralIterator<matrix_tpl, entry_tpl>::idx_t
        GeneralIterator<matrix_tpl, entry_tpl>::getRow() const {
    return r;
}

template <class matrix_tpl, class entry_tpl>
typename GeneralIterator<matrix_tpl, entry_tpl>::idx_t
        GeneralIterator<matrix_tpl, entry_tpl>::getCol() const {
    return m->colInd.at(i);
}

template <class matrix_tpl, class entry_tpl>
void GeneralIterator<matrix_tpl, entry_tpl>::stepAlongCol_(idx_t col) {
    for (r = r + 1; r < m->rowPtr.size() - 1; ++r) {
        i = m->goToCol(r, col);
        if ((i < m->rowPtr.at(r + 1)) && (m->colInd.at(i) == col))
            return;
    }
    m = 0;
}

template <class entry_tpl>
entry_tpl const* SMconstIterator<entry_tpl>::getVal() const {
    assert(m);
    return &(m->val.at(i));
}
