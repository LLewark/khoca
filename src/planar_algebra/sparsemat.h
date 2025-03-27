/*
 *
 *    src/planar_algebra/sparsemat.h --- Part of khoca, a knot homology calculator
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

template <class matrix_tpl, class entry_tpl> class GeneralIterator;
template <class entry_tpl> class SMIterator;
template <class entry_tpl> class SMconstIterator;
template <class entry_tpl> class ThreadSafeIterator;

template <class ctr_tpl, class index_tpl>
void makeSure(ctr_tpl &d, index_tpl a, bool add);

template <class entry_tpl>
class SparseMat {
    public:
        typedef typename std::vector<entry_tpl>::size_type idx_t;
        static_assert(sizeof(idx_t) >=
                sizeof(typename std::vector<idx_t>::size_type),
                "idx_t is too small.");
        typedef std::vector<idx_t> uintCont_t;

#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        static SparseMat<entry_tpl> setToDual(
                const SparseMat<entry_tpl> &other);
        /** Constructs an empty matrix with rows rows and cols columns.
         * @complexity O(1).
         */
        SparseMat(idx_t rows, idx_t cols);
        SparseMat() : SparseMat(0,0) { }
        explicit SparseMat(std::ifstream &f, bool intCoefficients);
        
        bool hasInvertibles() const;
        idx_t numberOfInvertibles() const;
        bool isZero() const;
        void writeToBin(std::ofstream &f) const;
        /** @brief Sets the entry at coordinate (row, col) to v (making a copy).
         * If a non-zero entry exists at that coordinate, it is lost
         * (and 0 is returned).
         * If v itself is zero, eraseEntry is called (and -1 is returned).
         * Otherwise, 1 is returned.
         * @pre row and col are valid indices.
         * @post v is set, invertibles is up to date.
         * @complexity O(log(entryCount) + rowCount).
         */
        int setEntry(idx_t row, idx_t col, entry_tpl v, bool isInv,
                idx_t *minIdx = nullptr);
        /** @brief High-speed version of setEntry: sets the last entry (highest
         * column index in last row). If isInvKnown, then the function does not
         * check itself whether v is invertible, but relies on isInv for that
         * information.
         * @pre v is non-zero; col is a valid index, and there are no elements
         *      at or after column col in the last row.
         * @post The entry is added, and invertibles is up to date.
         * @complexity O(1).
         */
        void setLastEntry(idx_t col, entry_tpl v, bool isInv);
        /** @brief Another version of setLastEntry: this one allows a different
         * row than the last; moreover, it accepts an rvalue for the entry,
         * and it check itself for invertibility.
         * @pre v is non-zero; col is a valid index, there is no element before
         * the position (row, col), and the last row has an entry.
         * @post The entry is added, and invertibles is up to date.
         * @complexity O(1).
         */
        void setLastEntry(idx_t col, idx_t row, entry_tpl &&v, bool isInv);
        /** Return entry at (row, col).
         * @pre row and col are valid indices.
         * @return a pointer to the value, nullptr if there is no entry at that
         *         position.
         * @post the matrix is unchanged.
         * @complexity O(colcount).
         */
        entry_tpl* getEntry(idx_t row, idx_t col, idx_t *minIdx = nullptr);
        const entry_tpl* getEntry(idx_t row, idx_t col, idx_t *minIdx = nullptr) const;
        /** Erases the entry at (row, col). If there is none, nothing happens.
         * @post No entry at (row, col), invertibles is up to date.
         * @complexity O(entryCount).
         */
        bool eraseEntry(idx_t row, idx_t col, idx_t *minIdx = nullptr);
        /** @complexity O(entryCount).
         */
        void eraseEntryByIdx(idx_t idx, idx_t row);
        /** Checks for sanity. This should always be true if called from the
         * outside (i.e. not from a member function of this).
         * @pre None.
         * @return Returns true if the matrix is consistent and
         * the deque of invertibles contains all and only the invertibles.
         */
        bool isSane() const;
        /** Erases row row.
         * @pre row is a valid row index.
         * @post the row row is deleted, invertibles is up to date.
         */
        void eraseRow(idx_t row);
        /** Erases column col.
         * @pre col is a valid column index.
         * @post the col column is deleted, invertibles is up to date.
         * @complexity O(entryCount).
         */
        void eraseCol(idx_t col);
        /** Returns whether the matrix is fully simplified, i.e. has no entries
         *  that are invertible.
         * @complexity O(entryCount).
         */
        bool isFullyGaussed() const { return invertibles.empty(); }
        /** Deletes all elements whose indices are in the delDeque;
         * delDeque is assumed to be sorted ascendingly (i.e. the smallest
         * index is at the beginning).
         */
        //void removeFromDeque(const std::deque<idx_t> &delDeque);
        /** Add n empty rows at the bottom of the matrix.
         * @pre None.
         * @post There are n more rows.
         * @complexity O(n).
         */
        void addRows(idx_t n = 1);
        /** Add n empty cols at the right end of the matrix.
         * @pre None.
         * @post There are n more cols.
         * @complexity O(1).
         */
        void addCols(idx_t n = 1);
        /** Appends copies copies of the row row at matrix' bottom.
         * @pre row is a valid row index.
         * @post copies copies are made, invertibles is up to date.
         * @complexity O(colCount * copies).
         */
        void copyRow(idx_t row, idx_t copies);
        /** Appends copies copies of the column col at the far right of the
         *  matrix
         * @pre row is a valid column index.
         * @post copies copies are made, invertibles is up to date.
         * @complexity O(entryCount + rowCount * copies).
         */
        void copyCol(idx_t row, idx_t copies);
        /** @complexity O(log(colCount)).
         */
        bool stepToNextInv(idx_t &invIdx, idx_t &row, idx_t &col,
                entry_tpl *&iso);

        void extractRow(idx_t row, idx_t avoidCol,
                std::vector<entry_tpl> &rVal, uintCont_t &rColInd) const;
        void extractCol(idx_t col, idx_t avoidRow,
                std::vector<entry_tpl> &cVal, uintCont_t &cRowInd) const;

        void setRowNumber(idx_t r);

        void noLongerInvertible(idx_t idx);
        void isNowInvertible(idx_t idx);

        const uintCont_t& getInvertibles() const { return invertibles; }
        
        idx_t getColCount() const;
        idx_t getRowCount() const;
        idx_t getEntryCount() const;

        void swap(SparseMat<entry_tpl> &other);
        void reserve(int sVal, int sRow, int sInv);

        std::ostream& detailedOutput(std::ostream &os) const;
    private:
        typedef typename std::vector<entry_tpl> valCont_t;

        idx_t colCount;
        /** Contains all entries ordered line-by-line, from left to right in
         * each line. "n-th entry" refers to this ordering.
         */
        valCont_t val;
        /** Saves in which column the n-th entry is. This vector has thus the
         * same size as val.
         */
        uintCont_t colInd;
        /* One entry more than number of rows. The r-th entry (for r = 0,1,...)
         * is the index of the first entry that is in a row with number at least
         * r. If no such * entry exists (and in particular for r = number of
         * rows), the r-th entry is 1 + maximal index.
         * Examples:
         * The empty matrix has rowPtr = {0}.
         * A matrix with r rows, c cols, but no entries has rowPtr consisting of
         * (r + 1) many zeroes.
         * The matrix
         * 0 0
         * * 0
         * 0 0
         * has rowPtr = {1, 1, 2, 2}.
         */
        uintCont_t rowPtr;
        /** Saves all indices of invertible entries. Size equals the number of
         * invertible entries.
         */
        uintCont_t invertibles;

        const entry_tpl* getEntry_(idx_t row, idx_t col, idx_t *minIdx = nullptr) const;

#ifndef notests
        void checkCol(idx_t col) const;
        void checkRow(idx_t row) const;
#endif
        /** @complexity O(log(colCount)).
         */
        idx_t goToCol(idx_t row, idx_t col, idx_t *minIdx = nullptr) const;
    friend class GeneralIterator<SparseMat<entry_tpl>, entry_tpl>;
    friend class GeneralIterator<const SparseMat<entry_tpl>, const entry_tpl>;
    friend class ThreadSafeIterator<entry_tpl>;
    friend class SMIterator<entry_tpl>;
    friend class SMconstIterator<entry_tpl>;
};

template <class entry_tpl>
std::ostream& operator<<(std::ostream &, const SparseMat<entry_tpl>&);

template <class matrix_tpl, class entry_tpl>
class GeneralIterator {
    public:
        typedef typename matrix_tpl::idx_t idx_t;
        void setToRowBegin(matrix_tpl &sm, idx_t row);
        void setToColBegin(matrix_tpl &sm, idx_t col);
        void setToMatBegin(matrix_tpl &sm);

        void stepAlongRow();
        void stepAlongCol();
        void stepAlongMat();

        bool isOn() const;

        void correctIdx(idx_t x);

        idx_t getIdx() const;
        idx_t getRow() const;
        idx_t getCol() const;
    protected:
        idx_t i, r;
        matrix_tpl *m;

        void stepAlongCol_(idx_t col);
};

template <class entry_tpl>
class SMIterator : public GeneralIterator<SparseMat<entry_tpl>, entry_tpl> {
    public:
        void stepAlongRow(bool deleteEntry = false);
        void stepAlongCol(bool deleteEntry = false);
        void stepAlongMat(bool deleteEntry = false);

        typedef typename SparseMat<entry_tpl>::idx_t idx_t;

        entry_tpl* getVal();

    private:
        typedef GeneralIterator<SparseMat<entry_tpl>, entry_tpl> base_t;
        using base_t::i;
        using base_t::r;
        using base_t::m;
        using base_t::stepAlongCol_;
};

template <class entry_tpl>
class SMconstIterator : public GeneralIterator<const SparseMat<entry_tpl>,
        const entry_tpl> { 
    public:
        const entry_tpl* getVal() const;
    private:
        typedef GeneralIterator<const SparseMat<entry_tpl>,
                const entry_tpl> base_t;
        using base_t::i;
        using base_t::r;
        using base_t::m;
        using base_t::stepAlongCol_;
};
