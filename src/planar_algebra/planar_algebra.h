/*
 *
 *    src/planar_algebra/planar_algebra.h --- Part of khoca, a knot homology calculator
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

/*! \mainpage plakho
 *
 * \section intro_sec Mathematical Background
 * \subsection can_subsec Canopolis
 * A canopolis is a collection of objects and functions satisfying some axioms.
 * The objects are:
 * <ol>
 * <li>A 0-cell is just a natural number, thought of as a finite tuple of points
 * living on the boundary circle of a disc.</li>
 *
 * <li>For each 0-cell, there is a collection of 1-cells with that 0-cell as
 * (so-called) boundary; a 1-cell is thought of as some 1-dimensional geometric
 * object living in the disc, such as a 1-manifold or a graph, whose
 * intersection with the circle has the 0-cell as cardinality.</li>
 *
 * <li>For each pair of 1-cells with the same boundary, there is a collection of
 * 2-cells from the first to the second 1-cell; these 1-cells are called domain
 * and codomain of the 2-cell. A 2-cell is thought of as some 2-dimensional
 * geometric object living in a cylinder, such that the intersection of the
 * 2-cell with the bottom and top disc is the first and second 1-cell,
 * respectively, and the intersection with the side of the cylinder consists of
 * vertical lines.</li>
 *
 * <li>Every 1-cell has a identity 2-cell, which is thought of as the 1-cell
 * times an interval.</li>
 * </ol>
 *
 * These objects come with the following structure, i.e. functions:
 * <ol>
 * <li>There is the disjoint union of two 1-cells, which is another 1-cell,
 * whose boundary is the sum of the boundaries of the two 1-cells. This
 * operation is not commutative, because 0-cells are tuples, not sets, of
 * points, and in the boundary of the disjoint union, the boundary of the first
 * 1-cell comes first; but it is associative.</li>
 *
 * <li>There is the disjoint union of two 2-cells, which is a 2-cell from the
 * disjoint union of the domains of the two 2-cells to the disjoint union of
 * their codomains. This operation is associative, and the disjoint union of two
 * identities is an identity 2-cell.</li>
 * </ol>
 *
 *
 */

template <class cobordism_tpl> class Complex;
template <class tangle_tpl> class VecTangles;
template <class cobordism_tpl> class MatLCCobos;
template <class tangle_tpl> class Tangle;
template <class cobordism_tpl> class LCCobos;
class Boundary;

/** \class Tangle
 * Virtual class, not to be used itself, only to be inherited from.
 * Contains all methods a <i>tangle class</i> must support, except for:
 * <ol>
 * <li> friend std::ostream& operator<<(std::ostream& os, const tangle_tpl&);
 * Prints to the given os, returns os.
 * </li>
 * </ol>
 */
template <class tangle_tpl>
class Tangle {
    public:
#ifndef getsize
        virtual void printSize(std::vector<word64> &s) const = 0;
#endif
        /** @brief Sets this tangle to the disjoint union of the two tangles o1
         * and o2.
         * @pre None.
         * @post The old value of this is lost.
         */
        virtual void setToUnion(const tangle_tpl &o1, const tangle_tpl &o2) = 0;
        /** @brief Glues together two boundary points (i.e. connects them by an
         * arc using an input diagram of the planar algebra.
         * @pre gluePoint[0] and gluePoint[1] must be two different numbers
         * which are smaller than the boundary size of this tangle.
         * @post qShift is unchanged, the tangle now contains the gluing.
         */
        virtual void glue(const boundary_t gluePoint[2]) = 0;

        static tangle_tpl setToDual(const tangle_tpl &other);

        /** @brief "Removes" the first loop (loops have an order), i.e.
         * calculates the sum of tangles which is isomorphic to this one.
         * @pre This tangle contains a loop, newCopies is empty.
         * @post This is now the first summand, the other summands are added to
         * newCopies.
         */
        virtual void deloop(std::vector<tangle_tpl> &newCopies) = 0;
        /** @brief Returns whether this tangle has a loop.
         * @pre None.
         */
        virtual bool hasLoop() const = 0;

        /** @brief Reads this tangle from a string (using some format specified
         * by the inheriting tangle class).
         * @param i An iterator pointing to the position in some string where
         * reading starts.
         * @param end Read until here (exclusively).
         * @param boundarySize The number of boundary points of the tangle.
         * @pre The string portion between i and end a contains a correctly
         * formatted tangle.
         * @return Returns whether the preconditions were satisfied.
         * @post The old value of this is lost. i points to end (in case of
         * success) or to the character where the error occurred.
         */
        /*virtual bool readFromString(std::string::const_iterator &i,
                std::string::const_iterator end, boundary_t boundarySize) = 0;*/

        virtual void writeToBin(std::ofstream &f) const = 0;
        virtual std::ostream& detailedOutput(std::ostream &os) const = 0;


        qShift_t getQShift() const { return qShift; }
#ifndef notests
        virtual bool isSane(const Boundary *b = nullptr) const = 0;
#endif
    protected:
        /** The quantum shift.*/
        qShift_t qShift;
};

/** \class Cobo
 * Virtual class, not to be used itself, only to be inherited from.
 * Contains all methods a <i>cobordism class</i> must support, except for:
 * <ol>
 * <li> friend std::ostream& operator<<(std::ostream& os, const cobordism_tpl&);
 * Prints to the given os, returns os.
 * </li>
 * <li> cobordism_tpl(const cobordism_tpl &idTangle);
 * Constructs the identity cobordism to the given Tangle.
 * </li>
 * </ol>
 */
template <class tangle_tpl, class cobordism_tpl, class coeff_tpl>
class Cobo {
    public:
#ifndef getsize
        virtual void printSize(std::vector<word64> &s) const = 0;
#endif
        /**
         * Normally, a cobordism lives in a linear combination of cobordisms
         * (LCCobos); coefficient is its coefficient.
         */
        coeff_tpl coefficient;

        Cobo(int x) : coefficient(x) { }
        Cobo(): Cobo(1) { }

        void switchSign() { coefficient.switchSign(); }
        void setToNegInv() { coefficient.switchSign(); coefficient.inv(); }

        virtual int reducify() = 0;

        virtual qShift_t degree(boundary_t) const = 0;
        virtual qShift_t relativeDegree() const = 0;
        virtual void compose(const cobordism_tpl &other,
                std::vector<cobordism_tpl> &result, const tangle_tpl &lower,
                const tangle_tpl &middle, const tangle_tpl &upper) const = 0;
        virtual void setToUnion(const tangle_tpl& lowerLeft,
                const tangle_tpl& lowerRight, const tangle_tpl& upperLeft,
                const tangle_tpl& upperRight, const cobordism_tpl &o1,
                const cobordism_tpl &o2) = 0;
        virtual void glue(const tangle_tpl& lower, const tangle_tpl& upper,
                const boundary_t gluePoint[2], boundary_t boundarySize,
                LCCobos<cobordism_tpl> & father) = 0;
        virtual void modifyDeloopCopy(int kind, bool left,
                std::vector<cobordism_tpl> &v, const tangle_tpl &lower,
                const tangle_tpl &upper) = 0;
        virtual bool operator<(const cobordism_tpl &other) const = 0;
        virtual bool operator==(const cobordism_tpl &other) const = 0;
        /** Returns true if the cobordism is an isomorphism
         * (but doesn't care about the coefficient)
         */
        virtual bool isInvertible(const tangle_tpl &lowerTangle,
                const tangle_tpl &upperTangle) const = 0;
        /** Returns whether the cobordism is empty (i.e. a cobordism from and to
         * the empty tangle).
         */
        virtual bool isEmpty() const = 0;
        virtual void writeToBin(std::ofstream &f) const = 0;
        virtual std::ostream& detailedOutput(std::ostream &os) const = 0;
#ifndef notests
        virtual bool isSane(const Boundary *b = nullptr,
                const tangle_tpl *domain = nullptr,
                const tangle_tpl *codomain = nullptr) const = 0;
#endif
};

/** Don't inherit from this class.
 */
class Boundary {
    public:
        Boundary() : size_(0) { }
        Boundary(boundary_t x) : size_(x) { }
        boundary_t size() const { return size_; };
        void setToSum(const Boundary &b1, const Boundary &b2);
        void setToSize(boundary_t val) { size_ = val; }
        void glue() { size_ -=2 ; }
    private:
        /* Orientations are not needed, nor planarity. */
        boundary_t size_;
};


class AbstractComplex {
    public:
        /** Not sure this is the right way to do it. But using
         * virtual ~AbstractComplex();
         * or
         * virtual ~AbstractComplex() = 0;
         * instead gives "undefined symbol" errors when executing.
         */
        virtual ~AbstractComplex() {};

        virtual void glue(const boundary_t gluePoints[2]) = 0;
        virtual void deleteNonIsos() = 0;
        virtual void reducify(int root) = 0;
        virtual int simplifyOnce(int qDiff, int numThreads, int progress) = 0;
        virtual int simplifyOnceAtTDegree(int qDiff, int tDegree,
                const std::string* s, int numThreads, int progress) = 0;
        virtual int isSimplified(int qDiff, int tDegree) = 0;
        virtual int isSimplified(int qDiff) = 0;
        virtual void writeToBin(std::ofstream &f) const = 0;
        virtual int size() const = 0;

        virtual void initialiseFrobenius(const std::vector<int> &F, int N) const = 0;
        virtual void printFrobenius(std::ostream& s) const = 0;

        virtual AbstractComplex* setToDual(const AbstractComplex *other) const
            = 0;
        virtual AbstractComplex* copy(const AbstractComplex *other) const = 0;
        virtual AbstractComplex* loadFromFile(std::ifstream &f) const = 0;
        virtual AbstractComplex* tensor(const AbstractComplex *cc1,
                const AbstractComplex *cc2) const = 0;
        virtual void calculateSmith(std::ostream& s, int progress) const = 0;
#ifndef getsize
        virtual void printSize(std::vector<word64> &s) const = 0;
#endif
        virtual std::ostream& output(std::ostream &os) const = 0;
        virtual std::ostream& detailedOutput(std::ostream &os) const = 0;
        friend inline std::ostream& operator<<(std::ostream &os,
                const AbstractComplex& c) { return c.output(os); }

        /** careful, modification is not threadsafe!
          */
        int simplificationsCount;
};

template <class cobordism_tpl>
void calculateSmithFriend(const Complex<cobordism_tpl> &, std::ostream&, int progress);

/*class alphTangle;
class alphCobo;*/
/** Don't inherit from this class.
 */
template <class cobordism_tpl>
class Complex : public AbstractComplex {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        typedef typename cobordism_tpl::tangle_t tangle_t;
        typedef VecTangles<tangle_t> vecTangles_t;
        typedef MatLCCobos<cobordism_tpl> matLCCobos_t;
        typedef std::vector<vecTangles_t> vecTanglesCont_t;
        typedef std::vector<matLCCobos_t> matLCCobosCont_t;
        typedef typename SparseMat<LCCobos<cobordism_tpl> >::idx_t idx_t;
        
        AbstractComplex* setToDual(const AbstractComplex *other) const {
            return setToDualConcrete(*((Complex<cobordism_tpl>*)other));
        }
        AbstractComplex* copy(const AbstractComplex *other) const {
            return new Complex(*((Complex<cobordism_tpl>*)other));
        }
        static Complex<cobordism_tpl>* setToDualConcrete(
                const Complex<cobordism_tpl> &other);
        AbstractComplex* loadFromFile(std::ifstream &f) const {
            return new Complex<cobordism_tpl>(f); }
        AbstractComplex* tensor(const AbstractComplex *cc1,
                const AbstractComplex *cc2) const {
            return new Complex<cobordism_tpl>(*((Complex<cobordism_tpl>*)cc1),
                    *((Complex<cobordism_tpl>*)cc2));
        }

        /* Applies a Gaussian elimination, if this is possible;
         * otherwise deloops, if this is possible.
         * Returns 0 for "did simplify", 1 for "did not, but there are non-zero
         * differentials" and 2 for "did not, and all matrices are zero."
         */
        Complex() : globalTShift(0) { }
        Complex(boundary_t b, qShift_t t) : boundary(b), globalTShift(t) { }
        int simplifyOnce(int qDiff, int numThreads, int progress);
        int simplifyOnceAtTDegree(int qDiff, int tDegree, const std::string* s,
                int numThreads, int progress);
        int isSimplified(int qDiff, int tDegree);
        int isSimplified(int qDiff);

        bool tryToGauss(int i, int qDiff, int numThreads);
        bool tryToDeloop(int i);
        void showProgressBar(const std::string* s);
        /**
         * Format: 1. global t-shift, 2. boundary size, 3. read vectors (until ;
         * - i.e. a double semicolon, since one semicolon only ends the vector
         *   reading), 4. read matrices (until end)
         */
        Complex(const Complex &cc1, const Complex &cc2);
        /** @brief Glues together two boundary points (i.e. connects them by an
         * arc using an input diagram of the planar algebra.
         * @pre gluePoint[0] and gluePoint[1] must be two different numbers
         * which are smaller than the boundary size of this tangle.
         * @post Boundary is smaller by two (boundary points are renumbered
         * automatically), all matrices and vectors have been changed.
         */
        void glue(const boundary_t gluePoints[2]);
        void writeToBin(std::ofstream &f) const;
        explicit Complex(std::ifstream &f);

        void initialiseFrobenius(const std::vector<int> &F, int N) const {
            cobordism_tpl::initialiseFrobenius(F, N);
        }
        void printFrobenius(std::ostream& s) const {
            const int N = cobordism_tpl::frobenius.size() - 1;
            cobordism_tpl::coeff_t::printRing(N, s);
            s << "[X] / (";
            for (auto i = cobordism_tpl::frobenius.rbegin(); i != cobordism_tpl::frobenius.rend(); ++i) {
                if (i->isNonZero()) {
                    if (i != cobordism_tpl::frobenius.rbegin())
                        s << " + ";
                    s << *i;
                    int exp = N - (i - cobordism_tpl::frobenius.rbegin());
                    if (exp > 1)
                        s << "*X^" << exp;
                    else if (exp == 1)
                        s << "*X";
                }
            }
            s << ")";
        }

        void deleteNonIsos();
        void reducify(int root);
        void calculateSmith(std::ostream& s, int progress) const {
            calculateSmithFriend(*this, s, progress);
	}
	friend void calculateSmithFriend<cobordism_tpl>(
                const Complex &that, std::ostream& s, int progress);

        int size() const { return vecTangles.size(); }

#ifndef notests
        bool isSane() const;
#endif
        std::ostream& output(std::ostream &os) const {
            os << "[";
            const char* comma = "";
            for (auto i = vecTangles.begin(); i != vecTangles.end(); ++i) {
                if (! i->size())
                    continue;
                os << comma;
                i->printFormatted(os, globalTShift + (i - vecTangles.begin()));
                comma = ",";
            }
            os << "]";
            return os;
        }
        std::ostream& detailedOutput(std::ostream &os) const {
            for (auto i = vecTangles.begin(); (i + 1) != vecTangles.end(); ++i) {
                matLCCobos.at(i - vecTangles.begin()).detailedOutput(os);
                os << std::endl;
            }
            return os;
        }
    private:
        /* Essential properties. */
        vecTanglesCont_t vecTangles;
        matLCCobosCont_t matLCCobos;
        Boundary boundary;
        qShift_t globalTShift;
};

/** Don't inherit from this class.
 */
template <class tangle_tpl>
class VecTangles {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        typedef typename std::vector<tangle_tpl> tangleCont_t;

        static VecTangles<tangle_tpl> setToDual(
                const VecTangles<tangle_tpl> &other);

        VecTangles() { }
        /** ...
         * @pre The complex is fully simplified, in particular, the deloopstack
         * is empty.
         * @post The deloopstack is still empty.
         */
        void appendTensorProduct(const VecTangles<tangle_tpl> &s1,
                const VecTangles<tangle_tpl> &s2);
        void glue(const boundary_t gluePoints[2]);
        int deloopsToBeDone() const;
        bool deloopingDone() const;
        word64 simplifyOnce();
        /**
         * Format: Series of tangles separated by semicola.
         */
        void erase(word64 idx) { tangles.erase(tangles.begin() + idx,
                tangles.begin() + idx + 1);
            shiftWhatIsHigher(deloopStack, idx); assert(isSane()); }

        const tangle_tpl& at(word64 idx) const { return tangles.at(idx); }
        word64 size() const { return tangles.size(); }
        void writeToBin(std::ofstream &f) const;
        VecTangles(std::ifstream &f, boundary_t size);
#ifndef notests
        bool isSane(const Boundary *b = nullptr) const;
#endif
        void printFormatted(std::ostream& os, qShift_t tShift) const {
            const char* comma = "";
            for (typename tangleCont_t::const_iterator i = tangles.begin();
                    i != tangles.end(); ++i) {
                os << comma << "[" << tShift << "," << i->getQShift()
                    << ",0,1]";
                comma = ",";
            }
        }
        std::ostream& detailedOutput(std::ostream &os) const {
            os << "(";
            for (typename tangleCont_t::const_iterator i = tangles.begin();
                    i != tangles.end(); ++i) {
                i->detailedOutput(os);
                if (i + 1 != tangles.end())
                    os << ",";
            }
            os << ")";
            return os;
        }
        const tangleCont_t& getTangles() const { return tangles; } 
        void getQs(qShift_t &qMax, qShift_t &qMin,
                std::vector<int> &idxTranslator,
                std::vector<int> &qDims) const;
    private:
        tangleCont_t tangles;

        /** Contains a list of all indices of tangles which have a loop.
         * The stack is correct and complete at all times. */
        std::vector<word64> deloopStack;
};

/** Don't inherit from this class.
 */
template <class cobordism_tpl>
class MatLCCobos {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        static MatLCCobos<cobordism_tpl> setToDual(
                const MatLCCobos<cobordism_tpl> &other);

        MatLCCobos() { }
        MatLCCobos(int rows, int cols) : morphisms(rows, cols) { }

        typedef typename cobordism_tpl::coeff_t coeff_t;
        typedef typename cobordism_tpl::tangle_t tangle_t;
        typedef LCCobos<cobordism_tpl> LCCobos_t;
        typedef SMIterator<LCCobos_t> iterator_t;
        typedef ThreadSafeIterator<LCCobos_t> threadSafeIterator_t;
        typedef SMconstIterator<LCCobos_t> constIterator_t;
        typedef std::vector<tangle_t> tangleCont_t;
        typedef VecTangles<tangle_t> vecTangles_t;
        typedef typename SparseMat<LCCobos<cobordism_tpl> >::idx_t idx_t;

        void isNowInvertible(typename SparseMat<LCCobos<cobordism_tpl> >::idx_t
                idx) { morphisms.isNowInvertible(idx); }

        bool isZero() const { return morphisms.isZero(); }
        bool hasInvertibles() const { return morphisms.hasInvertibles(); }
        int numberOfInvertibles() const {
            return morphisms.numberOfInvertibles();
        }
        iterator_t getIterator() {
            iterator_t x;
            x.setToMatBegin(this->morphisms);
            return x;
        }
        constIterator_t getConstIterator() const {
            constIterator_t x;
            x.setToMatBegin(this->morphisms);
            return x;
        }
        void glue(const boundary_t gluePoints[2],
                const VecTangles<typename cobordism_tpl::tangle_t> &lowerVec,
                const VecTangles<typename cobordism_tpl::tangle_t> &upperVec,
                boundary_t boundarySize);
        /** Somehow responsible for tensor product of complexes.
         * ...
         * @pre The complex is fully simplified, in particular, there are no
         * invertibles in the matrix.
         * @post There are no invertibles in the matrix.
         */
        int hossa(int colShift, const VecTangles<typename
                cobordism_tpl::tangle_t> *domainVecs[2],
                const VecTangles<typename cobordism_tpl::tangle_t> *vecs[2],
                const MatLCCobos<cobordism_tpl> *mats[2], bool switchSign);
        int getColCount() const { return morphisms.getColCount(); }
        void addCols(int t) { morphisms.addCols(t); assert(isSane()); }
        void eraseCol(int c) { morphisms.eraseCol(c); assert(isSane()); }
        void eraseRow(int r) {
            morphisms.eraseRow(r);
            assert(isSane());
        }
        LCCobos_t* gaussianEliminationDone(idx_t &row, idx_t &col,
                const vecTangles_t &domain, const vecTangles_t &codomain,
                int qDiff);
        bool gaussianElimination(idx_t &row, idx_t &col,
                const vecTangles_t &domain, const vecTangles_t &codomain,
                int qDiff, int numThreads);
        void gaussThread(const std::vector<LCCobos_t>
                &sameRowLCs, const std::vector<LCCobos_t> &sameColLCs,
                const std::vector<idx_t> sameRowIdxs, const std::vector<idx_t>
                sameColIdxs, bool compoIsTrivial, const LCCobos_t &isoInv,
                coeff_t y, const vecTangles_t &domains, const vecTangles_t
                &codomains, const tangle_t &domain, const tangle_t &codomain,
                int threadId, int numThreads, std::vector<LCCobos_t*> &result,
                std::mutex &resultMtx, std::condition_variable &compFinishedCv);

        /** To be called after a tangle in a complex was delooped. Makes copies
         * copies of the relevant row (i.e. the row with index idx) and appends
         * them to the bottom of the matrix. Then modifies the entries of the
         * relevant and the copies rows by calling modifyDeloopCopy.
         * @pre this is a matrix from lowerTangles to upperTangles, the tangle
         * at index idx in upperTangles has been delooped, i.e. copies summands
         * were appended to the bottom of upperTangles (and the matrix size does
         * not agree with the codomain size for now).
         * @param idx The index of the tangle in the codomain that was delooped.
         * @param copies The number of copies that were made (careful: the
         * tangle at index idx is reused, i.e. for the original Khovanov
         * homology copies is 1, since a tangle with a circle is isomorphic to
         * the sum of two tangles.
         * @param lowerTangles The domain of the matrix
         * @param upperTangles The codomain of the matrix
         * @post The matrix represents now the correct map from lowerTangles to
         * upperTangles.
         */
        void deloop(word64 idx, int copies, const tangleCont_t &lowerTangles,
                const tangleCont_t &upperTangles, bool left);

        void deleteNonIsos();

        void writeToBin(std::ofstream &f) const;
        explicit MatLCCobos(std::ifstream &f, bool intCoefficients) :
            morphisms(f, intCoefficients) { }
#ifndef notests
        bool isSane(const Boundary *b = nullptr,
                const tangleCont_t *domain = nullptr,
                const tangleCont_t *codomain = nullptr) const;
        bool multIsNonZero(const MatLCCobos<cobordism_tpl> &other,
                const tangleCont_t &lower,
                const tangleCont_t &middle, const tangleCont_t &upper) const;
#endif
        std::ostream& detailedOutput(std::ostream &os) const {
            return morphisms.detailedOutput(os);
        }
    private:
        SparseMat<LCCobos<cobordism_tpl> > morphisms;
};


/** \class LCCobos
 * Don't inherit from this class.
 */
template <class cobordism_tpl>
class LCCobos {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        typedef typename cobordism_tpl::coeff_t coeff_t;
        typedef std::vector<cobordism_tpl> cobosCont_t;
        typedef typename cobordism_tpl::tangle_t tangle_t;
        //gets Identity
        explicit LCCobos(const typename cobordism_tpl::tangle_t &idTangle) {
            cobordisms.emplace_back(idTangle);
            assert(isSane());
        }
        
        LCCobos() {};
        /** Composes with other LCCobo (which way?)
         * @pre other is a LCCobos whose domain is this' codomain (or the other
         * way round?)
         * @post is sortfactored.
         */
        void compose(const LCCobos<cobordism_tpl>&other, const tangle_t &lower,
                const tangle_t &middle, const tangle_t &upper);
        /** ...
         * @post is sortfactored.
         */
        void setToUnion(const typename cobordism_tpl::tangle_t& lowerLeft,
                const typename cobordism_tpl::tangle_t& lowerRight,
                const typename cobordism_tpl::tangle_t& upperLeft,
                const typename cobordism_tpl::tangle_t& upperRight,
                const LCCobos<cobordism_tpl>&xLeft,
                const LCCobos<cobordism_tpl>&xRight);
        /** ...
         * @post is sortfactored.
         */
        void glue(const typename cobordism_tpl::tangle_t& lower,
                const typename cobordism_tpl::tangle_t& upper,
                const boundary_t gluePoint[2], boundary_t boundarySize);
        /**
         * Appends a cobordism to the linear combination of this.
         * @pre this->cobordisms must be sorted.
         * @param x The cobordism to be appended, must have same domain and
         * codomain as this and have non-zero coefficient.
         * @post this->cobordisms is sortFactored.
         */
        void add(cobordism_tpl x);


        bool isZero() const { return cobordisms.empty(); }
        /** Returns whether the linear combination is invertible, i.e. has only
         * one summand with invertible coefficient who is an invertible
         * cobordism.
         * @pre Cobordisms must be sorted and factored, otherwise might return
         * false for invertible this.
        */
        bool obviouslyNotInvertible() const { return (cobordisms.size() != 1); }
        bool isInvertible(const tangle_t &lowerTangle,
                const tangle_t &upperTangle) const;
        /** @pre Is invertible.
         */
        bool compositionIsTrivial(const tangle_t &/*lowerTangle*/,
                const tangle_t &/*upperTangle*/) const { return true; }

        /**
         * Adds the other a cobordism to the linear combination of this.
         * @pre this->cobordisms must be sorted.
         * @param x The LCCobos to be appended, must have same domain and
         * codomain as this and have non-zero coefficient.
         * @post this->cobordisms is sortFactored.
         */
        void add(LCCobos<cobordism_tpl> other);

        /** Sets this to the negative inverse of other. The original this is
         * lost.  Undefined behaviour if this is not invertible.
         * @pre cobo is empty (temporary).
         */
        void setToNegInv(LCCobos<cobordism_tpl>);
        void switchSign() {
            assert(isSane());
            for (typename std::vector<cobordism_tpl>::iterator i =
                    cobordisms.begin(); i != cobordisms.end(); ++i)
                i->switchSign();
            assert(isSane());
        }
	bool reducify(int root);
        void getCoeff(coeff_t& y) const { y = cobordisms.front().coefficient; }
        const coeff_t* getCoeff() const {
            return &(cobordisms.front().coefficient);
        }
        void operator*=(coeff_t y) {
            for (typename std::vector<cobordism_tpl>::iterator i =
                    cobordisms.begin(); i != cobordisms.end(); ++i)
                i->coefficient *= y;
        }
#ifndef notests
        bool isSane(const Boundary *b = nullptr,
                const tangle_t *domain = nullptr,
                const tangle_t *codomain = nullptr) const;
#endif
        /** Is called by deloopLeft and deloopRight, changes all cobordisms in
         * the linear combination accordingly.
         * @pre This is an entry of a matrix with a tangle that was delooped in
         * the codomain (if left) or domain (if !left); more precisely, this is
         * in the original row (col if !left) if kind == 0, and in the kind-th
         * copy of the original row (col if !left). ...
         * @post Cobordisms now go from lower to upper. Cobordisms are
         * sortfactored.
         * @param kind is 0 for the tangle, 1 for the first copy etc.
         * @param left whether the delooped tangle is in the codomain of the
         * matrix that this is part of.
         * @param lower codomain of the matrix that this is part of.
         * @param upper domain of the matrix that this is part of.
         */
        void modifyDeloopCopy(int kind, bool left,
                const typename cobordism_tpl::tangle_t &lower,
                const typename cobordism_tpl::tangle_t &upper);

        void writeToBin(std::ofstream &f) const;
        explicit LCCobos(std::ifstream &f, bool intCoefficients);
        std::ostream& detailedOutput(std::ostream &os) const;
    private:
        std::vector<cobordism_tpl> cobordisms;
        /** Sorts cobordisms, and replaces a_1x + a_2x by (a_1 + a_2)x
         * (="factor"). May take some time (O(n*log(n)) for n = number of
         * summands), and should only be called if the cobordisms are in random
         * order (e.g. after an operation as modifyDeloopCopy which changes
         * every cobordism).
         * @pre None.
         * @post Cobordisms are sortfactored.
         */
        void sortAndFactor();
        void factor();
};
