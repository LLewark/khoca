/*
 *
 *    src/krasner/krasner.h --- Part of khoca, a knot homology calculator
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

template <class coeff_tpl, int bitSize> class KrasnerCobo;
template <int bitSize> class KrasnerCoboData;

/** \class KrasnerTangleData
  Encapsulates the tangle. This is a separate class so one can more
  easily change the data structured (if e.g. one finds a more compressed
  structure).
 */
class KrasnerTangleData {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        boundary_t size() const { return pairing.size(); }
        boundary_t at(boundary_t idx) const { return pairing.at(idx); }
        void set(boundary_t idx, boundary_t val) {
            pairing.at(idx) = val; }
        void push_back(boundary_t val) { pairing.push_back(val); }
        void append(const KrasnerTangleData &other) { assert(this != &other);
            pairing.insert(pairing.end(), other.pairing.begin(),
                    other.pairing.end()); }
        void erase(boundary_t idx) { pairing.erase(pairing.begin() + idx,
                pairing.begin() + idx + 1); }
        boundary_t getCircleCount() const { return circleCount; }
        void setCircleCount(boundary_t val) { circleCount = val; }
        void incCircleCount() { circleCount += 1; }
        void decCircleCount() { circleCount -= 1; }
        void shiftWhatIsHigherK(boundary_t x) { shiftWhatIsHigher(pairing, x); }

        bool operator==(const KrasnerTangleData &other) const {
            return ((circleCount == other.circleCount) &&
                    (pairing == other.pairing)); }

    private:
        /* This should then be more compressed! */
        boundary_t circleCount;
        std::vector<boundary_t> pairing;
};

class KrasnerTangle : public Tangle<KrasnerTangle> {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        KrasnerTangle() { }

        static int16_t N;

        boundary_t getCircleCount() const { return data.getCircleCount(); }

        void naiveEuler() const;
        void print() const;
        boundary_t boundarySize() const { return data.size(); }
        /** Boundary points are ordered. This induces an order on
         * the components of the tangle:
         * First come intervals, then come circles.
         * The intervals are sorted by the smaller of the numbers
         * of their two endpoints.
         *
         * This, in turn, induces an order on the facets of a cobordism
         * between two webs:
         * First come facets which intersects the tangle both above
         * and below, then come facets whose boundary is a circle in
         * the domain web, and finally those whose boundary is a circle
         * in the codomain web.
         * Facets which intersect both above and below are sorted by
         * the smallest number of an interval they intersect in the
         * domain web.
         */
        void glue(const boundary_t gluePoints[2]);
        void setToUnion(const KrasnerTangle &o1,
                const KrasnerTangle &o2);

        bool operator==(const KrasnerTangle &other) const {
            return (data == other.data); }
        bool operator!=(const KrasnerTangle &other) const {
            return !(*this == other); }

        //vidx_t deloopNumberOfCopies() const;
        void deloop(std::vector<KrasnerTangle> &newCopies);
        bool hasLoop() const;
        boundary_t connect(boundary_t i) const { return data.at(i); }

#ifndef notests
        bool isSane(const Boundary *b = nullptr) const;
#endif
        void writeToBin(std::ofstream &f) const;
        explicit KrasnerTangle(std::ifstream &f, boundary_t size);

        std::ostream& detailedOutput(std::ostream &os) const {
            os << *this;
            return os;
        }
        friend std::ostream& operator<<(std::ostream& os,
                const KrasnerTangle& x) {
            os << "(q^" << (int)(x.qShift) << ":";
            for (int i = 0; i < (int)(x.data.size()); ++i)
                os << (int)(x.data.at(i)) << ((i + 1 < (int)(x.data.size())) ? " " : "");
            if (x.data.getCircleCount() != 0)
                os << " [c=" << (int)(x.data.getCircleCount()) << "]";
            os << ")";
            return os;
        }
    //private:
        KrasnerTangleData data;
        void pairUp(boundary_t x, boundary_t y);
        void deleteAndShiftIdx(boundary_t x);

    //friend KrasnerCoboData;
};

/** \class see class KrasnerTangleData
 */
template <int bitSize> class KrasnerCoboData {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        /** standard vector operations */
        boundary_t dotsSize() const {
            assert(compressedIsCorrect());
#ifdef USEZIPDOTS
            return nbFacets;
#else
            return dots.size();
#endif
        }
        /** Won't be needed anymore */
        void resize(int s) {
            assert(compressedIsCorrect());
#ifdef USEOLDDOTS
            dots.resize(s, 0);
#endif
#ifdef USEZIPDOTS
            nbFacets = s;
#endif
            assert(compressedIsCorrect());
        }
        int dotsAt(boundary_t idx) const {
            assert(compressedIsCorrect());
#ifdef USEZIPDOTS
            assert(((int)idx + 1) * bitsPerDot <= bitSize);
            int result = ((compressedDots << (int)idx * bitsPerDot) >>
                    (bitSize - bitsPerDot)).to_ulong();
#ifdef USEOLDDOTS
            assert(result == dots.at(idx));
#endif
#else
            int result = dots.at(idx);
#endif
            return result;
        }
        void set(boundary_t i, int x, bool checkSanity = true) {
            if (checkSanity) assert(compressedIsCorrect());
#ifdef USEZIPDOTS
            assert(((int)i + 1) * bitsPerDot <= bitSize);
            assert((x >> bitsPerDot) == 0);
            // sets the i-th slot to zeroes
            compressedDots &= ~(std::bitset<bitSize>((1 << bitsPerDot) - 1)
                    << (bitSize - ((int)i + 1) * bitsPerDot));
            compressedDots |= (std::bitset<bitSize>(x) <<
                    (bitSize - bitsPerDot)) >> ((int)i * bitsPerDot);
#endif
#ifdef USEOLDDOTS
            dots.at(i) = x;
#endif
            if (checkSanity) assert(compressedIsCorrect());
        }
        void insert(boundary_t i, int x) {
            assert(compressedIsCorrect());
#ifdef USEZIPDOTS
            assert((nbFacets + 1) * bitsPerDot <= bitSize);
            compressedDots = ((compressedDots >>
                        (bitSize - (int)i * bitsPerDot))
                    << (bitSize - (int)i * bitsPerDot))
                | ((std::bitset<bitSize>(x) << (bitSize - bitsPerDot))
                    >> ((int)i * bitsPerDot))
                | ((compressedDots << ((int)i * bitsPerDot))
                    >> (((int)i + 1) * bitsPerDot));
            nbFacets++;
#endif
#ifdef USEOLDDOTS
            dots.insert(dots.begin() + i, x);
#endif
            assert(compressedIsCorrect());
        }
        void erase(boundary_t i) {
            assert(compressedIsCorrect());
#ifdef USEOLDDOTS
            dots.erase(dots.begin() + i);
#endif
#ifdef USEZIPDOTS
            compressedDots = ((compressedDots >>
                        (bitSize - (int)i * bitsPerDot))
                    << (bitSize - (int)i * bitsPerDot))
                | ((compressedDots << (((int)i + 1) * bitsPerDot))
                     >> ((int)i * bitsPerDot));
            nbFacets--;
#endif
            assert(compressedIsCorrect());
        }
        bool operator<(const KrasnerCoboData<bitSize> &other) const {
#ifdef USEZIPDOTS
            assert(other.nbFacets == nbFacets);
            bool result = false;
            for (int i = bitSize - 1; i >= bitSize - nbFacets * bitsPerDot; i--)
                if (compressedDots[i] ^ other.compressedDots[i]) {
                    result = other.compressedDots[i];
                    break;
                }
#ifdef USEOLDDOTS
            assert(result == (dots < other.dots));
#endif
#else
            bool result = (dots < other.dots);
#endif
            return result;
        }
        bool operator==(const KrasnerCoboData<bitSize> &other) const {
#ifdef USEZIPDOTS
            assert(other.nbFacets == nbFacets);
            const int shift = bitSize - nbFacets * bitsPerDot;
            bool result = (other.compressedDots >> shift)
                == (compressedDots >> shift);
#ifdef USEOLDDOTS
            assert(result == (dots == other.dots));
#endif
#else
            bool result = (dots == other.dots);
#endif
            return result;
        }

        /** Rearrangements, used for tensor. */
        void mashTogether(const KrasnerTangle& lowerLeft,
                          const KrasnerTangle& lowerRight,
                          const KrasnerTangle& upperLeft,
                          const KrasnerTangle& upperRight,
                          const KrasnerCoboData<bitSize> &coboLeft,
                          const KrasnerCoboData<bitSize> &coboRight);
        /** Rearrangements, used for composition. */
        void combine(boundary_t facetCountNew, boundary_t lowerCircleCount,
            boundary_t middleCircleCount, boundary_t upperCircleCount,
            const KrasnerCoboData<bitSize> &first,
            const KrasnerCoboData<bitSize> &other);

        void readFromBin(std::ifstream &f, boundary_t nbFacets_) {
#ifdef USEZIPDOTS
            nbFacets = nbFacets_;
#endif
#ifdef USEOLDDOTS
            dots.resize(nbFacets_);
#endif
            for (boundary_t i = 0; i < nbFacets_; ++i) {
                uint16_t newVal;
                readFromBinTpl(f, newVal);
                set(i, newVal, false);
            }
            assert(compressedIsCorrect());
        }
        static void setBitsPerDot(int
#ifdef USEZIPDOTS
                bitsPerDot_
#endif
                ) {
#ifdef USEZIPDOTS
            bitsPerDot = bitsPerDot_;
#endif
        }
#ifdef USEZIPDOTS
        KrasnerCoboData() : nbFacets(0) { }
#else
        KrasnerCoboData() { }
#endif
    private:
        /** @pre: sanity, dots-vector is taken care of.
         */
#ifdef USEZIPDOTS
        void insertZeroes(boundary_t idx, boundary_t count) {
            compressedDots = ((compressedDots >> (bitSize - idx * bitsPerDot))
                    << (bitSize - idx * bitsPerDot))
                | ((compressedDots << (idx * bitsPerDot))
                    >> ((idx + count) * bitsPerDot));
            nbFacets += count;
            assert(bitSize >= nbFacets);
        }
        /** @pre: sanity; dots-vector must be taken care of.
         */
        void insertRange(boundary_t idx, const KrasnerCoboData<bitSize> &other,
                boundary_t from, boundary_t to) {
            compressedDots = ((compressedDots >> (bitSize - idx * bitsPerDot))
                    << (bitSize - idx * bitsPerDot))
                | (((other.compressedDots << (from * bitsPerDot))
                    >> (bitSize + from * bitsPerDot - to * bitsPerDot))
                    << (bitSize - to * bitsPerDot + from * bitsPerDot
                        - idx * bitsPerDot))
                | ((compressedDots << (idx * bitsPerDot))
                    >> ((idx + to - from) * bitsPerDot));
            nbFacets = nbFacets + to - from;
            assert(bitSize >= nbFacets);
        }
#endif
        bool compressedIsCorrect() const {
#ifdef USEZIPDOTS
#ifdef USEOLDDOTS
            if ((int)dots.size() != nbFacets) return false;
            for (int i = 0; i < nbFacets; ++i) {
                int result = ((compressedDots << i * bitsPerDot) >>
                        (bitSize - bitsPerDot)).to_ulong();
                if (result != dots.at(i)) return false;
            }
#endif
#endif
            return true;
        }
        static int bitsPerDot;

#ifdef USEOLDDOTS
        std::vector<int16_t> dots;
#endif
#ifdef USEZIPDOTS
        boundary_t nbFacets;
        /** The first entry is at the left, i.e. at the most significant bit,
         * which has the highest index (for the []-operator of std::bitset).*/
        std::bitset<bitSize> compressedDots;
#endif
};

template <class coeff_tpl, int bitSize>
class KrasnerCobo : 
public Cobo<KrasnerTangle, KrasnerCobo<coeff_tpl, bitSize>, coeff_tpl> {
    public:
        /** The index of sl_N. */
        typedef int16_t N_t;
        static std::vector<coeff_tpl> frobenius;
        static N_t N() { return frobenius.size() - 1; }
        static void guaranteeMultVector(int x) {
            for (int i = 0; i <= x - (int)N() + 1 - (int)multVector.size();
                    ++i) {
                if (multVector.empty()) {
                    multVector.resize(1);
                    for (auto j = frobenius.begin();
                            j != frobenius.end() - 1; ++j) {
                        multVector.back().push_back(*j);
                        multVector.back().back().switchSign();
                    }
                } else {
                    multVector.push_back(multVector.back());
                    multVector.back().insert(multVector.back().begin(),
                            coeff_tpl(0));
                    for (auto j = frobenius.begin();
                            j != frobenius.end() - 1; ++j) {
                        coeff_tpl toSubstract = multVector.back().back();
                        toSubstract *= *j;
                        toSubstract.switchSign();
                        multVector.back().at(j - frobenius.begin()) +=
                            toSubstract;
                    }
                    multVector.back().pop_back();
                }
            }
        }
        static void initialiseFrobenius(const std::vector<int> &F, int N) {
            frobenius.clear();
            multVector.clear();
            coeff_tpl::initialiseFrobenius(frobenius, F, N);
        }
        // multVector contains the value of X^n for all (needed) n.
        static std::vector<std::vector<coeff_tpl> > multVector;

        static void frobX(std::vector<coeff_tpl> &x);
        static void frobXn(std::vector<coeff_tpl> &x, int n);
        static void frobMult(std::vector<coeff_tpl> &x, const std::vector<coeff_tpl> &y);
        static void frobGenus(std::vector<coeff_tpl> &x);
        static void frobComult(std::vector<coeff_tpl> &x, int d);
        static void composeHelper(boundary_t facetCount,
                const KrasnerTangle &middle, std::vector<int> &groups,
                const std::vector<boundary_t> &belongsToNew,
                std::vector<std::vector<coeff_tpl> > &newDots,
                const std::vector<boundary_t> &belongsToLower,
                const KrasnerCoboData<bitSize>& data,
                std::vector<int> &genusCorrection);

        using Cobo<KrasnerTangle, KrasnerCobo<coeff_tpl, bitSize>,
              coeff_tpl>::coefficient;
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        typedef coeff_tpl coeff_t;
        typedef KrasnerTangle tangle_t;
        typedef KrasnerCobo<coeff_tpl, bitSize> this_t;

        KrasnerCobo() { }
        /** Constructs this as identity cobordism of idTangle.
         * @pre The idTangle doesn't have any circles (it would be possible to
         *      allow it, but it isn't necessary).
         */
        explicit KrasnerCobo(const KrasnerTangle &idTangle) :
                Cobo<KrasnerTangle, this_t, coeff_tpl>(1) {
            assert(idTangle.getCircleCount() == 0);
            data.resize(idTangle.boundarySize() / 2);
        }

        /** This is correct.
        */
        qShift_t degree(boundary_t boundarySize) const;
        /** Careful: The result has to be corrected by
         *  - (# boundary points) / 2.
        */
        qShift_t relativeDegree() const;
        void compose(const KrasnerCobo<coeff_tpl, bitSize> &other,
                std::vector<this_t> &result, const tangle_t &lower,
                const tangle_t &middle, const tangle_t &upper) const;
        void setToUnion(const KrasnerTangle&, const KrasnerTangle&,
                const KrasnerTangle&, const KrasnerTangle&,
                const this_t &o1, const this_t &o2);

        void glue(const KrasnerTangle& lower, const KrasnerTangle& upper,
                const boundary_t gluePoint[2], boundary_t boundarySize,
                LCCobos<this_t> &father);
        void print() const;
        /** Format: 
         * Numbers are encoded as usual, then one number per dot - that's it.
         */
        /*bool readFromString(std::string::const_iterator &i,
                std::string::const_iterator end);*/
        void setToRandom() {}

        void modifyDeloopCopy(int kind, bool left, std::vector<this_t> &v,
                const KrasnerTangle&/* lower*/, const KrasnerTangle& upper);

        bool operator<(const this_t &other) const;
        bool operator==(const this_t &other) const;
        bool operator!=(const this_t &other) const {
            return !(*this == other);
        }

	int reducify();

        bool isInvertible(const tangle_t &lowerTangle,
                const tangle_t &upperTangle) const;
        bool isEmpty() const;
#ifndef notests
        bool isSane(const Boundary *b = nullptr,
                const tangle_t *domain = nullptr,
                const tangle_t *codomain = nullptr) const;
#endif

        void writeToBin(std::ofstream &f) const;
        explicit KrasnerCobo(std::ifstream &f, bool intCoefficients);
        std::ostream& detailedOutput(std::ostream &os) const {
            os << *this;
            return os;
        }
        friend std::ostream& operator<<(std::ostream& os,
                const KrasnerCobo<coeff_tpl, bitSize>& x) {
            os << "(" << x.coefficient << ")";
            if (x.data.dotsSize() && (int)x.data.dotsAt(0)) {
                os << "*X";
                if ((int)x.data.dotsAt(0) > 1)
                    os << "^" << (int)x.data.dotsAt(0);
            }
            return os;
        }
    private:
        KrasnerCoboData<bitSize> data;
};

template <class coeff_tpl, int bitSize>
boundary_t countFacets(std::vector<boundary_t> &belongsTo,
        boundary_t boundarySize,
        const typename KrasnerCobo<coeff_tpl, bitSize>::tangle_t &domain,
        const typename KrasnerCobo<coeff_tpl, bitSize>::tangle_t &codomain);
