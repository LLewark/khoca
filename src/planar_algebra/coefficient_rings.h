/*
 *
 *    src/planar_algebra/coefficient_rings.h --- Part of khoca, a knot homology calculator
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

/** \class Monomial
 * Note that exponents may have varying sizes, and that the coefficient may be
 * zero.
 */
template <typename coeff_tpl>
class Monomial {
    public:
        typedef Monomial<coeff_tpl> Monomial_t;

        bool isSane(const Monomial_t *compare) const;

        void inv();
        void switchSign();
        void operator*=(const Monomial_t&);
        bool operator<(const Monomial_t&) const;
        void operator*=(int r);
        bool isInvertible() const;
        bool isNonZero() const;
        Monomial<coeff_tpl> operator*(const Monomial_t&) const;

        Monomial<coeff_tpl>() { throw; }
        Monomial<coeff_tpl>(int x) : coefficient(x) { assert(x != 0); }
        Monomial<coeff_tpl>(const Monomial_t &x, const Monomial_t &y)
            : exponents(x.exponents), coefficient(x.coefficient) {
                coefficient += y.coefficient;
            }
        Monomial<coeff_tpl>(const std::vector<int16_t>
                exponents, const coeff_tpl &coefficient_) :
            exponents(exponents), coefficient(coefficient_) { }
        friend std::ostream& operator<<(std::ostream &os,
                const Monomial<coeff_tpl> &x) {
            bool isConstant = true;
            for (auto i = x.exponents.begin(); i != x.exponents.end(); ++i)
                if (*i) {
                    if (isConstant) {
                        x.coefficient.writeAsCoefficient(os);
                        isConstant = false;
                    }
                    os << signedCharToAlph(i - x.exponents.begin());
                    if (*i != 1)
                        os << "^" << *i;
                }
            if (isConstant)
                os << x.coefficient;
            return os;
        }
    private:
        std::vector<int16_t> exponents;
        coeff_tpl coefficient;
};

/** \class Monomial
  * Contains no zero monomials.
 */
template <typename coeff_tpl>
class Polynomial {
    public:
        typedef Polynomial<coeff_tpl> Polynomial_t;
        typedef Monomial<coeff_tpl> Monomial_t;
        typedef coeff_tpl Coefficient_t;

        bool isSane() const;

        static void printRing(int N, std::ostream& s) {
            coeff_tpl::printRing(N, s);
            s << "[";
            for (int i = 0; i < N; ++i) {
                if (i)
                    s << ", ";
                s << signedCharToAlph(i);
            }
            s << "]";
        }
        static void initialiseFrobenius(std::vector<Polynomial_t> &frobenius, const std::vector<int> &, int N) {
            for (int i = 0; i < N; ++i) {
                std::vector<int16_t> exponents(N, 0);
                exponents.at(i) = 1;
                frobenius.push_back(Polynomial_t(Monomial_t(exponents, coeff_tpl(1))));
            }
            frobenius.push_back(Polynomial_t(Monomial_t(std::vector<int16_t>(N, 0), coeff_tpl(1))));
        }

        const mpz_t* getmpz_t() const { throw; }

        Polynomial<coeff_tpl>() { }
        explicit Polynomial<coeff_tpl>(const Monomial_t& r) {
            if (r.isNonZero())
                monoms.push_back(r);
        }
        Polynomial<coeff_tpl>(int x) {
           if (x)
               monoms.emplace_back(x);
        }

        static uint16_t coefficientTypeToUint() {
            throw;
            return coeff_tpl::coefficientTypeToUint();
        }
        void writeToBin(std::ofstream &) const { throw; }
        explicit Polynomial<coeff_tpl>(std::ifstream &) { throw; }

        void inv();
        void switchSign();
        void operator+=(const Polynomial_t& r);
        void operator*=(const Polynomial_t& r);
        Polynomial<coeff_tpl> operator*(const Monomial_t& r) const;
        void operator*=(int r);
        bool isInvertible() const;
        bool isNonZero() const;

        friend std::ostream& operator<<(std::ostream &os,
                const Polynomial<coeff_tpl> &x) {
            if (x.monoms.empty())
                os << "0";
            else
                for (auto i = x.monoms.begin(); i != x.monoms.end(); ++i) {
                    if (i != x.monoms.begin())
                        os << " + ";
                    os << *i;
                }
            return os;
        }

    private:
        std::vector<Monomial_t> monoms;
};

/** \class MRational
 * Encapsulates the mpq_t type of the MPIR library
 * (arbitrary precision rational number).
 */
class MRational {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        static uint16_t coefficientTypeToUint() { return 1; }
        static void printRing(int, std::ostream& s) {
            s << "Q";
        }
        static void initialiseFrobenius(std::vector<MRational> &frobenius, const std::vector<int> &F, int) {
            for (auto i = F.begin(); i != F.end(); ++i)
                frobenius.push_back(MRational(*i));
            frobenius.push_back(MRational(1));
        }

        const mpz_t* getmpz_t() const { throw; }

        explicit MRational(const int &other);
        MRational() : MRational(1) { }
        explicit MRational(std::ifstream &);

        ~MRational();
        MRational(const MRational &other);
        void operator=(const MRational &other);

        MRational(MRational &&other);
        void operator=(MRational &&other);

        void inv();
        void switchSign();
        void operator+=(const MRational& r);
        void operator*=(const MRational& r);
        void operator*=(int r);
        bool isInvertible() const;
        bool isNonZero() const;

        void writeToBin(std::ofstream &) const;

        std::ostream& writeAsCoefficient(std::ostream &os) const;
        friend std::ostream& operator<<(std::ostream &, const MRational &);
    private:
        mpq_t val;
};

class MInteger {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const;
#endif
        static void printRing(int, std::ostream& s) {
            s << "Z";
        }
        static uint16_t coefficientTypeToUint() { return 0; }
        static void initialiseFrobenius(std::vector<MInteger> &frobenius, const std::vector<int> &F, int) {
            for (auto i = F.begin(); i != F.end(); ++i)
                frobenius.push_back(MInteger(*i));
            frobenius.push_back(MInteger(1));
        }

        explicit MInteger(const int &other);
        MInteger() : MInteger(1) { }
        explicit MInteger(std::ifstream &);

        ~MInteger();
        MInteger(const MInteger &other);
        void operator=(const MInteger &other);

        MInteger(MInteger &&other);
        void operator=(MInteger &&other);

        void inv();
        void switchSign();
        void operator+=(const MInteger& r);
        void operator*=(const MInteger& r);
        void operator*=(int r);
        bool isInvertible() const;
        bool isNonZero() const;

        const mpz_t* getmpz_t() const { return &val; }

        void writeToBin(std::ofstream &) const;

        std::ostream& writeAsCoefficient(std::ostream &os) const;
        friend std::ostream& operator<<(std::ostream &, const MInteger &);
    private:
        mpz_t val;
};

uint64_t calcInverse(uint64_t i, uint64_t p);

template <typename val_tpl>
class FF {
    public:
#ifndef getsize
        void printSize(std::vector<word64> &s) const {
            s.at(6) += sizeof(*this);
        }
#endif
        static void printRing(int, std::ostream& s) {
            s << "F_" << (int)p;
        }
        static uint16_t coefficientTypeToUint() { return p; }
        static void initialiseFrobenius(std::vector<FF<val_tpl> > &frobenius, const std::vector<int> &F, int) {
            for (auto i = F.begin(); i != F.end(); ++i)
                frobenius.push_back(FF<val_tpl>(*i));
            frobenius.push_back(FF<val_tpl>(1));
        }

        const mpz_t* getmpz_t() const { throw; }

        FF<val_tpl>() : FF<val_tpl>(1) { }
        explicit FF<val_tpl>(const int &other) :
            val( ((other < 0) ? p : 0) + (other % p)) {
                assert(other < (int)(1 << 8*sizeof(val_tpl)));
            }
        explicit FF<val_tpl>(std::ifstream &f) { readFromBinTpl(f, val); }

        bool isNonZero() const { return val; }
        bool isInvertible() const { return isNonZero(); }

        static void setP(val_tpl p_) {
            p = p_;
            inverses.clear();
            inverses.reserve(p);
            inverses.push_back(0);
            for (val_tpl i = 1; i < p; ++i)
                inverses.push_back(calcInverse(i, p));
        }

        void operator+=(const FF<val_tpl> &x) {
            val = (((uint64_t)(val)) + ((uint64_t)(x.val))) % p;
        }
        void operator*=(const FF<val_tpl> &x) {
            val = (((uint64_t)(val)) * ((uint64_t)(x.val))) % p;
        }
        void operator*=(int r) { val = (p + (val * r) % p) % p; } 
        void switchSign() { val = val ? (p - val) : 0; }
        void writeToBin(std::ofstream &f) const { writeToBinTpl(f, val); }
        void inv() { assert(val); val = inverses.at(val); }

        std::ostream& writeAsCoefficient(std::ostream &os) const {
            if (val != 1)
                os << *this;
            return os;
        }
        friend std::ostream& operator<<(std::ostream& os,
                const FF<val_tpl>& x) { return os << (int)x.val; }
    private:
        static val_tpl p;
        static std::vector<val_tpl> inverses;
        val_tpl val;
};
