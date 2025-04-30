/*
 *
 *    src/shared.h --- Part of khoca, a knot homology calculator
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

/* For a release build the following switches should be defined: */
#define NDEBUG
#define notests
#define getsize
#define USEZIPDOTS
#define PRINTCOMPLEX

/* The following should only be defined for debugging. */
//#define USEOLDDOTS

#ifndef USEZIPDOTS
#ifndef USEOLDDOTS
At least one of USEZIPDOTS and USEOLDDOTS must be defined.
#endif
#endif

#if defined(_MSC_VER) || defined(__BORLANDC__)
    typedef __int64 word64;
    #define W64LIT(x) (word64) x
#else
    typedef long long word64;
    #define W64LIT(x) (word64) x
#endif


int insane();

class io {
    public:
        static std::ostream *cprogress_s;
        static std::ostream& cprogress() { return *cprogress_s; }
        static std::string markUp;
        static std::string noMarkUp;
};

/** \typedef Type to be used for a number that isn't greater than the number
  * of boundary points of a tangle.
  */
typedef int8_t boundary_t;
const boundary_t boundary_t_max = INT8_MAX;
typedef int16_t qShift_t;

/** Decreases by 1 all entries in v that are greater or equal to x.
 */
template<class E_tpl>
void shiftWhatIsHigher(std::vector<E_tpl> &v, const E_tpl &x);

char signedCharToAlph(signed char x);

template<class E_tpl> bool alphToChar(char c, E_tpl &x);

bool alphToSignedChar(char c, signed char &x);

/** Writes an arbitrary class to a binary file.
 */
template <class E> void writeToBinTpl(std::ofstream &f, const E& data) {
    f.write((char*)&data, sizeof(E));
}

/** Reads an arbitrary class from a binary file.
 */
template <class E> void readFromBinTpl(std::ifstream &f, E& data) {
    f.read((char*)&data, sizeof(E));
}

#ifndef MULTICALLER
#define MULTICALLER \
 INSTANTIATE_META(8) \
 INSTANTIATE_META(16) \
 INSTANTIATE_META(24) \
 INSTANTIATE_META(32) \
 INSTANTIATE_META(48) \
 INSTANTIATE_META(64) \
 INSTANTIATE_META(80) \
 INSTANTIATE_META(96) \
 INSTANTIATE_META(112) \
 INSTANTIATE_META(128) \
 INSTANTIATE_META(160) \
 INSTANTIATE_META(192) \
 INSTANTIATE_META(224) \
 INSTANTIATE_META(256)
#endif
