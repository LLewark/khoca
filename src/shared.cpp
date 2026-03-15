/*
 *
 *    src/shared.cpp --- Part of khoca, a knot homology calculator
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

#include <vector>
#include <inttypes.h>
#include <string.h>
#include <fstream>
#include <iostream>

#include "shared.h"

template void shiftWhatIsHigher<boundary_t>(std::vector<boundary_t> &v,
        const boundary_t &x);
template void shiftWhatIsHigher<word64>(std::vector<word64> &v,
        const word64 &x);

template bool alphToChar<uint8_t>(char c, uint8_t &x);
template bool alphToChar<uint16_t>(char c, uint16_t &x);

int insane() {
    return 0;
}

template<class E_tpl>
void shiftWhatIsHigher(std::vector<E_tpl> &v, const E_tpl &x) {
    for (auto i = v.begin(); i != v.end(); ++i)
        if (*i >= x)
            *i -= 1;
}

char signedCharToAlph(signed char x) {
    if (x < 0)
        return (char)(65 - x);
    else
        return (char)(97 + x);
}

template<class E_tpl> bool alphToChar(char c, E_tpl &x) {
    if ((c < 97) || (c > 122))
        return false;
    x = c - 97;
    return true;
}

bool alphToSignedChar(char c, signed char &x) {
    if ((c > 96) && (c < 123))
        x = c - 97;
    else if ((c > 64) && (c < 91))
        x = 65 - (int)c;
    else
        return false;
    return true;
}
