/*
 *
 *    src/planar_algebra/coefficient_rings_explicitTemplates.cpp --- Part of khoca, a knot homology calculator
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

#define INSTANTIATE_ALL(COEFF) \
template bool Monomial<COEFF>::isNonZero() const; \
template void Polynomial<COEFF>::inv();\
template void Polynomial<COEFF>::switchSign();\
template void Polynomial<COEFF>::operator+=(const Polynomial<COEFF>& r);\
template void Polynomial<COEFF>::operator*=(const Polynomial<COEFF>& r);\
template Polynomial<COEFF> Polynomial<COEFF>::operator*(const Monomial<COEFF>& r) const;\
template void Polynomial<COEFF>::operator*=(int r);\
template bool Polynomial<COEFF>::isInvertible() const;\
template bool Polynomial<COEFF>::isNonZero() const;

INSTANTIATE_ALL(MInteger)
INSTANTIATE_ALL(MRational)
INSTANTIATE_ALL(FF<uint8_t>)
INSTANTIATE_ALL(FF<uint16_t>)

#undef INSTANTIATE_ALL
