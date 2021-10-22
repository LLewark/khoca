/*
 *
 *    src/planar_algebra/coefficient_rings.cpp --- Part of khoca, a knot homology calculator
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

#include <array>
#include <stddef.h>
#include <gmp.h>
#include <cstdint>
#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>

#include "../shared.h"
#include "coefficient_rings.h"

#include "coefficient_rings_explicitTemplates.cpp"

template<> uint8_t FF<uint8_t>::p = 0;
template<> std::vector<uint8_t> FF<uint8_t>::inverses = std::vector<uint8_t>();
template<> uint16_t FF<uint16_t>::p = 0;
template<> std::vector<uint16_t> FF<uint16_t>::inverses =
        std::vector<uint16_t>();

#ifndef getsize
void MRational::printSize(std::vector<long long> &s) const {
    s.at(6) += sizeof(mp_limb_t) *
        (mpz_size(mpq_numref(val)) + mpz_size(mpq_denref(val)));
}

void MInteger::printSize(std::vector<long long> &s) const {
    s.at(6) += sizeof(mp_limb_t) * (mpz_size(val));
}
#endif

// Monomial

template <typename coeff_tpl>
void Monomial<coeff_tpl>::inv() {
    for (auto i = exponents.begin(); i != exponents.end(); ++i)
        *i = -*i;
    coefficient.inv();
}

template <typename coeff_tpl>
void Monomial<coeff_tpl>::switchSign() {
    coefficient.switchSign();
}

template <typename coeff_tpl>
void Monomial<coeff_tpl>::operator*=(const Monomial_t &other) {
    if (exponents.size() < other.exponents.size())
        exponents.resize(other.exponents.size(), 0);
    for (auto i = exponents.begin(); (i != exponents.end()) && (i - exponents.begin() < (int)other.exponents.size()); ++i)
        *i += other.exponents.at(i - exponents.begin());
    coefficient *= other.coefficient;
}

template <typename coeff_tpl>
Monomial<coeff_tpl> Monomial<coeff_tpl>::operator*(
        const Monomial_t &other) const {
    Monomial_t result = *this;
    result *= other;
    return result;
}

template <typename coeff_tpl>
void Monomial<coeff_tpl>::operator*=(int r) {
    coefficient *= r;
}

template <typename coeff_tpl>
bool Monomial<coeff_tpl>::isInvertible() const {
    return coefficient.isInvertible();
}

template <typename coeff_tpl>
bool Monomial<coeff_tpl>::isNonZero() const {
    return coefficient.isNonZero();
}

template <typename coeff_tpl>
bool Monomial<coeff_tpl>::operator<(const Monomial_t &other) const {
    int i = 0;
    for (; (i < (int)exponents.size()) && (i < (int)other.exponents.size()); ++i) {
        if (exponents.at(i) < other.exponents.at(i))
            return true;
        else if (exponents.at(i) > other.exponents.at(i))
            return false;
    }
    for (; i < (int)other.exponents.size(); ++i)
        if (other.exponents.at(i))
            return true;
    return false;
}

template <typename coeff_tpl>
bool Monomial<coeff_tpl>::isSane(const Monomial_t *compare) const {
    if (! coefficient.isNonZero())
        return false;
    if (compare) {
        bool isEqual = true;
        int i = 0;
        for (; (i < (int)exponents.size()) && (i < (int)compare->exponents.size()); ++i)
            if (exponents.at(i) != compare->exponents.at(i)) {
                isEqual = false;
                break;
            }
        if (isEqual) {
            for (; i < (int)exponents.size(); ++i)
                if (exponents.at(i))
                    isEqual = false;
            for (; i < (int)compare->exponents.size(); ++i)
                if (compare->exponents.at(i))
                    isEqual = false;
        }
        if (isEqual)
            return false;
    }
    return true;
}

// Polynomial

template <typename coeff_tpl>
bool Polynomial<coeff_tpl>::isSane() const {
    for (auto i = monoms.begin(); i != monoms.end(); ++i) {
        if (! i->isSane((i == monoms.begin()) ? nullptr : &*(i-1)))
            return false;
        if ((i != monoms.begin()) && (*i < *(i-1)))
            return false;
    }
    return true;
}

template <typename coeff_tpl>
void Polynomial<coeff_tpl>::inv() {
    assert(isInvertible());
    monoms.at(0).inv();
    assert(isSane());
}

template <typename coeff_tpl>
void Polynomial<coeff_tpl>::switchSign() {
    for (auto i = monoms.begin(); i != monoms.end(); ++i)
        i->switchSign();
    assert(isSane());
}

template <typename coeff_tpl>
void Polynomial<coeff_tpl>::operator+=(const Polynomial_t &other) {
    std::vector<Monomial_t> oldmonoms;
    monoms.swap(oldmonoms);
    auto i = oldmonoms.begin();
    auto j = other.monoms.begin();
    while ((i != oldmonoms.end()) || (j != other.monoms.end())) {
        if (j == other.monoms.end())
            monoms.push_back(*i++);
        else if (i == oldmonoms.end())
            monoms.push_back(*j++);
        else if (*i < *j)
            monoms.push_back(*i++);
        else if (*j < *i)
            monoms.push_back(*j++);
        else {
            monoms.emplace_back(*i++, *j++);
            if (! monoms.back().isNonZero())
                monoms.resize(monoms.size() - 1);
        }
    }
    assert(isSane());
}

template <typename coeff_tpl>
Polynomial<coeff_tpl> Polynomial<coeff_tpl>::operator*(
        const Monomial_t& r) const {
    if (r.isNonZero()) {
        Polynomial_t result = *this;
        for (auto i = result.monoms.begin(); i != result.monoms.end(); ++i)
            *i *= r;
        assert(result.isSane());
        return result;
    } else
        return Polynomial_t();
}

template <typename coeff_tpl>
void Polynomial<coeff_tpl>::operator*=(const Polynomial_t &other) {
    std::vector<Monomial_t> oldmonoms;
    monoms.swap(oldmonoms);
    for (auto i = oldmonoms.begin(); i != oldmonoms.end(); ++i)
        *this += other * (*i);
    assert(isSane());
}

template <typename coeff_tpl>
void Polynomial<coeff_tpl>::operator*=(int r) {
    if (r)
        for (auto i = monoms.begin(); i != monoms.end(); ++i)
            *i *= r;
    else
        monoms.resize(0);
    assert(isSane());
}

template <typename coeff_tpl>
bool Polynomial<coeff_tpl>::isInvertible() const {
    return (monoms.size() == 1) && (monoms.at(0).isInvertible());
}

template <typename coeff_tpl>
bool Polynomial<coeff_tpl>::isNonZero() const {
    return (monoms.size());
}

// Rest

uint64_t calcInverse(uint64_t i, uint64_t p) {
    assert(i);
    int64_t result = 0;
    int64_t prevResult = 1;
    while (i != 0) {
        uint64_t div = p / i;
        uint64_t mod = p % i;
        p = i;
        i = mod;

        uint64_t tmp = result;
        result = prevResult;
        prevResult = tmp - div * prevResult;
    }
    if (p != 1)
        throw;
    return ((result < 0) ? prevResult : 0) + result;
}

void saveMPIRint(std::ofstream &f, const mpz_t x, bool sign = true) {
    if (sign) {
	signed char s = (mpz_sgn(x) == -1) ? 1 : 0;
	writeToBinTpl(f, s);
    }
    const uint8_t size = mpz_sizeinbase(x, 256);
    writeToBinTpl(f, size);
    void* p = malloc(size);
    mpz_export(p, nullptr, 1, size, 1, 0, x);
    f.write((char *)p, size);
    free(p);
}

void loadMPIRint(std::ifstream &f, mpz_t x, bool sign = true) {
    bool isNeg = false;
    if (sign) {
	uint8_t s;
	readFromBinTpl(f, s);
	isNeg = s;
    }
    uint8_t size;
    readFromBinTpl(f, size);
    void* p = malloc(size);
    f.read((char *)p, size);
    mpz_import(x, size, 1, 1, 1, 0, p);
    if (isNeg)
	mpz_neg(x, x);
    free(p);
}



MRational::MRational(std::ifstream &f) {
    mpz_t num, denum;
    mpz_init(num);
    mpz_init(denum);
    loadMPIRint(f, num, true);
    loadMPIRint(f, denum, false);
    mpq_init(val);
    mpq_set_num(val, num);
    mpq_set_den(val, denum);
    mpz_clear(num);
    mpz_clear(denum);
}

void MRational::writeToBin(std::ofstream &f) const {
    saveMPIRint(f, mpq_numref(val), true);
    saveMPIRint(f, mpq_denref(val), false);
}

MRational::~MRational(){
    mpq_clear(val);
}

MRational::MRational(const int &other){
    mpq_init(val);
    mpq_set_si(val, other, 1);
}

MRational::MRational(const MRational &other){
    mpq_init(val);
    *this = other;
}

void MRational::operator=(const MRational &other){
    mpq_set(val, other.val);
}

MRational::MRational(MRational &&other){
    mpq_init(val);
    *this = other;
}

void MRational::operator=(MRational &&other){
    mpq_swap(val, other.val);
}

void MRational::inv(){
    mpq_inv(val, val);
}

void MRational::switchSign(){
    mpq_neg(val, val);
}

void MRational::operator+=(const MRational& r){
    mpq_add(val, val, r.val);
}

void MRational::operator*=(const MRational& r){
    mpq_mul(val, val, r.val);
}

void MRational::operator*=(int r){
    if (r) {
        unsigned long int g = mpz_gcd_ui(NULL, mpq_denref(val), std::abs(r));
        mpz_div_ui(mpq_denref(val), mpq_denref(val), g);
        mpz_mul_ui(mpq_numref(val), mpq_numref(val), std::abs(r) / g);
        if (r < 0)
            switchSign();
    } else
        mpq_set_si(val, 0, 1);
}

bool MRational::isInvertible() const{
    return isNonZero();
}

bool MRational::isNonZero() const{
    return (mpq_sgn(val) != 0);
}

std::ostream& operator<<(std::ostream &os, const MRational &x) {
    return (os << x.val);
}

std::ostream& MRational::writeAsCoefficient(std::ostream &os) const {
    if (mpq_cmp_si(val, 1, 1) != 0) {
        if (mpq_cmp_si(val, -1, 1) == 0)
            os << "-";
        else
            os << *this;
    }
    return os;
}

MInteger::MInteger(std::ifstream &f) {
    mpz_init(val);
    loadMPIRint(f, val, true);
}

void MInteger::writeToBin(std::ofstream &f) const {
    saveMPIRint(f, val, true);
}

MInteger::~MInteger(){
    mpz_clear(val);
}

MInteger::MInteger(const int &other){
    mpz_init(val);
    mpz_set_si(val, other);
}

MInteger::MInteger(const MInteger &other){
    mpz_init(val);
    *this = other;
}

void MInteger::operator=(const MInteger &other){
    mpz_set(val, other.val);
}

MInteger::MInteger(MInteger &&other){
    mpz_init(val);
    *this = other;
}

void MInteger::operator=(MInteger &&other){
    mpz_swap(val, other.val);
}

void MInteger::inv(){
    // Inverse integers (namely +-1) are self-inverse.
}

void MInteger::switchSign(){
    mpz_neg(val, val);
}

void MInteger::operator+=(const MInteger& r){
    mpz_add(val, val, r.val);
}

void MInteger::operator*=(const MInteger& r){
    mpz_mul(val, val, r.val);
}

void MInteger::operator*=(int r){
    mpz_mul_si(val, val, r);
}

bool MInteger::isInvertible() const{
    return (mpz_cmp_si(val, 1) == 0) || (mpz_cmp_si(val, -1) == 0);
}

bool MInteger::isNonZero() const{
    return (mpz_sgn(val) != 0);
}

std::ostream& MInteger::writeAsCoefficient(std::ostream &os) const {
    if (mpz_cmp_si(val, 1) != 0) {
        if (mpz_cmp_si(val, -1) == 0)
            os << "-";
        else
            os << *this;
    }
    return os;
}


std::ostream& operator<<(std::ostream &os, const MInteger &x) {
    return (os << x.val);
}
