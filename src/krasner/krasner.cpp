/*
 *
 *    src/krasner/krasner.cpp --- Part of khoca, a knot homology calculator
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
#include <iostream>
#include <vector>
#include <deque>
#include <bitset>
#include <inttypes.h>
#include <stdlib.h>
#include <string>
#include <sstream>
#include <algorithm>
// Next line is necessary because of mpir
#include <stddef.h>
#include <gmp.h>
#include <fstream>

#include "../shared.h"

#include <assert.h>

#include "../planar_algebra/coefficient_rings.h"
#include "../planar_algebra/sparsemat.h"
#include "../planar_algebra/planar_algebra.h"
#include "krasner.h"

#include "krasnerExplicitTemplates.cpp"

template <class coeff_tpl, int bitSize> std::vector<coeff_tpl>
    KrasnerCobo<coeff_tpl, bitSize>::frobenius = std::vector<coeff_tpl>();
template <class coeff_tpl, int bitSize> std::vector<std::vector<coeff_tpl> >
    KrasnerCobo<coeff_tpl, bitSize>::multVector =
    std::vector<std::vector<coeff_tpl> >();
int16_t KrasnerTangle::N = 0;

template<int bitSize> int KrasnerCoboData<bitSize>::bitsPerDot = -1;

#ifndef getsize
void KrasnerTangleData::printSize(std::vector<word64> &s) const {
    s.at(2) += sizeof(boundary_t) * (1 + pairing.capacity()) +
        sizeof(std::vector<boundary_t>) + sizeof(qShift_t);
}

void KrasnerTangle::printSize(std::vector<word64> &s) const {
    data.printSize(s);
}

template <int bitSize>
void KrasnerCoboData<bitSize>::printSize(std::vector<word64> &s) const {
#ifdef USEOLDDOTS
    s.at(5) += sizeof(dots) + sizeof(int16_t) * dots.capacity();
#endif
#ifdef USEZIPDOTS
    s.at(7) += sizeof(nbFacets) + sizeof(compressedDots);
#endif
}

template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::printSize(
        std::vector<word64> &s) const {
    data.printSize(s);
    s.at(6) += sizeof(coeff_tpl);
}
#endif

template <class coeff_tpl, int bitSize>
KrasnerCobo<coeff_tpl, bitSize>::KrasnerCobo(std::ifstream &f,
        bool intCoefficients) {
    if (intCoefficients) {
        uint8_t sgn;
        uint32_t c;
        readFromBinTpl(f, sgn);
        readFromBinTpl(f, c);
        coefficient = std::move(coeff_tpl((sgn ? 1 : -1) * c));
    } else
        coefficient = std::move(coeff_tpl(f));
    boundary_t dotsSize;
    readFromBinTpl(f, dotsSize);
    data.readFromBin(f, dotsSize);
}

KrasnerTangle::KrasnerTangle(std::ifstream &f, boundary_t size) {
    readFromBinTpl(f, qShift);
    {
        qShift_t qShiftTimesN;
        readFromBinTpl(f, qShiftTimesN);
        qShift += N * qShiftTimesN;
    }
    boundary_t circleCount;
    readFromBinTpl(f, circleCount);
    data.setCircleCount(circleCount);
    for (boundary_t i = 0; i < size; ++i) {
        boundary_t newVal;
        readFromBinTpl(f, newVal);
        data.push_back(std::move(newVal));
    }
}

template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::writeToBin(std::ofstream &f) const {
    coefficient.writeToBin(f);
    writeToBinTpl(f, data.dotsSize());
    for (int i = 0; i < data.dotsSize(); ++i)
        writeToBinTpl(f, data.dotsAt(i));
}

void KrasnerTangle::writeToBin(std::ofstream &f) const {
    writeToBinTpl(f, qShift);
    // the times N part of the qShift -- cannot be guessed.
    writeToBinTpl(f, qShift_t(0)); 
    writeToBinTpl(f, data.getCircleCount());
    for (int i = 0; i < data.size(); ++i)
        writeToBinTpl(f, data.at(i));
}


/** Enumerates the facets of a cobordism from domain to codomain,
 * and saves for each boundary point the number of the facet to which it belongs.
 */
template <class coeff_tpl, int bitSize>
boundary_t countFacets(std::vector<boundary_t> &belongsTo,
        boundary_t boundarySize,
        const typename KrasnerCobo<coeff_tpl, bitSize>::tangle_t &domain,
        const typename KrasnerCobo<coeff_tpl, bitSize>::tangle_t &codomain) {

    belongsTo.resize(boundarySize, boundary_t_max);
    boundary_t facetCount = 0;
    for (int i = 0; i < (int)belongsTo.size(); ++i) {
        if (belongsTo.at(i) != boundary_t_max)
           continue;
        bool onTop = false;
        boundary_t j = i;
        do {
            belongsTo.at(j) = facetCount;
            const KrasnerTangle &t = onTop ? codomain : domain;
            j = t.connect(j);
            onTop = ! onTop;
        } while (j != i);
        facetCount += 1;
    }
    return facetCount + domain.getCircleCount() + codomain.getCircleCount();
}

#ifndef notests

template <class coeff_tpl, int bitSize>
bool KrasnerCobo<coeff_tpl, bitSize>::isSane(const Boundary *b,
        const tangle_t *domain, const tangle_t *codomain) const {
    if (b || domain || codomain) {
        assert(b && domain && codomain);

        boundary_t facetCount;
        {
            std::vector<boundary_t> belongsTo;
            facetCount = countFacets<coeff_tpl, bitSize>(belongsTo, b->size(),
                    *domain, *codomain);
        }

        if (facetCount != data.dotsSize()) {
            insane();
            return false;
        }
    }
    return true;
}

bool KrasnerTangle::isSane(const Boundary *b) const {
    if ((b) && (data.size() != b->size())) {
    insane();
        return false;
    }
    if (data.size() % 2) {
    insane();
        return false;
    }
    for (boundary_t idx = 0; idx < data.size(); ++idx) {
        if (data.at(idx) >= data.size()) {
        insane();
            return false;
    }
        if ((data.at(idx) == idx) || (data.at(data.at(idx)) != idx)) {
        insane();
            return false;
    }
    }
    return true;
}


#endif

// Multiplies the given element of the frobenius algebra by x
template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::frobX(std::vector<coeff_tpl> &x) {
    x.insert(x.begin(), coeff_tpl(0));
    for (auto j = frobenius.begin(); j != frobenius.end() - 1; ++j) {
        coeff_tpl toSubstract = x.back();
        toSubstract *= *j;
        toSubstract.switchSign();
        x.at(j - frobenius.begin()) += toSubstract;
    }
    x.pop_back();
}

// Multiplies the given element of the frobenius algebra by a power of x
template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::frobXn(std::vector<coeff_tpl> &x, int n) {
    for (int i = 0; i < n; ++i)
        frobX(x);
}

// Multiplies two elements of the frobenius algebra with each other
template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::frobMult(std::vector<coeff_tpl> &x, const std::vector<coeff_tpl> &y) {
    std::vector<coeff_tpl> tmp(N(), coeff_tpl(0));
    tmp.swap(x);
    for (auto i = y.begin(); i != y.end(); ++i) {
        for (auto j = tmp.begin(); j != tmp.end(); ++j) {
            coeff_tpl toAdd(*j);
            toAdd *= *i;
            x.at(j - tmp.begin()) += toAdd;
        }
        frobX(tmp);
    }
}

template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::frobGenus(std::vector<coeff_tpl> &x) {
    std::vector<coeff_tpl> g(N(), coeff_tpl(0));
// Why is this done every time again? You should make it static!
    for (auto i = g.begin(); i != g.end(); ++i) {
        *i = coeff_tpl(i - g.begin() + 1);
        *i *= frobenius.at(i - g.begin() + 1);
    }
    frobMult(x, g);
}

// The comultiplication of 1 is a linear combination of X^i \otimes X^j. Write it as
// X^0 \oplus p_0 + ... + X^{N-1} \oplus p_{N-1}.
// This returns the product of x with p_d.
template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::frobComult(std::vector<coeff_tpl> &x, int d) {
    std::vector<coeff_tpl> c(N(), coeff_tpl(0));
    for (auto i = c.begin(); i != c.begin() + N() - d; ++i)
        *i = frobenius.at(i - c.begin() + d + 1);
    frobMult(x, c);
}


template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::composeHelper(boundary_t facetCount,
        const KrasnerTangle &middle,
        std::vector<int> &groups, const std::vector<boundary_t> &belongsToNew,
        std::vector<std::vector<coeff_tpl> > &newDots,
        const std::vector<boundary_t> &belongsToLower,
        const KrasnerCoboData<bitSize>& data,
        std::vector<int> &genusCorrection) {

    std::vector<bool> visitedLower(facetCount, false);
    // Traverse every edge of the middle tangle once
    for (int i = 0; i < middle.boundarySize(); ++i) {
        boundary_t j = middle.connect(i);
        if (j > i) continue;
        const int group1 = groups.at(belongsToNew.at(i));
        const int group2 = groups.at(belongsToNew.at(j));
        // If it wasn't clear before that the two facets containing boundary point i and j
        // are connected (i.e. in the same group), connect them now:
        // group1 is deleted, the coefficients multiplied, the genus added.
        if (group1 != group2) {
            std::replace(groups.begin(), groups.end(), group1, group2);
            frobMult(newDots.at(group2), newDots.at(group1));
            genusCorrection.at(group2) += genusCorrection.at(group1) - 4;
            genusCorrection.at(group1) = 0;
        }
        genusCorrection.at(group2) += 1;
        if (! visitedLower.at(belongsToLower.at(i))) {
            visitedLower.at(belongsToLower.at(i)) = true;
            frobXn(newDots.at(group2), data.dotsAt(belongsToLower.at(i)));
            genusCorrection.at(group2) -= 2;
        }
    }
}

int power(int base, int exp) {
    int result = 1;
    for (int i = 0; i < exp; ++i)
        result *= base;
    return result;
}

template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::compose(
        const KrasnerCobo<coeff_tpl, bitSize> &other,
        std::vector<KrasnerCobo<coeff_tpl, bitSize> > &result,
        const tangle_t &lower, const tangle_t &middle,
        const tangle_t &upper) const {

    {
        Boundary b(lower.boundarySize());
        assert(isSane(&b, &lower, &middle));
        assert(other.isSane(&b, &middle, &upper));
    }

    // r is the resulting cobordism
    // in fact, the result may of course be a linear combination of cobordisms
    // in that case, there will be several copies r, which will then be modified
    KrasnerCobo<coeff_tpl, bitSize> r;
    r.coefficient *= coefficient;
    r.coefficient *= other.coefficient;

    // Traverse the spheres which form around the circles of the middle tangle
    for (int i = 0; i < middle.getCircleCount(); ++i) {
        // d is number of dots on the i-th sphere
        int d = (data.dotsAt(data.dotsSize() - middle.getCircleCount() + i)
                      + other.data.dotsAt(other.data.dotsSize()
                      - upper.getCircleCount() - middle.getCircleCount() + i));
        // fewer than (d - 1) dots evaluate to 0 (for all frobenius algebras),
        // and (d - 1) dots evaluate to 1.
        if (d < KrasnerTangle::N - 1)
            return;
        if (d == KrasnerTangle::N - 1)
            continue;
        guaranteeMultVector(d);
        // to evaluate the sphere with d dots, write X^d as sum a_i X^i with
        // i between 0 and N - 1; the result is then a_{N-1}.
        r.coefficient *= multVector.at(d - N()).back();
        if (! r.coefficient.isNonZero())
            return;
    }

    // Count the number of facets of the lower, the upper, and the new
    // (composed) cobordism
    // belongsTo... saves to which facet # each boundary point belongs on each
    // of the three levels.
    std::vector<boundary_t> belongsToNew, belongsToLower, belongsToUpper;
    boundary_t facetCountNew = countFacets<coeff_tpl, bitSize>(belongsToNew,
            lower.boundarySize(), lower, upper)
        - upper.getCircleCount() - lower.getCircleCount();
    boundary_t facetCountLower = countFacets<coeff_tpl, bitSize>(
            belongsToLower, lower.boundarySize(), lower, middle)
        - middle.getCircleCount() - lower.getCircleCount();
    boundary_t facetCountUpper = countFacets<coeff_tpl, bitSize>(
            belongsToUpper, lower.boundarySize(), middle, upper)
        - upper.getCircleCount() - middle.getCircleCount();

    // several new facets may still be joined (will be cut by surgery later).
    // such joined facets form a "group", and groups[i] is the number of the
    // group the i-th facet belongs to
    std::vector<int> groups(facetCountNew);
    for (int i = 0; i < facetCountNew; ++i)
        groups.at(i) = i;

    // This is the coefficient each facet carries (is set to 1 at first, but may
    // be modified by genus).
    std::vector<std::vector<coeff_tpl> > newDots(facetCountNew,
            std::vector<coeff_tpl>(KrasnerTangle::N, coeff_tpl(0)));
    for (auto i = newDots.begin(); i != newDots.end(); ++i)
        i->at(0) = coeff_tpl(1);

    // Facets may have genus; this is now cut, which changes the coefficient
    // Also, newDots is changed if facets are in groups.
    {
        std::vector<int> genusCorrection(facetCountNew, 2);
        composeHelper(facetCountLower, middle, groups, belongsToNew, newDots,
                belongsToLower, data, genusCorrection);
        composeHelper(facetCountUpper, middle, groups, belongsToNew, newDots,
                belongsToUpper, other.data, genusCorrection);
        for (int i = 0; i < facetCountNew; ++i) {
            assert(genusCorrection.at(i) >= 0);
            assert((genusCorrection.at(i) % 4) == 0);
            for (int j = 0; j < (genusCorrection.at(i) / 4); ++j)
                frobGenus(newDots.at(i));
        }
    }

    // needed for bookkeeping
    std::vector<int> lastOfGroup(facetCountNew, 0);
    std::vector<int> groupSizes(facetCountNew, 0);
    int numberOfGroups = 0;
    for (int i = 0; i < facetCountNew; ++i) {
        lastOfGroup.at(groups.at(i)) = i;
        if ((groupSizes.at(groups.at(i)) += 1) == 1)
            numberOfGroups += 1;
    }

    int product = power(KrasnerTangle::N, facetCountNew);

    // this is the "raw" cobordism; for each configuration of dots on r (of
    // which there are product many), the coefficient will now be computed, and
    // a copy of r will be placed into result if the coefficient is non-zero.
    r.data.combine(facetCountNew, lower.getCircleCount(),
            middle.getCircleCount(), upper.getCircleCount(), data, other.data);

    // How does this work exactly? It looks strange...
    const int oldResultSize = result.size();
    result.reserve(oldResultSize + product);
    for (int i = 0; i < product; ++i) {
        std::vector<std::vector<coeff_tpl> > pushVector = newDots;
        // c is the coefficient of the cobordism with this particular
        // configuration of dots.
        coeff_tpl c = coeff_tpl(1);
        for (int j = 0; j < facetCountNew; ++j) {
            const int g = groups.at(j);
            // number of dots on that facet
            const int d = (i / power(KrasnerTangle::N, j))
                                % KrasnerTangle::N;
            r.data.set(j, d);
            // like this, the coefficient of the facet (coming from genus on
            // that facet) is factored in only once; furthermore, for each group
            // member, there is a comultiplication.
            if (j == lastOfGroup.at(g))
                c *= pushVector.at(g).at(d);
            else
                frobComult(pushVector.at(g), d);
        }
        if (c.isNonZero()) {
            result.push_back(r);
            result.back().coefficient *= c;
        }
    }
    {
        Boundary b(lower.boundarySize());
        for (int i = oldResultSize; i < (int)result.size(); ++i)
            assert(result.at(i).isSane(&b, &lower, &upper));
    }
}

template <class coeff_tpl, int bitSize>
bool KrasnerCobo<coeff_tpl, bitSize>::isInvertible(
        const tangle_t & lowerTangle, const tangle_t &upperTangle) const {
    //return (data.dotsSize() == 0);
    if ((lowerTangle != upperTangle) || (lowerTangle.getCircleCount() != 0))
        return false;
    for (int i = 0; i < data.dotsSize(); ++i)
        if (data.dotsAt(i) != 0)
            return false;
    return true;
}

template <class coeff_tpl, int bitSize>
bool KrasnerCobo<coeff_tpl, bitSize>::isEmpty() const {
    return (data.dotsSize() == 0);
}

template <class coeff_tpl, int bitSize>
bool KrasnerCobo<coeff_tpl, bitSize>::operator<(
        const KrasnerCobo<coeff_tpl, bitSize> &other) const {
    return (data < other.data);
}

template <class coeff_tpl, int bitSize>
bool KrasnerCobo<coeff_tpl, bitSize>::operator==(
        const KrasnerCobo<coeff_tpl, bitSize> &other) const {
    return (data == other.data);
}

template <class coeff_tpl, int bitSize>
int KrasnerCobo<coeff_tpl, bitSize>::reducify() {
    assert(data.dotsSize() == 1);
    int result = data.dotsAt(0);
    data.set(0,0);
    return result;
}

/** @post May render itself unusable ("move").
 */
template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::modifyDeloopCopy(int kind, bool left,
        std::vector<KrasnerCobo<coeff_tpl, bitSize> > &v,
        const KrasnerTangle&/* lower*/, const KrasnerTangle& upper)  {

    const int idx = data.dotsSize() - 1 -
        ((left) ? 0 : upper.getCircleCount()); 
    const int d = data.dotsAt(idx);
    guaranteeMultVector(kind + d);
    if ((left && (d == kind)) ||
            ((!left) && (((kind + d + 1) == N()) ||
               (((kind + d + 1) >= N()) &&
                (multVector.at(kind + d - N()).back().isNonZero()))))) {
        data.erase(idx);
        if ((!left) && (kind + d + 1 > N()))
            coefficient *= multVector.at(kind + d - N()).back();
        v.push_back(std::move(*this));
    }
}

void KrasnerTangle::print() const {
    std::cout << "q^(" << qShift << ")";
}

void KrasnerTangle::setToUnion(const KrasnerTangle &o1,
        const KrasnerTangle &o2) {
    const boundary_t glueShift = o1.boundarySize();

    data.setCircleCount(o1.data.getCircleCount() + o2.data.getCircleCount());
    qShift = o1.qShift + o2.qShift;

    data.append(o1.data);
    data.append(o2.data);

    for (boundary_t i = glueShift; i < data.size(); ++i)
        data.set(i, data.at(i) + glueShift);
}

void KrasnerTangle::glue(const boundary_t gluePoints[2]) {
    boundary_t connectedTo[2] = { data.at(gluePoints[0]),
        data.at(gluePoints[1]) };
    if (connectedTo[0] == gluePoints[1])
        data.incCircleCount();
    else
        pairUp(connectedTo[0], connectedTo[1]);
    deleteAndShiftIdx(std::max(gluePoints[0], gluePoints[1]));
    deleteAndShiftIdx(std::min(gluePoints[0], gluePoints[1]));
}

void KrasnerTangle::pairUp(boundary_t x, boundary_t y) {
    data.set(x, y);
    data.set(y, x);
}

void KrasnerTangle::deleteAndShiftIdx(boundary_t x) {
    data.shiftWhatIsHigherK(x);
    data.erase(x);
}

void KrasnerTangle::deloop(std::vector<KrasnerTangle> &newCopies) {
    data.decCircleCount();
    newCopies.resize(KrasnerTangle::N - 1, *this);
    qShift += 1 - KrasnerTangle::N;
    for (int i = 0; i < (int)newCopies.size(); ++i)
        newCopies.at(i).qShift += 2 * i + 3 - KrasnerTangle::N;
}

bool KrasnerTangle::hasLoop() const {
    return (data.getCircleCount() > 0);
}

template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::print() const {
    std::cout << "Coefficient: " << coefficient << ", relative degree: "
        << relativeDegree() << ".";
    for (boundary_t i = 0; i < data.dotsSize(); ++i)
        std::cout << " " << (int)(data.dotsAt(i));
    std::cout << "\n";
}

template <class coeff_tpl, int bitSize>
qShift_t KrasnerCobo<coeff_tpl, bitSize>::degree(boundary_t boundarySize)
        const {
    return relativeDegree() - (boundarySize / 2) * (N() - 1);
}

template <class coeff_tpl, int bitSize>
qShift_t KrasnerCobo<coeff_tpl, bitSize>::relativeDegree() const {
    qShift_t result = 0;
    for (boundary_t i = 0; i < data.dotsSize(); ++i)
        result -= 2 * (int)(data.dotsAt(i));
    return result + data.dotsSize() * (N() - 1);
}

template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::setToUnion(
        const KrasnerTangle& lowerLeft, const KrasnerTangle& lowerRight,
        const KrasnerTangle& upperLeft, const KrasnerTangle& upperRight,
        const KrasnerCobo<coeff_tpl, bitSize> &coboLeft,
        const KrasnerCobo<coeff_tpl, bitSize> &coboRight) {

    coefficient = coboLeft.coefficient;
    coefficient *= coboRight.coefficient;
    data.mashTogether(lowerLeft, lowerRight, upperLeft, upperRight,
            coboLeft.data, coboRight.data);
    assert(data.dotsSize() ==
            coboLeft.data.dotsSize() + coboRight.data.dotsSize());
}

template <int bitSize>
void KrasnerCoboData<bitSize>::mashTogether(const KrasnerTangle& lowerLeft,
        const KrasnerTangle& lowerRight, const KrasnerTangle& upperLeft,
        const KrasnerTangle& upperRight, const KrasnerCoboData &coboLeft,
        const KrasnerCoboData &coboRight) {

    assert(compressedIsCorrect());

    assert(this != &coboLeft);
    assert(this != &coboRight);
    assert(&coboLeft != &coboRight);

    const boundary_t lowerLeftCircles = lowerLeft.data.getCircleCount();
    const boundary_t lowerRightCircles = lowerRight.data.getCircleCount();
    const boundary_t upperLeftCircles = upperLeft.data.getCircleCount();
    const boundary_t upperRightCircles = upperRight.data.getCircleCount();

#ifdef USEOLDDOTS
    dots.insert(dots.end(), coboLeft.dots.begin(),
            coboLeft.dots.end() - lowerLeftCircles - upperLeftCircles);
    dots.insert(dots.end(), coboRight.dots.begin(),
            coboRight.dots.end() - lowerRightCircles - upperRightCircles);
    dots.insert(dots.end(), coboLeft.dots.end() - lowerLeftCircles
            - upperLeftCircles, coboLeft.dots.end() - upperLeftCircles);
    dots.insert(dots.end(), coboRight.dots.end() - lowerRightCircles
            - upperRightCircles, coboRight.dots.end() - upperRightCircles);
    dots.insert(dots.end(), coboLeft.dots.end() - upperLeftCircles,
            coboLeft.dots.end());
    dots.insert(dots.end(), coboRight.dots.end() - upperRightCircles,
            coboRight.dots.end());
#endif

#ifdef USEZIPDOTS
    insertRange(nbFacets, coboLeft, 0,
            coboLeft.nbFacets - lowerLeftCircles - upperLeftCircles);
    insertRange(nbFacets, coboRight, 0,
            coboRight.nbFacets - lowerRightCircles - upperRightCircles);
    insertRange(nbFacets, coboLeft,
            coboLeft.nbFacets - lowerLeftCircles - upperLeftCircles,
            coboLeft.nbFacets - upperLeftCircles);
    insertRange(nbFacets, coboRight,
            coboRight.nbFacets - lowerRightCircles - upperRightCircles,
            coboRight.nbFacets - upperRightCircles);
    insertRange(nbFacets, coboLeft,
            coboLeft.nbFacets - upperLeftCircles, coboLeft.nbFacets);
    insertRange(nbFacets, coboRight,
            coboRight.nbFacets - upperRightCircles, coboRight.nbFacets);
#endif

    assert(compressedIsCorrect());
}

/** @post: Dots now look like this:
           facetCountNew many 0s, then the facets bounding lower Circles (caps),
           then the facets bounding upper circles (cups).
 */
template <int bitSize>
void KrasnerCoboData<bitSize>::combine(boundary_t facetCountNew,
        boundary_t lowerCircleCount, boundary_t middleCircleCount,
        boundary_t upperCircleCount, const KrasnerCoboData &first,
        const KrasnerCoboData &other) {

    assert(compressedIsCorrect());
#ifndef NDEBUG
    const int newNbDots = (int)facetCountNew + upperCircleCount
        + lowerCircleCount;
#endif
    assert(newNbDots * bitsPerDot <= bitSize);

#ifdef USEOLDDOTS
    dots.reserve((int)facetCountNew + upperCircleCount + lowerCircleCount);
    dots.resize(facetCountNew, 0);
    dots.insert(dots.end(),
        first.dots.end() - middleCircleCount - lowerCircleCount,
            first.dots.end() - middleCircleCount);
    dots.insert(dots.end(), other.dots.end() - upperCircleCount,
            other.dots.end());
#endif

#ifdef USEZIPDOTS
    insertZeroes(0, facetCountNew);
    insertRange(nbFacets, first, first.nbFacets - middleCircleCount -
            lowerCircleCount, first.nbFacets - middleCircleCount);
    insertRange(nbFacets, other, other.nbFacets - upperCircleCount,
            other.nbFacets);
#endif

    assert(compressedIsCorrect());
}

template <class coeff_tpl, int bitSize>
void KrasnerCobo<coeff_tpl, bitSize>::glue(const KrasnerTangle& lower,
        const KrasnerTangle& upper, const boundary_t gluePoints[2],
        boundary_t boundarySize,
        LCCobos<KrasnerCobo<coeff_tpl, bitSize> > &father) {

    std::deque<boundary_t> minima;
    boundary_t idx[2] = { boundary_t_max, boundary_t_max };
    boundary_t minimaAfterDeletion[2] = { boundary_t_max, boundary_t_max };
    bool circleIsUpper[2] = { false, false };
    /* The following block determines:
     * - idx[k] is the index of the facet containing the k-th GluePoint
     * - minima is a list of the boundary point with the lowest number
     *   contained by each facet that is not a gluepoint,
     *   and boundary_t_max if no such point exists
     * - minimaAfterDeletion is only used if one facet is glued to itself; in
     *   that case, it contains the minimum boundary point of the two halves of
     *   the facet boundary (and bondary_t_max if the half is empty,
     *   respectively).
     * - circleIsUpper is only used if one facet is glued to itself; if
     *   additionally, the k-th half of the facet gives a circle, it says
     *   whether this circle is in the upper tangle.
     */
    {
        std::vector<bool> visited(boundarySize, false);
        int i = 0;
        while (i < (int)visited.size()) {
            boundary_t minimum = boundary_t_max;
            boundary_t minimaAfterDeletionTmp[2] =
                { boundary_t_max, boundary_t_max };
            boundary_t j = i;
            bool onTop = false, is[2];
            is[0] = is[1] = false;
            do {
                visited.at(j) = true;
                if ((j != gluePoints[0]) && (j != gluePoints[1])) {
                    if (j < minimum)
                        minimum = j;
                    const int l = (is[0] == is[1]) ? 0 : 1;
                    if (j < minimaAfterDeletionTmp[l])
                        minimaAfterDeletionTmp[l] = j;
                } else
                    for (int k = 0; k < 2; ++k)
                        if (j == gluePoints[k]) {
                            idx[k] = minima.size();
                            is[k] = true;
                        }
                const boundary_t oldJ = j;
                const KrasnerTangle &t = onTop ? upper : lower;
                j = t.connect(j);
                if (onTop && (((oldJ == gluePoints[1]) && (j == gluePoints[0]))
                          || ((oldJ == gluePoints[0]) &&
                              (j == gluePoints[1]))))
                    circleIsUpper[(is[0] == is[1]) ? 0 : 1] = true;
                onTop = ! onTop;
            } while (j != i);
            if (is[0] && is[1]) {
                minimaAfterDeletion[0] = minimaAfterDeletionTmp[0];
                minimaAfterDeletion[1] = minimaAfterDeletionTmp[1];
            }
            minima.push_back(minimum);
            while ((i < (int)visited.size()) && visited.at(i))
                ++i;
        }
    }

    //!!i
    if ((idx[0] == boundary_t_max) || (idx[1] == boundary_t_max))
        throw;
    /* Glowing a facet to itself (if), or to another one (else). */
    if (idx[0] == idx[1]) {
        int d = data.dotsAt(idx[0]);
        data.erase(idx[0]);
        minima.erase(minima.begin() + idx[0]);
        boundary_t newIdx[2] = { 0, 0 };
        for (int k = 0; k < 2; ++k) {
            if (minimaAfterDeletion[k] == boundary_t_max)
                newIdx[k] = circleIsUpper[k] ? (data.dotsSize()) :
                    (data.dotsSize() - upper.getCircleCount());
            else while ((newIdx[k] < (int)minima.size())
                    && (minima.at(newIdx[k])
                        < minimaAfterDeletion[k]))
                newIdx[k] += 1;
        }
        if (newIdx[0] > newIdx[1])
            std::swap(newIdx[0], newIdx[1]);
        data.insert(newIdx[1], 0);
        data.insert(newIdx[0], 0);

        coeff_tpl oldCoefficient = coefficient;
        for (typename std::vector<coeff_tpl>::const_iterator j =
                frobenius.begin(); j != frobenius.end(); ++j) {
            const int jb = j - frobenius.begin();
            if ((jb == d) || (!j->isNonZero()))
                continue;
            coefficient *= *j;
            if (jb < d)
                coefficient.switchSign();
            for (int i = std::min(jb, d); i < std::max(jb, d); ++i) { 
                KrasnerCobo<coeff_tpl, bitSize> x = *this;
                x.data.set(newIdx[0], i);
                x.data.set(newIdx[1] + 1, jb + d - 1 - i);

                father.add(x);
            }
            coefficient = oldCoefficient;
        }
    } else {
        int d = data.dotsAt(idx[0]) + data.dotsAt(idx[1]);
        data.erase(std::max(idx[0], idx[1]));
        data.erase(std::min(idx[0], idx[1]));
        boundary_t key = std::min(minima.at(idx[0]), minima.at(idx[1]));
        assert (key != boundary_t_max);
        minima.erase(minima.begin() + std::max(idx[0], idx[1]));
        minima.erase(minima.begin() + std::min(idx[0], idx[1]));
        boundary_t newIdx = 0;
        while ((newIdx < (int)minima.size()) && (minima.at(newIdx) < key))
            newIdx += 1;
        if (d < KrasnerTangle::N) {
            data.insert(newIdx, d);
            father.add(*this);
        } else {
            guaranteeMultVector(d);
            for (typename std::vector<coeff_tpl>::const_iterator i =
                    multVector.at(d - N()).begin();
                    i != multVector.at(d - N()).end(); ++i)
                if (i->isNonZero()) {
                    KrasnerCobo<coeff_tpl, bitSize> x = *this;
                    x.data.insert(newIdx, i -
                        multVector.at(d - N()).begin());
                    x.coefficient *= *i;
                    father.add(x);
                }
        }
    }
}
