/*
 *
 *    src/planar_algebra/smith.cpp --- Part of khoca, a knot homology calculator
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
#include <vector>
#include <deque>
#include <bitset>
#include <stdlib.h>
#include <inttypes.h>
#include <sstream>
#include <algorithm>
// Next line is necessary due to a bug of mpir
#include <stddef.h>
#include <gmp.h>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <pari/pari.h>

#include "../shared.h"

#include <assert.h>

#include "coefficient_rings.h"
#include "sparsemat.h"
#include "planar_algebra.h"
#include "../krasner/krasner.h"

/** In some pari versions, it's called "CATCH", in others "pari_CATCH"...
  */
#ifndef pari_CATCH
#define pari_CATCH CATCH
#endif
#ifndef pari_ENDCATCH
#define pari_ENDCATCH ENDCATCH
#endif
#ifndef pari_TRY
#define pari_TRY TRY
#endif

#define INSTANTIATE_META(BITSIZE) \
template void calculateSmithFriend(const Complex<KrasnerCobo<MInteger, BITSIZE> >&, std::ostream&, int progress);\
template void calculateSmithFriend(const Complex<KrasnerCobo<MRational, BITSIZE> >&, std::ostream&, int progress);\
template void calculateSmithFriend(const Complex<KrasnerCobo<FF<uint8_t>, BITSIZE> >&, std::ostream&, int progress);\
template void calculateSmithFriend(const Complex<KrasnerCobo<FF<uint16_t>, BITSIZE> >&, std::ostream&, int progress);\
template void calculateSmithFriend(const Complex<KrasnerCobo<Polynomial<MInteger>, BITSIZE> >&, std::ostream&, int progress);\
template void calculateSmithFriend(const Complex<KrasnerCobo<Polynomial<MRational>, BITSIZE> >&, std::ostream&, int progress);\
template void calculateSmithFriend(const Complex<KrasnerCobo<Polynomial<FF<uint8_t> >, BITSIZE> >&, std::ostream&, int progress);\
template void calculateSmithFriend(const Complex<KrasnerCobo<Polynomial<FF<uint16_t> >, BITSIZE> >&, std::ostream&, int progress);

MULTICALLER

#undef INSTANTIATE_META


bool withPari();
GEN mpz2GEN(mpz_t X);

bool withPari() {
    return true;
}


/** Variables and functions to protect the Pari-stack of external usage
 */
bool backup_needed  = true; // false if cypari / cypari2 don't use the same Pari stack
void* mem_block; // pointer to save memory block of externally used Pari stack
size_t mem_block_size = 0; // size of externally used Pari stack
pari_sp av = 0; // avma of Pari-stack from external usage

void pari_backup() {
    av = avma;
    if (av == 0) {backup_needed = false;}
    if (backup_needed) {
        mem_block_size = (size_t) (pari_mainstack->top - pari_mainstack->bot);
        mem_block = malloc(mem_block_size);
        memcpy(mem_block, (void*) pari_mainstack->bot, mem_block_size);
    }
}

void pari_rollback() {
    if (backup_needed) {
        avma = av;
        memcpy((void*) pari_mainstack->bot, mem_block, mem_block_size);
        free(mem_block);
    }
}

void init_pari(word64 pariStackSize) {
    if (backup_needed) {
        size_t size = std::max((size_t) pariStackSize, pari_mainstack->rsize);
        size_t maxsize = std::max((size_t) pariStackSize, pari_mainstack->vsize);
        paristack_setsize(size, maxsize);
    }
    else {
        pari_init(pariStackSize, 0);
    }
}

void close_pari() {
    if (!backup_needed) {
        pari_close();
    }
}

/** Converts an MPIR-arbitrary precision integer to a PARI-arbitrary precision
 * integer.
 */
GEN mpz2GEN(const mpz_t X) {
    int l = X->_mp_size;
    int lx = labs(l)+2;
    GEN x = cgeti(lx);
    x[1] = evalsigne(l > 0? 1: -1) | evallgefint(lx);
    for (int i = 2; i < (int)lx; i++)
        x[i] = X->_mp_d[i-2];
    return x;
}

template <class cobordism_tpl>
void calculateSmithFriend(const Complex<cobordism_tpl>& that,
        std::ostream& s, int progress) {
    /* Get minimal and maximal q-degree */
    std::vector<std::vector<int> > idxTranslators, qtDims;
    idxTranslators.resize(that.vecTangles.size());
    qtDims.resize(that.vecTangles.size());
    std::vector<qShift_t> qMins, qMaxs;
    qMins.resize(that.vecTangles.size());
    qMaxs.resize(that.vecTangles.size());
    for (auto k = that.vecTangles.cbegin(); k != that.vecTangles.cend(); ++k) {
        k->getQs(qMaxs.at(k - that.vecTangles.begin()),
                qMins.at(k - that.vecTangles.begin()),
                idxTranslators.at(k - that.vecTangles.begin()),
                qtDims.at(k - that.vecTangles.begin()));
        qMaxs.at(k - that.vecTangles.begin());
        qMins.at(k - that.vecTangles.begin());
    }
    s << "[";
    const char* comma = "";
    for (int i = 0; i < (int)qtDims.size(); ++i)
        for (int j = 0; j < (int)qtDims.at(i).size(); ++j) {
            s << comma << "[" << (int)i + that.globalTShift << "," <<
                2 * (int)j + qMins.at(i) << ",0," << qtDims.at(i).at(j) << "]";
            comma = ",";
        }

    static word64 pariStackSize;
    pariStackSize = 16 * 1024 * 1024;
    while (true) {
        init_pari(pariStackSize);
        pari_CATCH(CATCH_ALL) {
            close_pari();
            pariStackSize *= 2;
            std::cerr << "Pari stack overflows. Doubling stack to "
                << pariStackSize / (1024 * 1024) << "MB and retrying.\n";
        } pari_TRY {
            int numberOfSnf = 0;
            for (int i = 0; i < (int)that.matLCCobos.size(); ++i)
                if ((qtDims.at(i).size() != 0) &&
                        (qtDims.at(i + 1).size() != 0))
                    numberOfSnf += 1 + (qMaxs.at(i) - qMins.at(i)) / 2;
            std::stringstream ssLoc;
            int snfCalculated = 0;
            for (int i = 0; i < (int)that.matLCCobos.size(); ++i) {
                if ((qtDims.at(i).size() == 0) ||
                        (qtDims.at(i + 1).size() == 0))
                    continue;
                std::vector<GEN> vSM;
                vSM.reserve(1 + (qMaxs.at(i) - qMins.at(i)) / 2);
                for (int j = 0; j < 1 + (qMaxs.at(i) - qMins.at(i)) / 2; ++j)
                    if (((2*j + qMins.at(i)) <= qMaxs.at(i + 1)) &&
                            ((2*j + qMins.at(i)) >= qMins.at(i + 1)))
                        vSM.push_back(zeromatcopy(qtDims.at(i + 1).at(
                                        j + (qMins.at(i) - qMins.at(i + 1))/2),
                                    qtDims.at(i).at(j)));
                    else
                        vSM.push_back(zeromatcopy(0, 0));
                for (typename Complex<cobordism_tpl>::matLCCobos_t::\
                        constIterator_t j = that.matLCCobos.at(i).\
                        getConstIterator(); j.isOn(); j.stepAlongMat()) {
                    qShift_t q = ( that.vecTangles.at(i).getTangles().at(
                                j.getCol()).getQShift() - qMins.at(i)) / 2;
                    if (2 * q + qMins.at(i) != that.vecTangles.at(
                                i + 1).getTangles().at(j.getRow()).getQShift())
                        continue;
                    GEN mft = mpz2GEN(*(j.getVal()->getCoeff()->getmpz_t()));
                    gcoeff(vSM.at(q), 1 + idxTranslators.at(i + 1).at(
                                j.getRow()), 1 + idxTranslators.at(i).at(
                                j.getCol())) = mft;
                }
                for (int j = 0; j < 1 + (qMaxs.at(i) - qMins.at(i)) / 2; ++j) {
                    int r = 0;
                    if (progress) {
                        if (((2*j + qMins.at(i)) <= qMaxs.at(i + 1)) &&
                                ((2*j + qMins.at(i)) >= qMins.at(i + 1)))
                            std::cerr << " Calculating Smith Normal form at t = "
                                << (int)i - that.globalTShift << ", q = "
                                << 2 * (int)j + qMins.at(i) << " (" << qtDims.at(
                                            i + 1).at(j + (qMins.at(i) - qMins.at(
                                                        i + 1))/2) << " x "
                                            << qtDims.at(i).at(j) << ")\033[K";
                        std::cerr << "\n";
                    }
                    GEN y = ZM_snf(vSM.at(j));
                    snfCalculated += 1;
                    if (progress) {
                        int percent = 100 * snfCalculated / numberOfSnf;
                        std::cerr << " " << std::setw(3) << percent <<
                            "% [" << std::string(percent/2, '*') <<
                            std::string(50 - percent/2,'.') << "] " << snfCalculated
                            << "/" << numberOfSnf << "\033[K\x1b[A\r" << std::flush;
                    }
                    
                    for (int k = lg(y) - 1; (k >= 1) &&
                            (signe(gel(y, k)) != 0); --k) {
                        ++r;
                        GEN factors = Z_factor(gel(y, k));
                        for (int x = 1; x < lg(gel(factors,1)); ++x) {
                            int primePowerFactor =
                                itos(powii(gmael2(factors, 1, x),
                                            gmael2(factors, 2, x)));
                            ssLoc << ",[" << (int)i + 1 + that.globalTShift
                                << "," << 2 * j + qMins.at(i) << ","
                                << primePowerFactor << "," << 1 << "]";
                        }
                    }
                    // dim kernel...
                    ssLoc << ",[" << (int)i + that.globalTShift << ","
                        << 2 * j + qMins.at(i) << ",0," << - (int)r << "]";
                    // ...minus dim image
                    ssLoc << ",[" << (int)i + 1 + that.globalTShift << ","
                        << 2 * j + qMins.at(i) << ",0," << -(int)r << "]";
                }
            }
            s << ssLoc.str();
            break;
        } pari_ENDCATCH
    }
    if (progress)
        std::cerr << "\n\n";
    close_pari();
    s << "]";
}
