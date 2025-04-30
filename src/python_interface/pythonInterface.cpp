/*
 *
 *    src/python_interface/pythonInterface.cpp --- Part of khoca, a knot homology calculator
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

/** @file pythonInterface.cpp
    @brief Is called by pui.pyx, calls planar_algebra.cpp

There is a stack of complexes, which is manipulated.
Everything can be multithreaded (in the future).
*/

#include <thread>
#include <mutex>
#include <condition_variable>

#include <exception>
#include <iostream>
#include <vector>
#include <deque>
#include <bitset>
#include <inttypes.h>
#include <string>
#include <sstream>
#include <stdlib.h>
#include <algorithm>
// stddef.h included for mpir
#include <stddef.h>
#include <gmp.h>
#include <fstream>
#include <math.h>

#include "../shared.h"
#include <assert.h>

#include "../planar_algebra/coefficient_rings.h"
#include "../planar_algebra/sparsemat.h"
#include "../planar_algebra/planar_algebra.h"
#include "../planar_algebra/smith.h"
#include "../krasner/krasner.h"
#include "pythonInterface.h"

#include "../compilerFlagsInfo.cpp"

std::mutex statusMtx;
std::condition_variable thereIsWorkCv;

int activeThreads = 0;

std::ostream * io::cprogress_s = &std::cerr;
std::string io::markUp = "\033[1m";
std::string io::noMarkUp = "\033[m";

void ComplexStack::outputDetailed(int idx) {
    ((AbstractComplex*)complexStack.at(idx))->detailedOutput(std::cout);
}

void ComplexStack::printCompileInfo() {
    std::cerr << COMPILERFLAGSINFO << " version. ";
#ifdef NDEBUG
    std::cerr << "Assertions are off. ";
#else
    std::cerr << "Assertions are on. ";
#endif
#ifdef notests
    std::cerr << "Tests are off. ";
#else
    std::cerr << "Tests are on. ";
#endif
#ifdef getsize
    std::cerr << "Size functions are off. ";
#else
    std::cerr << "Size functions are on. ";
#endif
#ifdef USEOLDDOTS
    std::cerr << "Old dots are on. ";
#else
    std::cerr << "Old dots are off. ";
#endif
#ifdef USEZIPDOTS
    std::cerr << "New dots are on. ";
#else
    std::cerr << "New dots are off. ";
#endif
    std::cerr << "Supported bitsizes:";
#define INSTANTIATE_META(b) std::cerr << " " << b;
    MULTICALLER
#undef INSTANTIATE_META
    std::cerr << ".\n";
}

ComplexStack::ComplexStack(int mod_, std::vector<int> F, int N, int girth, int verbose) :
    mod(mod_), tokenComplex(nullptr), page(0) {

    if ((mod > 1) && (mod < 256))
        FF<uint8_t>::setP(mod);
    else if ((mod >= 256) && (mod < 65536))
        FF<uint16_t>::setP(mod);
    else if ((mod != 0) && (mod != 1)) {
        std::cerr << "Parameter mod is " << mod_ << ", but currently, only "
            "finite fields with up to 2^16 elements are supported.";
        throw;
    }

    const bool equivariant = (F.size() == 0);
    KrasnerTangle::N = equivariant ? N : F.size();

    /** This dirty bit of code creates the tokenComplex with the smallest
      supported bitsite that is big enough.*/
    const int bitsPerDot = (int)ceil(log2(KrasnerTangle::N));
    const int bitSize = ((bitsPerDot * (2 * girth) + 7) / 8) * 8;
    if (verbose)
        std::cerr << "Each dot needs " << bitsPerDot << " bits; maximal number" \
        " of facets is " << (2 * girth) << ", so minimal bitSize is " <<
        bitSize << ". ";

#define INSTANTIATE_META(BITSIZE) \
    if ((tokenComplex == nullptr) && (BITSIZE >= bitSize)) { \
        if (verbose) \
            std::cerr << "Taking bitSize " << BITSIZE << ".\n"; \
        KrasnerCoboData<BITSIZE>::setBitsPerDot(bitsPerDot); \
        if (equivariant) { \
            if (mod == 0) \
                tokenComplex = (void *)(new \
                Complex<KrasnerCobo<Polynomial<MInteger>, BITSIZE> >); \
            else if (mod == 1) \
                tokenComplex = (void *)(new \
                Complex<KrasnerCobo<Polynomial<MRational>, BITSIZE> >); \
            else if (mod < 256) \
                tokenComplex = (void *)(new \
                Complex<KrasnerCobo<Polynomial<FF<uint8_t> >, BITSIZE > >); \
            else \
                tokenComplex = (void *)(new \
                Complex<KrasnerCobo<Polynomial<FF<uint16_t> >, BITSIZE > >); \
        } else { \
            if (mod == 0) \
                tokenComplex = (void *)(new \
                        Complex<KrasnerCobo<MInteger, BITSIZE> >); \
            else if (mod == 1) \
                tokenComplex = (void *)(new \
                        Complex<KrasnerCobo<MRational, BITSIZE> >); \
            else if (mod < 256) \
                tokenComplex = (void *)(new \
                        Complex<KrasnerCobo<FF<uint8_t>, BITSIZE > >); \
            else \
                tokenComplex = (void *)(new \
                        Complex<KrasnerCobo<FF<uint16_t>, BITSIZE > >); \
        } \
    }

    MULTICALLER
#undef INSTANTIATE_META

    if (tokenComplex == nullptr) {
        std::cerr << "bitSize " << bitSize << " exceeds limit.\n";
        throw;
    }

    pari_backup();
    ((AbstractComplex*)tokenComplex)->initialiseFrobenius(F, N);
    std::cout << "Frobenius algebra: ";
    ((AbstractComplex*)tokenComplex)->printFrobenius(std::cout);
    std::cout << "." << std::endl;
}

#ifndef getsize
void ComplexStack::outputTotalSize() const {
    std::vector<word64> s(8, 0);
    for (std::deque<void*>::const_iterator i = complexStack.cbegin();
            i != complexStack.cend(); ++i)
        if (*i)
            ((AbstractComplex*)(*i))->printSize(s);
    std::cout << "\nSizedata:";
    for (int i = 0; i < 8; ++i)
        std::cout << s.at(i) << "\t";
    std::cout << "\n";
}
#endif

ComplexStack::~ComplexStack() {
    pari_rollback();
    delete ((AbstractComplex*)tokenComplex);
    for (auto i = complexStack.begin(); i != complexStack.end(); ++i)
        deleteComplex(i - complexStack.begin());
}

void ComplexStack::printHomology(int) {
}

void ComplexStack::calculateHomology(int idx, std::string& result) {
    printHomology(idx);
    std::stringbuf buffer;
    std::ostream os(&buffer);
    os << *((AbstractComplex*)complexStack.at(idx));
    result = buffer.str();
}

void ComplexStack::calculateIntegralHomology(int idx, std::string& result, int progress) {
    assert(mod == 0);
    std::stringbuf buffer;
    std::ostream os(&buffer);
    ((AbstractComplex*)(complexStack.at(idx)))->calculateSmith(os, progress);
    result = buffer.str();
}

bool ComplexStack::guaranteeSize(int size) {
    if ((int)complexStack.size() < size)
        complexStack.resize(size, nullptr);
    return true;
}

int ComplexStack::loadComplexFromFile(int idx, std::string fileName,
        int fileFormat) {
    if (! guaranteeSize(idx + 1))
        return 1;
    if (fileFormat == 0) {
	std::ifstream f(fileName, std::ios::binary);
        complexStack.at(idx) =
            (void*)(((AbstractComplex*)tokenComplex)->loadFromFile(f));
    } else return 2;
    return 0;
}

void ComplexStack::deleteComplex(int idx) {
    delete (AbstractComplex*)complexStack.at(idx);
    complexStack.at(idx) = nullptr;
}

int ComplexStack::tensorComplexes(int idx, int firstIdx, int secondIdx, int progress) {
    if (progress)
        io::cprogress() << io::markUp << "T(" << firstIdx << "," << secondIdx
            << ")->" << idx<< ": " << io::noMarkUp;
    if (! guaranteeSize(idx + 1))
        return 1;
    complexStack.at(idx) = (void*)(((AbstractComplex*)tokenComplex)->tensor(
                (AbstractComplex*)complexStack.at(firstIdx),
                (AbstractComplex*)complexStack.at(secondIdx)));
#ifndef PRINTCOMPLEX
    ((AbstractComplex*)complexStack.at(idx))->detailedOutput(std::cout);
#endif
#ifndef getsize
    outputTotalSize();
#endif
    return 0;
}

void ComplexStack::glueComplex(int idx, int gluePoint1, int gluePoint2) {
    boundary_t gluePoints[2] =
        { (boundary_t)gluePoint1, (boundary_t)gluePoint2 };
    ((AbstractComplex*)complexStack.at(idx))->glue(gluePoints);
#ifndef getsize
    outputTotalSize();
#endif
}

void ComplexStack::reducify(int idx) {
    ((AbstractComplex*)complexStack.at(idx))->reducify(root);
}

void ComplexStack::deleteNonIsos(int idx) {
    ((AbstractComplex*)complexStack.at(idx))->deleteNonIsos();
}

int ComplexStack::simplifyComplexOnce(int idx, int numThreads, int progress) {
    int success;
    success = ((AbstractComplex*)(complexStack.at(idx)))->simplifyOnce(
            (int)page * (-2), numThreads, progress);
    if ((success != 0) && (progress))
        io::cprogress() << "\n";
    return success;
}

void advanceToNextFree(bool *done, int numJobs, int &myJob, int* status) {
    while ((myJob < numJobs) && ((status[myJob]) || done[myJob]))
        myJob += 1;
}

int ComplexStack::allDone(int idx) {
    return ((AbstractComplex*)complexStack.at(idx))->isSimplified(
            (int)page * (-2));
}

int ComplexStack::doneAtTDegree(int idx, int tDegree) {
    return ((AbstractComplex*)complexStack.at(idx))->isSimplified(
            (int)page * (-2), tDegree);
}

bool allTrue(bool* done, int numJobs) {
    for (int i = 0; i < numJobs; ++i)
        if (! done[i])
            return false;
    return true;
}

void ComplexStack::startThread(int numJobs, int* status,
        int idx, int id, bool* changed, bool* done, int numThreads, int progress) {
    while (true) {
        int myJob;
        /* Find a job or go to sleep.
           If there is more than one job, awaken somebody else.
         */
        while (true) {
            std::unique_lock<std::mutex> lck (statusMtx);
            if (allTrue(done, numJobs)) {
                thereIsWorkCv.notify_all();
                return;
            } 
            myJob = 0;
            advanceToNextFree(done, numJobs, myJob, status);
            if (myJob < numJobs) {
                if (myJob > 1) status[myJob - 2] = 1;
                if (myJob > 0) status[myJob - 1] = 1;
                status[myJob] = id + 2;
                if (myJob < numJobs - 1) status[myJob + 1] = 1;
                if (myJob < numJobs - 2) status[myJob + 2] = 1;
                {
                    int mySecondJob = myJob + 2;
                    advanceToNextFree(done, numJobs, mySecondJob, status);
                    if (mySecondJob < numJobs)
                        thereIsWorkCv.notify_one();
                }
                activeThreads += 1;
                break;
            }
            thereIsWorkCv.wait(lck);
        }
        // Do the job.
        int success = 0;
        while (success == 0) {
            std::stringstream ss;
            ss << "  {";
            for (int i = 0; i < numJobs; ++i) {
                if (done[i]) ss << "_";
                else if (status[i] >= 2) ss << (char)(63 + status[i]);
                else if (status[i] == 1) ss << "-";
                else if (status[i] == 0) ss << ".";
            }
            ss << "}, T=" << activeThreads << ".";
            std::string myStr = ss.str();
            success = ((AbstractComplex*)complexStack.at(idx))->\
                simplifyOnceAtTDegree((int)page * (-2), myJob, &myStr,
                        numThreads, progress);
            if (success == 0) {
#ifndef PRINTCOMPLEX
                ((AbstractComplex*)complexStack.at(idx))->detailedOutput(std::cout);
#endif
                *changed = true;
            }
        }
        // Log what you did and free yourself.
        {
            std::unique_lock<std::mutex> lck (statusMtx);
            activeThreads -= 1;

            if ((myJob > 1) && ((myJob < 3) || (status[myJob - 3] < 2)) &&
                    ((myJob < 4) || (status[myJob - 4] < 2)))
                status[myJob - 2] = 0;
            if ((myJob > 0) && ((myJob < 2) || (status[myJob - 2] < 2)) &&
                    ((myJob < 3) || (status[myJob - 3] < 2)))
                status[myJob - 1] = 0;
            status[myJob] = 0;
            if ((myJob < numJobs - 1) && ((myJob > numJobs - 3) ||
                        (status[myJob + 2] < 2)) && ((myJob > numJobs - 4) ||
                            (status[myJob + 3] < 2)))
                status[myJob + 1] = 0;
            if ((myJob < numJobs - 2) && ((myJob > numJobs - 4) ||
                        (status[myJob + 3] < 2)) && ((myJob > numJobs - 5) ||
                            (status[myJob + 4] < 2)))
                status[myJob + 2] = 0;

            // Gaussing doesn't undo neighbours; delooping may undo
            // (by creating an iso) the predecessor.
            if (myJob > 0) done[myJob - 1] = false;
            done[myJob] = true;
        }
    }
}

int ComplexStack::simplifyComplexParallely(int idx, int numThreads, int progress) {
    int numJobs;
    numJobs = ((AbstractComplex*)complexStack.at(idx))->size();

    // Arrays are used because several threads may write to different
    // indices of one array at the same time.
    int* status = new int[numJobs];
    bool* done = new bool[numJobs];

    bool changed = false;
    std::fill_n(status, numJobs, 0);
    std::fill_n(done, numJobs, false);
    /*
 1: only higher-degree isos left
 2: all matrices are zero
 4: only higher-degree isos left && nothing changed
 5: all matrices are zero && nothing changed
*/
    std::thread* t = new std::thread[numThreads];
    activeThreads = 0;
    for (int i = 0; i < numThreads; ++i)
        t[i] = std::thread(&ComplexStack::startThread, this, numJobs,
                (int*)status, idx, i, &changed, (bool *)done, numThreads, progress);
    for (int i = 0; i < numThreads; ++i)
        t[i].join();
    if (progress)
        io::cprogress() << "\n";

    delete[] status;
    delete[] done;
    delete[] t;

    return allDone(idx) + (changed ? 0 : 3);
}

void ComplexStack::saveComplexToFile(int idx, std::string fileName,
        int /*fileFormat*/) const {
    std::ofstream f(fileName, std::ofstream::binary);
    ((AbstractComplex*)complexStack.at(idx))->writeToBin(f);
    f.close();
}

void ComplexStack::setRoot(int r) {
    root = r;
}

void ComplexStack::stepPage() {
    page += 1;
}

void ComplexStack::resetSimplificationsCounter(int idx) {
    ((AbstractComplex*)complexStack.at(idx))->simplificationsCount = 0;
}

int ComplexStack::getPage() const {
    return page;
}

int ComplexStack::firstFreeIdx() {
    for (auto i = complexStack.begin(); i != complexStack.end(); ++i)
        if (*i == nullptr)
            return (i - complexStack.begin());
    guaranteeSize(complexStack.size() + 1);
    return complexStack.size() - 1;
}

int ComplexStack::dualizeComplex(int fromIdx, int toIdx) {
    if (! guaranteeSize(toIdx + 1))
        return 1;
    complexStack.at(toIdx) = (void*)(((AbstractComplex*)tokenComplex)->\
            setToDual((AbstractComplex*)complexStack.at(fromIdx)));
    return 0;
}

int ComplexStack::copyComplex(int fromIdx, int toIdx) {
    if (! guaranteeSize(toIdx + 1))
        return 1;
    complexStack.at(toIdx) = (void*)(((AbstractComplex*)tokenComplex)->\
            copy((AbstractComplex*)complexStack.at(fromIdx)));
    return 0;
}

void ComplexStack::resetPage() {
    page = 0;
}
