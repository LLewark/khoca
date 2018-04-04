/*
 *
 *    src/python_interface/pythonInterface.h --- Part of khoca, a knot homology calculator
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

class ComplexStack {
    public:
        // If F is empty, do equivariant of sl_N.
        ComplexStack(int mod_, std::vector<int> F, int N, int, int);
        ~ComplexStack();
#ifndef getsize
        void outputTotalSize() const;
#endif
        void printHomology(int idx);
        void calculateHomology(int idx, std::string& result);
        void calculateIntegralHomology(int idx, std::string& result, int progress);
        bool guaranteeSize(int size);
	/** Loads a complex from file (formatted as specified), and appends it
         * to the stack. */
        int loadComplexFromFile(int idx, std::string fileName, int fileFormat);
        void deleteComplex(int idx);
        int tensorComplexes(int idx, int firstIdx, int secondIdx, int progress);
	void deleteNonIsos(int idx);
        void glueComplex(int idx, int gluePoint1, int gluePoint2);
        int simplifyComplexOnce(int idx, int numThreads, int progress);
        int simplifyComplexParallely(int idx, int numThreads, int progress);
        void saveComplexToFile(int idx, std::string fileName,
                int /*fileFormat*/) const;
        void setRoot(int r);
	int copyComplex(int fromIdx, int toIdx);
        int dualizeComplex(int fromIdx, int toIdx);
	int firstFreeIdx();
	void reducify(int idx);

	void resetSimplificationsCounter(int idx);
        void stepPage();
	void resetPage();
        int getPage() const;

        void printCompileInfo();
        void outputDetailed(int idx);
    private:
        void startThread(int numJobs, int* status, int idx, int id,
                bool* changed, bool* done, int numThreads, int progress);
        int allDone(int idx);
        int doneAtTDegree(int idx, int tDegree);

        int mod;
        /** Is only used to copy from etc.
         */
        void *tokenComplex;
        int page;
        int root;
        std::deque<void*> complexStack;
};
