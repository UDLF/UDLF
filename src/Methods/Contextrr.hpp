/* <Contextrr.hpp>
 *
 * Contextrr method class header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 ******************************************************************************************
 *
 * ContextRR is presented in the paper:
 *   D. C. G. Pedronette and R. da S. Torres.
 *   "Exploiting contextual information for image re-ranking."
 *   Iberoamerican Congress on Pattern Recognition (CIARP’2010), pages 541–548, 2010.
 *   http://dl.acm.org/citation.cfm?id=1948207.1948291
 *
 ******************************************************************************************
 *
 * This file is part of Unsupervised Distance Learning Framework (UDLF).
 *
 * UDLF is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * UDLF is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License along
 * with UDLF.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#ifndef CONTEXTRR_HPP
#define CONTEXTRR_HPP

#include <memory>
#include <map>
#include <vector>
#include <array>

#include "Udl.hpp"

class Contextrr: public Udl {
    public:
            Contextrr();

    private:
            void runUdlMethod() override;
            void runFusionMethod() override;
            void initDataStructuresUdl() override;
            void initDataStructuresFusion() override;
            void loadParameters() override;
            void checkParameters() override;
            void prepareInput() override;
            void prepareOutput() override;
            void prepareInputFusion();

            void setMatrixW(float value);
            void contextImgProc();
            void execContextImgProc();
            void execContextImgProcNoSyncW();
            void kernelCtxImageProc(int pN, int pK);
            void kernelCtxImageProcNoSyncW(int pN, int pK);
            void processMatrixWUpdates();
            void searchCurrentMaxDistance();
            void execComputeNewDists();
            void execComputeNewDistsNoMin();
            void kernelComputeNewDistancesNoMin(int posX);
            void kernelComputeNewDistances(int posX, int posY);
            void normalizeMinDistances();
            void execSortRankedLists();
            void kernelSortRankedLists(int curRL);

            //method parameters
            int l;
            int k;
            int t;
            int nByK;
            bool hasOptimizations;

            //method core values
            long int n2;
            unsigned int l2;
            unsigned int incSz;
            unsigned int nullImg;
            unsigned int maxHip;

            //global variables
            int curIteration = 0;
            float maxDist = 0;

            //method constants
            const int& cN = n;
            const int& cL = l;
            const int h = 3;
            const int h2 = 9;
            const int white = 1;
            const int black = 0;

            //method structures
            float* matrixW = NULL;
            std::vector<int> incX;
            std::vector<int> incY;
            std::vector<float> incValue;

};

#endif // CONTEXTRR_HPP
