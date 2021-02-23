/* <RDPAC.hpp>
 *
 * RDPAC header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *******************************************************************************************************************
 *
 * Rank-based Diffusion Process with Assured Convergence (RDPAC) is presented in the paper:
 *    PEDRONETTE, D. C. G.; VALEM, L. P.; LATECKI, L. J. .
 *    Efficient Rank-based Diffusion Process with Assured Convergence.
 *    Journal of Imaging (ISSN 2313-433X).
 *    (Under minor revision)
 *
 *******************************************************************************************************************
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

#ifndef RDPAC_HPP
#define RDPAC_HPP

#include "Udl.hpp"

class RDPAC: public Udl {
    public:
            RDPAC();

    private:
            void runUdlMethod() override;
            void runFusionMethod() override;
            void initDataStructuresUdl() override;
            void initDataStructuresFusion() override;
            void loadParameters() override;
            void checkParameters() override;
            void prepareInput() override;
            void prepareOutput() override;

            //Copy two matrices
            void copyMatrices(float* src, float* dst, long int lim);

            //Easily get the element from rk
            inline int getRKElem(long int query, long int pos);
            inline float getMatrixElem(float* m, long int i, long int j);
            inline void setMatrixElem(float* m, long int i, long int j, float value);

            //Sorting method
            void execSortRankedLists();
            void kernelSortRankedLists(int rk);

            //RDPAC Methods
            void performSingleExecution();
            void execFillPosMatrix();
            void kernelFillPosMatrix(int rk);
            void reciprocalReferences();
            void generateFullLMatrixW();
            void sumTranspostL();
            void copyingMatrixWToDBL();
            void cleanSimMatrix();
            void cleanRWLMatrix();
            void generatekNNMatrixW(int k);
            void normalizeMatrix(float* matrix, int lim);
            void multiplyByTranspostL(int k);
            void multiplyByTranspostLStoreList(int k);
            void sumIdentityAlpha(float alpha);
            void multiplyL1(float* m1, float* m2, float* dst);
            void multiplyL2(float* m1, float* m2, float* dst);
            //Methods for fusion
            void copyMatrixWtoFuseWL();
            void copyMatrixFuseWtoDBL();

            //RDPAC Structures
            float* W = NULL;
            float* RWL = NULL; // also known as matrix P
            float* tmpLine = NULL;
            std::map<int, std::vector<int>> nonZero;
            //Fusion
            float* fuseW = NULL;
            std::map<int, std::vector<int>> setnz;

            //RDPAC Parameters (read from the configuration file)
            int l;
            double l_mult; // value to multiply L
            double p;
            double pl;
            int k_start;
            int k_inc;
            int k_end;
            bool k_auto;

};

#endif // LATECKI_HPP
