/* <Cprr.hpp>
 *
 * CPRR method class header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *******************************************************************************************************************
 *
 * Cartesian Product of Ranking References (CPRR) is presented in the paper:
 *   L. P. Valem and D. C. G. Pedronette.
 *   "Unsupervised similarity learning through cartesian product of ranking references for image retrieval tasks."
 *   2016 29th SIBGRAPI Conference on Graphics, Patterns and Images (SIBGRAPI).
 *   http://dx.doi.org/10.1109/SIBGRAPI.2016.042
 *   http://sibgrapi.sid.inpe.br/col/sid.inpe.br/sibgrapi/2016/07.22.13.56/doc/PaperSIBGRAPI-2016_vFinal.pdf
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

#ifndef CPRR_HPP
#define CPRR_HPP

#include <memory>
#include <map>
#include <vector>
#include <array>

#include "Udl.hpp"

class Cprr: public Udl {
    public:
            Cprr();

    private:
            void runUdlMethod() override;
            void runFusionMethod() override;
            void initDataStructuresUdl() override;
            void initDataStructuresFusion() override;
            void loadParameters() override;
            void checkParameters() override;
            void prepareInput() override;
            void prepareOutput() override;

            void execFillPosMatrix();
            void execSortRankedLists();
            void execSortRankedListsZero();
            void execSortRankedListsAgg(std::vector<std::vector<int>>& rk);
            void execCartProd();
            void execReverseCartProd();
            void execCleanRef();
            void kernelFillPosMatrix(int rk);
            void kernelSortRankedLists(int rk);
            void kernelSortRankedListsZero(int rk);
            void kernelSortRankedListsAgg(int curRL, std::vector<std::vector<int>>& rk);
            void kernelCartProd(int rk);
            void kernelReverseCartProd(int rk);
            void kernelCleanRef(int rk);
            bool search(int img, const std::vector<int>& rk);

            //method parameters
            int l; //depth of the ranked lists
            int k; //number of nearest neighbors
            int t; //number of iterations

            //lists used to store references
            std::vector<unsigned int> imgRef;
            std::vector<unsigned int> posRef;
            std::vector<unsigned int> nRef;

            //aggregation structures
            float* matrixAgg = NULL;
            std::vector<std::vector<int>> rkListsAgg;

};

#endif // CPRR_HPP
