/* <RlRecom.hpp>
 *
 * RL-Recommendation method class header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 **************************************************************************************************
 *
 * RL-Recommendation is presented in the paper:
 *   L. P. Valem, D. C. G. Pedronette, R. d. S. Torres, E. Borin, and J. Almeida.
 *   "Effective, efficient, and scalable unsupervised distance learning in image retrieval tasks."
 *   ICMR, 2015
 *   http://dx.doi.org/10.1145/2671188.2749336
 *
 **************************************************************************************************
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

#ifndef RLRECOM_HPP
#define RLRECOM_HPP

#include <memory>
#include <map>
#include <vector>
#include <array>

#include "Udl.hpp"

class RlRecom: public Udl {
    public:
            RlRecom();

    private:
            void runUdlMethod() override;
            void runFusionMethod() override;
            void initDataStructuresUdl() override;
            void initDataStructuresFusion() override;
            void loadParameters() override;
            void checkParameters() override;
            void prepareInput() override;
            void prepareOutput() override;

            void execFillMatrix();
            void execFillPosMatrix();
            void execCalcCohesion();
            void execPerformRecommendations();
            void execSortRankedLists();
            void kernelFillMatrix(int rk);
            void kernelFillPosMatrix(int rk);
            void kernelCalcCohesion(int rk);
            void kernelPerformRecommendations(int rk);
            void kernelSortRankedLists(int rk);

            //auxiliary functions
            float weightNormCohesion(int position);
            float getAcumWeightCohesion(int size);
            bool isIn(int rk, int image, int size);
            float getCurrentDistance(int i, int j);
            void setNewDistance(int i, int j, float newDist);

            //method parameters
            int l; //depth of the ranked lists
            int k; //initial number of nearest neighbors
            double lambda;
            double epsilon;

            //vector to store the cohesions
            std::vector<float> cohesionVector;

};

#endif // RLRECOM_HPP
