/* <RkGraph.hpp>
 *
 * RkGraph method header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *************************************************************************************************
 *
 * Ranked List Graph Distance is presented in the paper:
 *   D. C. G. Pedronette, J. Almeida, and R. da S. Torres.
 *   "A graph-based ranked-list model for unsupervised distance learning on shape retrieval."
 *   Pattern Recognition Letters, 83, Part 3:357 – 367, 2016. Efficient Shape Representation,
 *   Matching, Ranking, and its Applications.
 *   http://dx.doi.org/10.1016/j.patrec.2016.05.021
 *
 *************************************************************************************************
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

#ifndef RKGRAPH_HPP
#define RKGRAPH_HPP

#include "Udl.hpp"

class RkGraph: public Udl {
    public:
            RkGraph();

    private:
            void runUdlMethod() override;
            void runFusionMethod() override;
            void initDataStructuresUdl() override;
            void initDataStructuresFusion() override;
            void loadParameters() override;
            void checkParameters() override;
            void prepareInput() override;
            void prepareOutput() override;

            void initGraphStructures();
            void initStructure(double*& structure);

            void computeMutualRankDists();
            void runIteration();
            void fillSimMetricForImage(int img);
            void setSimMetric(int img1, int img2);
            double calcRBO(int i1, int i2);
            void computeIncAdjForImage(int imgq);
            void computeNewDists();

            void execSortRankedLists();
            void kernelSortRankedLists(int rk);

            int getPosition(int qImg, int img);

            //method parameters
            int k;
            int t;
            int l;
            double p;

            //auxiliar structures
            double* smetric = NULL;
            double* adj = NULL;

            //aggregation matrix
            float* matrixAgg = NULL;

};

#endif // RKGRAPH_HPP
