/* <RlSim.hpp>
 *
 * RlSim method class header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *************************************************************************************************
 *
 * Original RL-Sim is presented in the paper:
 *   D. C. G. Pedronette and R. d. S. Torres.
 *   "Image re-ranking and rank aggregation based on similarity of ranked lists."
 *   Pattern Recognition, 46(8):2350–2360, 2013.
 *   http://dx.doi.org/10.1016/j.patcog.2013.01.004
 *
 * RL-Sim* is presented in the paper:
 *   C. Y. Okada, D. C. G. a. Pedronette, and R. da S. Torres.
 *   "Unsupervised distance learning by rank correlation measures for image retrieval."
 *   Proceedings of the 5th ACM on International Conference on Multimedia Retrieval, ICMR ’15,
 *   pages 331–338, New York, NY, USA, 2015. ACM.
 *   http://dx.doi.org/10.1145/2671188.2749335
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

#ifndef RLSIM_HPP
#define RLSIM_HPP

#include <memory>
#include <map>
#include <vector>
#include <array>

#include "Udl.hpp"

class RlSim: public Udl {
    public:
            RlSim();

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

            void execUpdateDistances();
            void execSortRankedLists();
            void kernelUpdateDistances(int i1);
            void kernelSortRankedLists(int rk);

            void resetTmpMatrix();
            void releaseTmpMatrix();
            void normalizeMinDistances();

            void prepareRankAggregationMatrix();

            //metrics
            float intersection(int i1, int i2);
            float rbo(int i1, int i2);
            float kendallTau(int i1, int i2);
            float kendallTauW(int i1, int i2);
            float jaccard(int i1, int i2);
            float jaccardK(int i1, int i2);
            float goodman(int i1, int i2);
            float spearman(int i1, int i2);

            //method parameters
            int topK;
            int cK;
            int t;
            std::string metric;

            //auxiliar structures
            float* tmpMatrix = NULL;

};

#endif // RLSIM_HPP
