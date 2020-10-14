/* <BFSTree.hpp>
 *
 * BFSTree header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *******************************************************************************************************************
 *
 * The BFSTree algorithm is presented in the paper:
 *    PEDRONETTE, D. C. G.; VALEM, L. P.; TORRES, R. S. .
 *    "A BFS-Tree of ranking references for unsupervised manifold learning."
 *    Pattern Recognition, v. 111, p. 107666, 2021.
 *    https://doi.org/10.1016/j.patcog.2020.107666
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

#ifndef BFSTree_HPP
#define BFSTree_HPP

#include "Udl.hpp"

class BFSTree: public Udl {
    public:
            BFSTree();

    private:
            void runUdlMethod() override;
            void runFusionMethod() override;
            void initDataStructuresUdl() override;
            void initDataStructuresFusion() override;
            void loadParameters() override;
            void checkParameters() override;
            void prepareInput() override;
            void prepareOutput() override;

            // Get and set elem
            inline int getRKElem(long int query, long int pos);
            inline float getMatrixElem(float* m, long int i, long int j);
            inline void setMatrixElem(float* m, long int i, long int j, float value);
            inline float minZero(float x, float y);
            inline void incSimW(int img1, int img2, float newSim);

            void normalizeWProb();
            void computeWDiffusion();

            // AcumTree
            void acumBFSTree();
            void computeIncTree(int qImg);

            // Normalization
            void computeRankNormalization(int type);
            void computeRankNormalization2(int type);

            // Sorting
            void execSortRankedLists();
            void kernelSortRankedLists(int rk);
            void execSortRankedLists2();
            void kernelSortRankedLists2(int rk);

            // Correlation Measures
            void  computeRankCorrelation();
            float rankCorrelationMeasure(int i1, int i2);
            float intersection(int i1, int i2);
            float rbo(int i1, int i2);
            float kendallTau(int i1, int i2);
            float kendallTauW(int i1, int i2);
            float jaccard(int i1, int i2);
            float jaccardK(int i1, int i2);
            float spearman(int i1, int i2);
            float goodman(int i1, int i2);

            // Auxiliar Structures
            std::vector<std::map<int, float>> imgCor;

            // BFSTree Parameters (read from the configuration file)
            int l;
            int k;
            int topK;
            std::string correlation_metric;

};

#endif // BFSTree_HPP
