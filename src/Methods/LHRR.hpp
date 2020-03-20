/* <HyperGraph.hpp>
 *
 * HyperGraph method class header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *******************************************************************************************************************
 *
 * Log-based Hypergraph of Ranking References (LHRR) is presented in the paper:
 *    PEDRONETTE, D. C. G.; VALEM, L. P.; ALMEIDA, J., TORRES, R. S. .
 *    Multimedia Retrieval Through Unsupervised Hypergraph-Based Manifold Ranking.
 *    IEEE Transactions on Image Processing, v. 28, p. 5824–5838, 2019.
 *    https://ieeexplore.ieee.org/document/8733193
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

#ifndef LHRR_HPP
#define LHRR_HPP

#include <memory>
#include <unordered_map>
#include <vector>
#include <array>
#include <algorithm>
#include <limits>

#include "Udl.hpp"

class LHRR: public Udl {
    public:
            LHRR();

    private:
            void runUdlMethod() override;
            void runFusionMethod() override;
            void initDataStructuresUdl() override;
            void initDataStructuresFusion() override;
            void loadParameters() override;
            void checkParameters() override;
            void prepareInput() override;
            void prepareOutput() override;

            void hyperGraphIteration();
            void initializeDataStructures();
            void loadHyperEdges();
            void loadRevHyperEdges();
            void createHyperEdge(int img);
            void includeHyperEdgeValue(int i, int j, double value);
            double weightPosition(int pos);

            void computeCartesianProductHyperEdges();
            void computeReciprocalHyperEdgesSimilarities();
            void computeDBBySimilarities();
            void computeHyperEdgesSimilarities();
            void compressHE();

            void resetDB(int value);
            int searchPairByKey(int key, std::vector<std::pair<int, double>>& hyperEdge);
            double multiplyMaps(std::vector<std::pair<int, double>>& he1, std::vector<std::pair<int, double>>& he2);

            void execFillPosMatrix();
            void execSortRankedLists();
            void kernelFillPosMatrix(int rk);
            void kernelSortRankedLists(int rk);

            void sortTmpList(int qimg);
            void joinRks(int qimg);

            void sortAll(int qimg);

            //method parameters
            int l; //depth of the ranked lists
            int k; //number of nearest neighbors
            int t; //number of iterations

            std::vector<std::vector<std::pair<int, double>>> hyperEdges;
            std::vector<std::vector<std::pair<int, double>>> revHyperEdges;
            std::vector<double> confid;

            std::vector<std::vector<int>> tmpList;

            std::vector<std::vector<int>> imgList;
            std::vector<std::vector<int>> imgListRev;
            std::vector<std::vector<float_t>> valList;
            std::vector<std::vector<float_t>> valListRev;

            float* matrixAgg = NULL;

            //time evaluation
            timeval startTimeFillPos;
            timeval startTimeInit;
            timeval startTimeLoadHE;
            timeval startTimeLoadRHE;
            timeval startTimeResetDB;
            timeval startTimeProdCart;
            timeval startTimeSimHE;
            timeval startTimeSimRHE;
            timeval startTimeSort;
            float totalTimeFillPos = 0;
            float totalTimeInit = 0;
            float totalTimeLoadHE = 0;
            float totalTimeLoadRHE = 0;
            float totalTimeResetDB = 0;
            float totalTimeProdCart = 0;
            float totalTimeSimHE = 0;
            float totalTimeSimRHE = 0;
            float totalTimeSort = 0;
};

#endif // LHRR_HPP
