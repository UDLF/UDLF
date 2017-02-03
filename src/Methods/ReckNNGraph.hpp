/* <ReckNNGraph.hpp>
 *
 * ReckNNGraph method class header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *****************************************************************************************************************
 *
 * ReckNNGraph is presented in the paper:
 *   D. C. G. Pedronette, O. A. Penatti, and R. da S. Torres.
 *   "Unsupervised manifold learning using Reciprocal kNN Graphs in image re-ranking and rank aggregation tasks."
 *   Image and Vision Computing, 32(2):120 – 130, 2014.
 *   http://dx.doi.org/10.1016/j.imavis.2013.12.009
 *
 *****************************************************************************************************************
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

#ifndef RECKNNGRAPH_HPP
#define RECKNNGRAPH_HPP

#include <memory>
#include <map>
#include <vector>
#include <array>

#include "Udl.hpp"

class ReckNNGraph: public Udl {
    public:
            ReckNNGraph();

    private:
            void runUdlMethod() override;
            void runFusionMethod() override;
            void initDataStructuresUdl() override;
            void initDataStructuresFusion() override;
            void loadParameters() override;
            void checkParameters() override;
            void prepareInput() override;
            void prepareOutput() override;

            //main functions
            void initializeVariables();
            void initializeMatrix();
            void initializePositions();
            void initializeDB();
            void iterationReciprocalkNNGraph();
            double evaluateRankedList(int imgId, float scoreThreshold);
            double evaluateSetByCliqueNonWeithted(const std::vector<unsigned int>& set, int size);
            void computeIncrements(double score, const std::vector<unsigned int>& set, double scoreThreshold, int size);
            void computeNewDB();
            void updateModels();
            void sortRankedList(int curRL);

            //aggregation (fusion) functions
            void runRankAggregation();
            void copyDB();
            void aggregate();
            float checkAggregateDistance(float curDist);

            //auxiliary functions
            unsigned int getImgAtPos(unsigned int qImg, unsigned int iPos);
            float getDistance(unsigned int imgId1, unsigned int imgId2);
            void setNewDistance(unsigned int imgId1, unsigned int imgId2, float newDist);
            float getTempDistance(unsigned int imgId1, unsigned int imgId2);
            void setTempNewDistance(unsigned int imgId1, unsigned int imgId2, float newDist);
            unsigned int getPosition(unsigned int imgId1, unsigned int imgId2);

            //method parameters
            int l;
            int k;
            double epsilon;

            //method core/global variables
            int knn;
            int curIterationReckNN;
            double convScore, prevConvScore, diffConvScore;

            //method structures
            float* tmpMatrix = NULL;
            unsigned int* posMatrix = NULL;

};

#endif // RECKNNGRAPH_HPP
