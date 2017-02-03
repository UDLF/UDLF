/* <CorrelationGraph.hpp>
 *
 * CorrelationGraph method class header file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *************************************************************************************************
 *
 * Correlation Graph Manifold Learning is presented in the paper:
 *   D. C. G. Pedronette and R. da S. Torres.
 *   "A correlation graph approach for unsupervised manifold learning in image retrieval tasks."
 *   Neurocomputing, 208:66 – 79, 2016. SI: BridgingSemantic.
 *   http://dx.doi.org/10.1016/j.neucom.2016.03.081
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

#ifndef CORRELATION_GRAPH_HPP
#define CORRELATION_GRAPH_HPP

#include <memory>
#include <map>
#include <set>
#include <vector>

#include "Udl.hpp"

class CorrelationGraph: public Udl {
    public:
            CorrelationGraph();

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

            //general functions
            void runCorrelationGraphReRanking();
            void cleanCorrelationValues();
            void computeCorrelationForKNN();
            float rbo(int i1, int i2, int paramK);
            float pearson(int i1, int i2, int paramK);
            void cleanGraphStructures();
            void buildCorrelationGraph();
            void initializeGraph();
            void computeIncrementsByGraph();
            void buildSCC();
            void computeIncrementsBySCC();
            void incDistance(int i, int j, double inc);
            void computeFinalDistances();

            //tarjan
            void tarjan();
            void strongConnect(int v);
            void initTarjanStructures();

            //sorting (heapsort)
            void execSortRankedLists();
            void reRankingImageHeapSort(unsigned int rk);
            void heapsort(std::vector<float>& distances, std::vector<int>& curRk, int n);
            void exchange(std::vector<float>& distances, std::vector<int>& curRk, int i, int j);
            void downheap(std::vector<float>& distances, std::vector<int>& curRk, int n, int v);
            void buildheap(std::vector<float>& distances, std::vector<int>& curRk, int n);

            //rank-aggregation (fusion)
            void prepareRankAggregationMatrix();
            void normalizeProbDB();
            float* tmpMatrix = NULL;

            //graph representation
            struct Edge {
                int destination;
                float weight;
            };
            std::vector<std::vector<Edge>> graph;

            //algorithm structures
            double curThreshold;
            std::vector<std::vector<float>> correlationList;

            //tarjan structures
            int index;
            std::vector<int>  indexes;
            std::vector<int>  lowlink;
            std::vector<int>  stack;
            std::vector<bool> onStack;
            std::vector<std::vector<int>> scc; //stores tarjan's algorithm output

            //parameters
            int curK;
            double thStart;
            double thEnd;
            double thInc;
            std::string correlationFunc;
            int l;

};

#endif // CORRELATION_GRAPH_HPP
