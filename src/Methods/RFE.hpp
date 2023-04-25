/* <RFE.hpp>
 *
 * RFE header file
 *
 ***********************************************************************************
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

#ifndef RFE_HPP
#define RFE_HPP

#include "Udl.hpp"
#include <algorithm>
#include <numeric>
#include <set>

#include <map>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

class RFE: public Udl {
    public:
            RFE();

    private:
            void runUdlMethod() override;
            void runFusionMethod() override;
            void initDataStructuresUdl() override;
            void initDataStructuresFusion() override;
            void loadParameters() override;
            void checkParameters() override;
            void prepareInput() override;
            void prepareOutput() override;

            void rankNorm();
            float rankScoreSigmoid(int x, int pk);

            void execSortRankedLists();
            void kernelSortRankedLists(int rk);

            void execEmbeddingRR();
            void doHyperRREmbedding();
            void createEmb();
            void createReverseEmb();
            void combineReverseEmb();
            float rankKScore(int pos);
            void reRankingByEmbSim();

            void execCartesianProduct();
            void doHyperEffecEstimation();
            void reRankingByCPRR();

            // Effectiveness Estimations
            void sortEmb();
            void effectEstimationEmb();

            // Utils Embeddings: Inc e Sim
            float simEmb(const std::map<int, float>& di, const std::map<int, float>& dj);
            void incEmb(std::map<int, float>& d, int img, float value);

            void incSimByECP(int qi, int qj, float wi, float wj);

            void kernelSortRankedLists_extra(int rk, std::set<int> extra_elems);

            void execRRCCs();
            void doCCsEmbedding();
            void createEdgesCCS();
            void processCCs();
            void initEdgesCCs();
            void computeCCsEmb();
            void computeEmbsByCCs();
            void reRankingByCCs();

            void execRRByEmbCCs();

            void unionCCs(int o1, int o2);
            bool sameCC(int o1, int o2);

            unsigned int getPosition(unsigned int imgId1, unsigned int imgId2);

            void sortCCEmb();

            void exportEmbeddingsFile();
            void exportCCsFile();
            std::vector<int> exportCCs();

            void execSortRankedListsAgg(std::vector<std::vector<int>>& rk);
            void kernelSortRankedListsAgg(int curRL, std::vector<std::vector<int>>& rk);
            bool search(int img, const std::vector<int>& rk);

            void acumEmb();

            void reRankingByEmbSimAgg();

            float cosineSimCCs(const std::map<int, float>& di, const std::map<int, float>& dj);

            void reRankingByEmbCC();

            std::map<int, std::map<int, float>> emb;   // embeddings
            std::vector<std::vector<int>> sorted_emb;    // sorted embeddings indexes
            std::map<int, std::map<int, float>> remb;  // reverse embeddings
            std::vector<float> ee; // effectiveness estimations
            std::vector<int> sorted_ee;  // sorted effectiveness estimations indexes

            std::vector<std::set<int>> incElems;
            std::vector<std::set<int>> incElemsL;

            std::vector<std::pair<float, std::pair<int, int>>> edgesCCs;

            std::map<int, std::set<int>> rcc;
            std::map<int, int> icc;

            std::map<int, std::map<int, float>> ccemb;
            std::vector<std::vector<int>> sorted_ccemb;    // sorted cc_embeddings indexes

            std::map<int, std::map<int, float>> embycc;

            //RFE Parameters (read from the configuration file)
            int l;
            int k;
            int t;
            double pa;
            double th_cc;
            bool performCCs;
            bool rerankByEmb;
            bool exportEmbeddings;
            std::string embeddingsPath;
            std::string ccsPath;

            //aggregation structures
            float* matrixAgg = NULL;
            std::vector<std::vector<int>> rkListsAgg;

};

#endif // RFE_HPP
