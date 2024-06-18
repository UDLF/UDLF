/* <RFE.cpp>
 *
 * RFE implementation file
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

#include <iostream>
#include <math.h>
#include <omp.h>

#include "RFE.hpp"

/* Constructor */
RFE::RFE() {

}

void RFE::loadParameters() {
    Exec exec = Exec::getInstance();

    //RFE Parameters (read from the configuration file)
    exec.getConfigVariable(k, "PARAM_RFE_K");
    exec.getConfigVariable(t, "PARAM_RFE_T");
    exec.getConfigVariable(l, "PARAM_RFE_L");
    exec.getConfigVariable(pa, "PARAM_RFE_PA");
    exec.getConfigVariable(th_cc, "PARAM_RFE_TH_CC");
    exec.getConfigVariable(performCCs, "PARAM_RFE_PERFORM_CCS");
    exec.getConfigVariable(rerankByEmb, "PARAM_RFE_RERANK_BY_EMB");
    exec.getConfigVariable(exportEmbeddings, "PARAM_RFE_EXPORT_EMBEDDINGS");
    exec.getConfigVariable(embeddingsPath, "PARAM_RFE_EMBEDDINGS_PATH");
    exec.getConfigVariable(ccsPath, "PARAM_RFE_CCS_PATH");
}

void RFE::checkParameters() {
    if (l > n) {
        std::cout << "L can't be greater than N" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
}

void RFE::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*l);
    initSparseMatrix(matrix);

    std::cout << "Initialized successfully!" << std::endl;
}

void RFE::initDataStructuresFusion() {
    std::cout << "Initializing data structures..." << std::endl;

    initDataStructuresUdl(); //init the main structures

    //structures used for aggregation (fusion) tasks
    initMatrix(matrixAgg);
    for (long int q = 0; q < ((long int) n)*n; q++) {
        matrixAgg[q] = 0;
    }
    rkListsAgg.resize(n);  //make sure it's empty

    std::cout << "Initialized successfully!" << std::endl;
}

void RFE::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
    }
}

void RFE::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void RFE::runUdlMethod() {
    std::cout << "\n Executing RFE!\n\n";

    rankNorm();

    for (int ct = 0; ct < t; ct++) {
        std::cout << "\n\t [+] Executing Embedding RR, iteration t=" << ct << "\n";
        execEmbeddingRR();
    }

    // Cartesian Product of Hyperedges
    execCartesianProduct();

    // Connected Components Re-Ranking
    if (performCCs) {
        execRRCCs();
    }

    // Export Embeddings
    if (rerankByEmb || exportEmbeddings) {
        execRRByEmbCCs();
    }
}

void RFE::runFusionMethod() {
    for (std::string const& file : fusionFiles) {
        //prepare the ranked lists to run
        gettimeofday(&startTimeToDecrement, NULL);
            readInputFile(file);
            prepareInput();
        totalTimeToDecrement = Time::addTime(startTimeToDecrement, totalTimeToDecrement);

        //all the images that appear in the main ranked list are stored in a secondary one (union of ranked lists)
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                int img = rkLists[l*i + j];
                if (!search(img, rkListsAgg[i])) {
                    rkListsAgg[i].push_back(img);
                }
            }
        }

        //execute the main method
        rankNorm();

        //store values in a secondary matrix
        for (long int i = 0; i < n; i++) {
            for (int j = 0; j < l; j++) {
                int img = rkLists[l*i + j];
                matrixAgg[((long int) n)*i + img] += matrix[((long int) n)*i + img] + 1;
            }
        }

        //reinitialize the main matrix for the next iteration
        initSparseMatrix(matrix);
    }

    std::cout << std::endl << " Fusion Steps: " << std::endl;

    delete [] matrix;   //the main values are not necessary anymore
    matrix = matrixAgg; //the current matrix is now the one that has the accumulated values

    //execute a sorting to prepare for the final iterations
    execSortRankedListsAgg(rkListsAgg);

    //copy the content between ranked lists
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            rkLists[l*i + j] = rkListsAgg[i][j];
        }
    }
    rkListsAgg.clear(); //these values are not necessary anymore

    //final processing
    std::cout << std::endl << " Final Processing: " << std::endl;
    runUdlMethod();
}

void RFE::rankNorm() {
    std::cout << "\t [+] Normalizing Rankings ...\n";
    #pragma omp parallel for
    for (int q = 0; q < n; q++) {
        for (int i = 0; i < l; i++) {
            int img_rl_i = rkLists[l*q + i];
            int posqi = i+1;
            int posiq = getPosition(img_rl_i, q);
            float newSim = (pow(rankScoreSigmoid(posqi,k/2), 2) * rankScoreSigmoid(posiq,k/2)); // better
            matrix[((long int) n)*q + img_rl_i] = newSim;
            //matrix[n*img_rl_i + q] = newSim;
        }
    }
    execSortRankedLists();
}

float RFE::rankScoreSigmoid(int x, int pk) {
    if (x==1)
        return 1;

    float sigmoidv =  1/ (1+exp((-pa) * (x-pk) ));
    return 1-sigmoidv;
}

void RFE::execEmbeddingRR() {
    doHyperRREmbedding();
    reRankingByEmbSim();
}

void RFE::doHyperRREmbedding() {
    createEmb();
    createReverseEmb();
    combineReverseEmb();
}

void RFE::execCartesianProduct() {
    doHyperRREmbedding();
    doHyperEffecEstimation();
    reRankingByCPRR();
}

void RFE::doHyperEffecEstimation() {
    sortEmb();
    effectEstimationEmb();
}

void RFE::execRRCCs() {
    doHyperRREmbedding();
    doHyperEffecEstimation();
    doCCsEmbedding();
    reRankingByCCs();
}

void RFE::doCCsEmbedding() {
    createEdgesCCS();
    processCCs();
    computeCCsEmb();
    computeEmbsByCCs();
}

void RFE::execRRByEmbCCs() {
    doHyperRREmbedding();
    doHyperEffecEstimation();
    doCCsEmbedding();
    if (rerankByEmb) {
        reRankingByEmbCC();
    }
    if (exportEmbeddings) {
        exportEmbeddingsFile();
        exportCCsFile();
    }
}

void RFE::reRankingByEmbCC() {
    std::cout << "\n\t  [+] Re-Ranking by Embeddings CC \n";
    initSparseMatrix(matrix);
    for (int q = 0; q < n; q++) {
        for (int pos = 0; pos < l; pos++) {
            int img = rkLists[q*l + pos];
            float newSim = (1+cosineSimCCs(embycc[q], embycc[img])) / (pos+1);
            matrix[((long int) n)*q + img] = newSim;
        }
    }
    execSortRankedLists();
}

void RFE::exportEmbeddingsFile() {
    std::cout << "\n\t [+] Exporting embeddings to " << embeddingsPath << "\n";
    std::ofstream fileOutStream;
    fileOutStream.open(embeddingsPath);
    for (int q = 0; q < n; q++) {
        for (auto it : embycc[q]) {
            fileOutStream << it.second << " ";
        }
        fileOutStream << "\n";
    }
    fileOutStream.close();
    std::cout << "\t\t Done exporting!\n";
}

void RFE::exportCCsFile() {
    std::cout << "\n\t [+] Exporting CCs to " << ccsPath << "\n";
    std::vector<int> ttsets = exportCCs();
    std::ofstream fileOutStream;
    fileOutStream.open(ccsPath);
    for (int q = 0; q < n; q++) {
        fileOutStream << ttsets[q] << "\n";
    }
    fileOutStream.close();
    std::cout << "\t\t Done exporting!\n";
}

std::vector<int> RFE::exportCCs() {
    int c = 0;
    std::vector<int> ttsets;
    ttsets.resize(n);
    for (int q = 0; q < n; q++) {
        ttsets[q] = -1;
    }
    std::map<int, int> dcc;
    for (int q = 0; q < n; q++) {
        int curRep = icc[q];
        if (rcc[curRep].size() > 1) {
            dcc[curRep] = 1;
            int curIdx = distance(dcc.begin(), dcc.find(curRep));
            ttsets[c] = curIdx;
        }
        c++;
    }
    return ttsets;
}

void RFE::createEdgesCCS() {
    std::cout << "\t\t [*] Creating Edges Candidates\n";

    for (int q = 0; q < n; q++) {
        for (int i = 0; i < k; i++) {
            int img_nn = rkLists[l*q + i];
            if (q == img_nn) {  // skip query
                continue;
            }
            float confValue = (simEmb(emb[q], emb[img_nn]) * ee[q] * ee[img_nn]);
            std::pair<int, int> img_pair = {q, img_nn};
            std::pair<float, std::pair<int, int>> tup = {confValue, img_pair};
            edgesCCs.push_back(tup);
        }
    }
    // Using simple sort() function to sort
    std::sort(edgesCCs.rbegin(), edgesCCs.rend());
}

void RFE::processCCs() {
    std::cout << "\t\t [*] Creating CCs\n";
    initEdgesCCs();
    int i = 0;
    for (auto t : edgesCCs) {
        if (i >= ((int) edgesCCs.size()*th_cc)) {
            break;
        }
        std::pair<int, int> img_pair = t.second;
        if (!(sameCC(img_pair.first, img_pair.second))) {
            unionCCs(img_pair.first, img_pair.second);
        }
        i++;
    }
}

void RFE::initEdgesCCs() {
    // CC Threshold
    if (th_cc == 0) {
        float avg = 0;
        for (int q : sorted_ee) {
            avg += ee[q];
        }
        avg = avg / n;
        th_cc = avg / 2;
    }
    // CCs Structures
    rcc.clear();
    icc.clear();
    for (int q = 0; q < n; q++) {
        rcc[q] = std::set<int>{q};
        icc[q] = q;
    }
}

void RFE::unionCCs(int o1, int o2) {
    int r1 = icc[o1];
    int r2 = icc[o2];

    std::set_union(rcc[r2].cbegin(), rcc[r2].cend(), rcc[r1].cbegin(), rcc[r1].cend(), std::inserter(rcc[r1], rcc[r1].cend()));
    for (auto o : rcc[r2]) {
        icc[o] = r1;
    }

    rcc.erase(r2);
}

bool RFE::sameCC(int o1, int o2) {
    return (icc[o1] == icc[o2]);
}

void RFE::computeCCsEmb() {
    ccemb.clear();
    for (auto it1 : rcc) {
        int r = it1.first;
        std::set<int> d = it1.second;
        if (d.size() > 1) {
            std::map<int, float> curemb;
            for (auto o : rcc[r]) {
                for (auto it2 : emb[o]) {
                    int key = it2.first;
                    curemb[key] = curemb[key] + emb[o][key]; //original
                }
            }
            ccemb[r] = curemb;
        }
    }
    sortCCEmb();
}

void RFE::computeEmbsByCCs() {
    embycc.clear();
    for (int q = 0; q < n; q++) {
        auto qemb = emb[q];
        std::map<int, float> curCCEmb;
        for (auto it : ccemb) {
            int r = it.first;
            auto ccd = it.second;
            float vccemb = simEmb(qemb, ccd); // original
            curCCEmb[r] = vccemb;
        }
        embycc[q] = curCCEmb;
    }
}

void RFE::reRankingByCCs() {
    incElems.clear();
    incElemsL.clear();
    incElems.resize(n);
    incElemsL.resize(n);

    std::cout << "\n\n\t [+] Re-Ranking by CCs \n";
    initSparseMatrix(matrix);
    for (auto it : ccemb) {
        int r = it.first;
        std::map<int, float> d = it.second;
        for (int i = 0; i < std::min(k, ((int) d.size())); i++) {
            int k1 = sorted_ccemb[r][i];
            for (int j = 0; j < std::min(k, ((int) d.size())); j++) {
                int k2 = sorted_ccemb[r][j];
                float newSim = (1+sqrt((i+1)*(j+1)) * simEmb (embycc[k1], embycc[k2])) / getPosition(k1, k2);
                matrix[((long int) n)*k1 + k2] = newSim;
                incElems[k1].insert(k2);
            }
        }
    }

    for (int q = 0; q < n; q++) {
        for (int i = 0; i < l; i++) {
            int img_i = rkLists[l*q + i];
            incElemsL[q].insert(img_i);
        }
    }

    std::cout << "\n\t\t [*] Sort Ranked Lists ... \n";
    for (int q = 0; q < n; q++) {
        std::set<int> result;
        std::set_difference(incElems[q].cbegin(), incElems[q].cend(), incElemsL[q].cbegin(), incElemsL[q].cend(), std::inserter(result, result.cend()));
        kernelSortRankedLists_extra(q, result);
        //reRankingImageHeapSort(q, result);
    }
}

unsigned int RFE::getPosition(unsigned int qImg, unsigned int img) {
    for (int i = 0; i < l; i++) {
        if (rkLists[l*qImg + i] == img) {
            return (i+1);
        }
    }
    return (l+2);
}

// Effectiveness Estimation

// Comparator function to sort pairs according to second value
bool cmp(std::pair<int, float>& a, std::pair<int, float>& b) {
    return a.second > b.second;
}

void RFE::sortEmb() {
    sorted_emb.clear();
    sorted_emb.resize(n);
    #pragma omp parallel for
    for (int q = 0; q < n; q++) {
        // Declare vector of pairs
        std::vector<std::pair<int, float>> A;

        // Copy key-value pair from Map to vector of pairs
        for (auto& it : emb[q]) {
            A.push_back(it);
        }

        // Sort using comparator function
        std::sort(A.begin(), A.end(), cmp);

        // Add index to vector
        for (auto& it : A) {
            sorted_emb[q].push_back(it.first);
        }
    }
}

void RFE::sortCCEmb() {
    sorted_ccemb.resize(n);
    for (auto iterator : ccemb) {
        int q = iterator.first;
        // Declare vector of pairs
        std::vector<std::pair<int, float>> A;

        // Copy key-value pair from Map to vector of pairs
        for (auto& it : ccemb[q]) {
            A.push_back(it);
        }

        // Sort using comparator function
        std::sort(A.begin(), A.end(), cmp);

        // Add index to vector
        for (auto& it : A) {
            sorted_ccemb[q].push_back(it.first);
        }
    }
}

template<typename T>
std::vector<int> tag_sort(const std::vector<T>& v)
{
    std::vector<int> result(v.size());
    std::iota(std::begin(result), std::end(result), 0);
    std::sort(std::begin(result), std::end(result),
            [&v](const auto & lhs, const auto & rhs)
            {
                return v[lhs] > v[rhs];
            }
    );
    return result;
}

void RFE::effectEstimationEmb() {
    ee.clear();
    ee.resize(n);
    #pragma omp parallel for
    for (int q = 0; q < n; q++) {
        auto qemb = emb[q];
        float eeValue = 0;
        int pos = 0;
        for (int& key : sorted_emb[q]) {
            pos++;
            if (pos > k) {
                break;
            }
            eeValue += emb[q][key];
        }
        ee[q] = eeValue;
    }
    sorted_ee = tag_sort(ee);
    // Normalization [0,1]
    int key = sorted_ee[0];
    float grValue = ee[key];
    for (int q = 0; q < n; q++) {
        ee[q] = ee[q] / grValue;
    }
}

void RFE::reRankingByCPRR() {
    std::cout << "\n\t [+] Re-Ranking by CPRR\n";

    initSparseMatrix(matrix);

    incElems.resize(n);
    incElemsL.resize(n);
    for (int q : sorted_ee) {
        float value = ee[q];
        incSimByECP(q, q, value, 1);
    }

    #pragma omp parallel for
    for (int q = 0; q < n; q++) {
        int pos = 1;
        for (int i = 0; i < l; i++) {
            int img_i = rkLists[l*q + i];
            pos++;
            incElemsL[q].insert(img_i);
        }
        matrix[((long int) n)*q + q] = pow(2, 30);
    }

    std::cout << "\n\t\t [*] Sort Ranked Lists ... \n";
    #pragma omp parallel for
    for (int q = 0; q < n; q++) {
        std::set<int> result;
        std::set_difference(incElems[q].cbegin(), incElems[q].cend(), incElemsL[q].cbegin(), incElemsL[q].cend(), std::inserter(result, result.cend()));
        kernelSortRankedLists_extra(q, result);
    }
}

void RFE::createEmb() {
    std::cout << "\t\t [*] Creating Embeddings ... \n";
    emb.clear();
    #pragma omp parallel for
    for (int q = 0; q < n; q++) {
        std::map<int, float> qemb;
        for (int pos1 = 0; pos1 < k; pos1++) {
            int img_n = rkLists[l*q + pos1];
            float wn = rankKScore(pos1+1);
            incEmb(qemb, img_n, wn);
            for (int pos2 = 0; pos2 < k; pos2++) {
                float wnn = rankKScore(pos2+1);
                int img_nn = rkLists[l*img_n + pos2];
                incEmb(qemb, img_nn, wn*wnn);
            }
        }
        emb[q] = qemb;
    }
}

void RFE::acumEmb() {
    std::cout << "\t\t [*] Creating Embeddings ... \n";
    for (int q = 0; q < n; q++) {
        for (int pos1 = 0; pos1 < k; pos1++) {
            int img_n = rkLists[l*q + pos1];
            float wn = rankKScore(pos1+1);
            incEmb(emb[q], img_n, wn);
            for (int pos2 = 0; pos2 < k; pos2++) {
                float wnn = rankKScore(pos2+1);
                int img_nn = rkLists[l*img_n + pos2];
                incEmb(emb[q], img_nn, wn*wnn);
            }
        }
    }
}

void RFE::createReverseEmb() {
    std::cout << "\t\t [*] Creating Reverse Embeddings ... \n";
    remb.clear();
    for (int q = 0; q < n; q++) {
        std::map<int, float> qremb;
        remb[q] = qremb;
    }
    for (int q = 0; q < n; q++) {
        std::map<int, float> qemb = emb[q];
        for (auto it = qemb.cbegin(); it != qemb.cend(); it++) {
            auto key = it->first;
            incEmb(remb[key], q, qemb[key]);
        }
    }
}

void RFE::combineReverseEmb() {
    #pragma omp parallel for
    for (int q = 0; q < n; q++) {
            std::map<int, float> qemb = emb[q];
            for (auto it = qemb.cbegin(); it != qemb.cend(); it++) {
                auto key = it->first;
                emb[q][key] = (emb[q][key]) * (remb[q][key]);
            }
    }
}

float RFE::rankKScore(int pos) {
    return (1 - (log(pos)/log(k)) );
}

void RFE::reRankingByEmbSim() {
    std::cout << "\t\t [*] Computing Similarities between Embeddings ... \n";
    initSparseMatrix(matrix);
    #pragma omp parallel for
    for (int q = 0; q < n; q++) {
        int pos = 1;
        for (int i = 0; i < l; i++) {
            int img_i = rkLists[l*q + i];
            float newSim  = (simEmb(emb[q], emb[img_i])) / pos;
            matrix[((long int) n)*q + img_i] = newSim;
            pos += 1;
        }
    }
    execSortRankedLists();
}

void RFE::reRankingByEmbSimAgg() {
    std::cout << "\t\t [*] Computing Similarities between Embeddings (Agg) ... \n";
    initSparseMatrix(matrix);
    #pragma omp parallel for
    for (int q = 0; q < n; q++) {
        for (int i = 0; i < rkListsAgg[q].size(); i++) {
            int img_i = rkListsAgg[q][i];
            float newSim  = (simEmb(emb[q], emb[img_i]));
            matrix[((long int) n)*q + img_i] = newSim;
        }
    }
    execSortRankedListsAgg(rkListsAgg);
}

// Utils Embeddings: Inc e Sim
float RFE::simEmb(const std::map<int, float>& di, const std::map<int, float>& dj) {
    float sim = 0;
    for (auto& pair : di) {
        auto dj_iter = dj.find(pair.first);
        if (dj_iter != dj.end()) {
            sim += pair.second * dj_iter->second;
        }
    }
    return sim;
}

void RFE::incEmb(std::map<int, float>& d, int img, float value) {
    if (d.find(img) != d.end()) {
        d[img] += value;
    } else {
        d[img] = value;
    }
}

void RFE::incSimByECP(int qi, int qj, float wi, float wj) {
    const auto& iemb = emb[qi];
    const auto& jemb = emb[qj];

    #pragma omp parallel for
    for (const int i : sorted_emb[qi]) {
        auto it_i = iemb.find(i);
        if (it_i == iemb.end()) continue;  // Ensure the key exists in the map
        float vi = it_i->second;

        for (const int j : sorted_emb[qj]) {
            auto it_j = jemb.find(j);
            if (it_j == jemb.end()) continue;  // Ensure the key exists in the map
            float vj = it_j->second;

            incElems[i].insert(j);
            long int idx = static_cast<long int>(n) * i + j;
            matrix[idx] += (vi * vj) * wi * wj;
        }
    }
}

void RFE::execSortRankedLists() {
    std::cout << "\n\t\t [*] Sort Ranked Lists ... \n";
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists(i);
    }
}

void RFE::kernelSortRankedLists(int rk) {
    long int LcurRL = l*rk;
    long int cNcurRL = ((long int) n)*rk;
    float a[l];

    long int index;
    for (int j = 0; j < l; j++) {
        index = cNcurRL + rkLists[LcurRL + j];
        a[j] = matrix[index];
    }

    //---------------------- INSERTION SORT --------------------------
    int i, j, keyR;
    float keyA;

    for (j = 1; j < l; j++) {
        keyA = a[j];
        keyR = rkLists[LcurRL + j];
        i = j - 1;
        while (i >= 0 && a[i] < keyA) {
            a[i + 1] = a[i];
            rkLists[LcurRL + i + 1] = rkLists[LcurRL + i];
            i--;
        }
        a[i + 1] = keyA;
        rkLists[LcurRL + i + 1] = keyR;
    }
    //----------------------------------------------------------------

    //Setting query image at first position
    i = 0;
    while ((rkLists[LcurRL + i] != rk)&&(i < l)) {
        i++;
    }
    if (i > 0) {
        int aux = rkLists[LcurRL + 0];
        rkLists[LcurRL + 0] = rkLists[LcurRL + i];
        rkLists[LcurRL + i] = aux;

        float auxA = a[0];
        a[0] = a[i];
        a[i] = auxA;
    }
}

void RFE::kernelSortRankedLists_extra(int rk, std::set<int> extra_elems) {
    long int cNcurRL = ((long int) n)*rk;
    float a[l+extra_elems.size()];

    std::vector<int> rkTmp;
    rkTmp.resize(l+extra_elems.size());

    int pos = 0;
    for (int i =0; i < l; i++) {
        rkTmp[i] = rkLists[l*rk + i];
        pos++;
    }
    for (int img : extra_elems) {
        rkTmp[pos] = img;
        pos++;
    }

    int new_size = rkTmp.size();


    long int index;
    for (int j = 0; j < new_size; j++) {
        index = cNcurRL + rkTmp[j];
        a[j] = matrix[index];
    }


    //---------------------- INSERTION SORT --------------------------
    int i, j, keyR;
    float keyA;

    for (j = 1; j < new_size; j++) {
        keyA = a[j];
        keyR = rkTmp[j];
        i = j - 1;
        while (i >= 0 && a[i] < keyA) {
            a[i + 1] = a[i];
            rkTmp[i + 1] = rkTmp[i];
            i--;
        }
        a[i + 1] = keyA;
        rkTmp[i + 1] = keyR;
    }
    //----------------------------------------------------------------

    //Setting query image at first position
    i = 0;
    while ((rkTmp[i] != rk)&&(i < new_size)) {
        i++;
    }
    if (i > 0) {
        int aux = rkTmp[0];
        rkTmp[0] = rkTmp[i];
        rkTmp[i] = aux;

        float auxA = a[0];
        a[0] = a[i];
        a[i] = auxA;
    }

    // copy back
    for (int i = 0; i < l; i++) {
        rkLists[l*rk + i] = rkTmp[i];
    }
}

void RFE::execSortRankedListsAgg(std::vector<std::vector<int>>& rk) {
    std::cout << "\n\t Sort Ranked Lists (Fusion) ... \n";
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        kernelSortRankedListsAgg(i, rk);
    }
}

void RFE::kernelSortRankedListsAgg(int curRL, std::vector<std::vector<int>>& rk) {
    int LcurRL = rk[curRL].size();
    long int cNcurRL = ((long int) n)*curRL;
    float a[rk[curRL].size()];

    for (int j = 0; j < rk[curRL].size(); j++) {
        a[j] = matrix[cNcurRL + rk[curRL][j]];
    }

    //---------------------- INSERTION SORT --------------------------
    int i, j, keyR;
    float keyA;

    for (j = 1; j < rk[curRL].size(); j++) {
        keyA = a[j];
        keyR = rk[curRL][j];
        i = j - 1;
        while (i >= 0 && a[i] < keyA) {
            a[i + 1] = a[i];
            rk[curRL][i + 1] = rk[curRL][i];
            i--;
        }
        a[i + 1] = keyA;
        rk[curRL][i + 1] = keyR;
    }
    //----------------------------------------------------------------

    //Setting query image at first position
    i = 0;
    while ((rk[curRL][i] != curRL)&&(i < rk[curRL].size())) {
        i++;
    }
    if (i > 0) {
        int aux = rk[curRL][0];
        rk[curRL][0] = rk[curRL][i];
        rk[curRL][i] = aux;

        float auxA = a[0];
        a[0] = a[i];
        a[i] = auxA;
    }
}

bool RFE::search(int img, const std::vector<int>& rk) {
    for (int i = 0; i < rk.size(); i++) {
        if (img == rk[i]) {
            return true;
        }
    }
    return false;
}

float RFE::cosineSimCCs(const std::map<int, float>& di, const std::map<int, float>& dj) {
    float fsim = 0;
    float qi = 0;
    float qj = 0;
    for (auto it : di) {
        int k = it.first;
        float vi = it.second;
        float vj = dj.at(k);
        qi += vi*vi;
        qj += vj*vj;
        fsim += vi*vj;
    }
    if ((qi*qj) == 0) {
        return 0;
    }
    return (fsim / (sqrt(qi)*sqrt(qj)));
}
