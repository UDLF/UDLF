/* <HyperGraph.cpp>
 *
 * HyperGraph method implementation file
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

#include <iostream>
#include <omp.h>

#include "LHRR.hpp"

/* Constructor */
LHRR::LHRR() {

}

void LHRR::loadParameters() {
    Exec exec = Exec::getInstance();
    exec.getConfigVariable(l, "PARAM_LHRR_L");
    exec.getConfigVariable(k, "PARAM_LHRR_K");
    exec.getConfigVariable(t, "PARAM_LHRR_T");
}

void LHRR::checkParameters() {
    if (l > n) {
        std::cout << "L can't be greater than N" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
    if (k > l) {
        std::cout << "K can't be greater than L" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
}

void LHRR::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*l);
    initSparseMatrix(matrix);

    std::cout << "Initialized successfully!" << std::endl;
}

void LHRR::initDataStructuresFusion() {
    std::cout << "Initializing data structures..." << std::endl;

    initDataStructuresUdl(); //init the main structures

    //structures used for aggregation (fusion) tasks
    initSparseMatrix(matrixAgg);
    for (long int i = 0; i < n*n; i++) {
        matrixAgg[i] = 0;
    }

    //clear rks
    rkLists.clear();
    rkLists.resize(n*l);

    std::cout << "Initialized successfully!" << std::endl;
}

void LHRR::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
        initSparseMatrix(matrix);
    }
}

void LHRR::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void LHRR::runUdlMethod() {
    gettimeofday(&startTimeFillPos, NULL);
        execFillPosMatrix();
        execSortRankedLists();
    totalTimeFillPos = Time::addTime(startTimeFillPos, totalTimeFillPos);

    int iteration = 1;
    while (iteration <= t) {
        std::cout << std::endl << "  -> Iteration " << iteration << std::endl;
        hyperGraphIteration();
        iteration++;
    }
}

void LHRR::runFusionMethod() {
    for (std::string const& file : fusionFiles) {
        initSparseMatrix(matrix);
        for (long int i = 0; i < n*n; i++) {
            matrix[i] = 0;
        }

        //Read matrix from descriptor
        gettimeofday(&startTimeToDecrement, NULL);
            readInputFile(file);
            prepareInput();
        totalTimeToDecrement = Time::addTime(startTimeToDecrement, totalTimeToDecrement);

        //Run method for the current data
        runUdlMethod();

        //Store values in aggregation matrix
        for (long int i = 0; i < n; i++) {
            for (long int j = 0; j < l; j++) {
                int img = rkLists[l*i + j];
                matrixAgg[n*i + img] += (log(l+2)/log(j+2))*(confid[i]+1);
            }
        }
    }

    for (long int i = 0; i < n*n; i++) {
        matrix[i] = matrixAgg[i];
    }
    genRksFromSimMatrix();

    //Final execution
    t = 1;
    runUdlMethod();
}

void LHRR::initializeDataStructures() {
    hyperEdges.clear();
    revHyperEdges.clear();
    confid.clear();

    hyperEdges.resize(n);
    revHyperEdges.resize(n);
    confid.resize(n);

    tmpList.clear();
    tmpList.resize(n);

    imgList.clear();
    imgListRev.clear();
    valList.clear();
    valListRev.clear();

    imgList.resize(n);
    imgListRev.resize(n);
    valList.resize(n);
    valListRev.resize(n);
}

void LHRR::hyperGraphIteration() {
    gettimeofday(&startTimeInit, NULL);
        initializeDataStructures();
    totalTimeInit = Time::addTime(startTimeInit, totalTimeInit);

    gettimeofday(&startTimeLoadHE, NULL);
        loadHyperEdges();
        compressHE();
    totalTimeLoadHE = Time::addTime(startTimeLoadHE, totalTimeLoadHE);

    gettimeofday(&startTimeLoadRHE, NULL);
        loadRevHyperEdges();
    totalTimeLoadRHE = Time::addTime(startTimeLoadRHE, totalTimeLoadRHE);

    gettimeofday(&startTimeResetDB, NULL);
        initSparseMatrix(matrix);
        resetDB(1);
    totalTimeResetDB = Time::addTime(startTimeResetDB, totalTimeResetDB);

    gettimeofday(&startTimeProdCart, NULL);
        computeCartesianProductHyperEdges();
    totalTimeProdCart = Time::addTime(startTimeProdCart, totalTimeProdCart);

    gettimeofday(&startTimeSimHE, NULL);
        computeHyperEdgesSimilarities();
    totalTimeSimHE = Time::addTime(startTimeSimHE, totalTimeSimHE);

    gettimeofday(&startTimeSimRHE, NULL);
        computeReciprocalHyperEdgesSimilarities();
    totalTimeSimRHE = Time::addTime(startTimeSimRHE, totalTimeSimRHE);

    gettimeofday(&startTimeSort, NULL);
        computeDBBySimilarities();
    totalTimeSort = Time::addTime(startTimeSort, totalTimeSort);
}

void LHRR::compressHE() {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        //init structures
        std::vector<int> tmpImgList;
        std::vector<float_t> tmpValList;
        float tmpArr[n];
        for (int j = 0; j < n; j++) {
            tmpArr[j] = 0;
        }

        //fill array
        int j = 0;
        for (const int& img : imgList[i]) {
            if ((tmpArr[img] == 0) && (valList[i][j] != 0)){
                tmpImgList.push_back(img);
            }
            tmpArr[img] += valList[i][j];
            j++;
        }

        //fill values
        for (const int& img : tmpImgList) {
            tmpValList.push_back(tmpArr[img]);
        }

        //bind new lists
        imgList[i] = tmpImgList;
        valList[i] = tmpValList;
    }
}

void LHRR::computeCartesianProductHyperEdges() {
    std::cout << "\n\t + Cartesian Product Hyper Edges ...\n";
    for (int img = 0; img < n; img++) {
        double conf = confid[img];
        auto curHEdge = hyperEdges[img];
        for (auto entry1 : curHEdge) { //for the first k hyper edges or for all of them?
            int img1 = entry1.first;
            double value1 = entry1.second;
            for (auto entry2 : curHEdge) { //for the first k hyper edges or for all of them?
                int img2 = entry2.first;
                double value2 = entry2.second;
                long int index = ((long int) n)*img1 + img2;
                double inc = conf * value1 * value2;
                if (matrix[index] == 0) {
                    tmpList[img1].push_back(img2);
                    matrix[index] = 1;
                }
                matrix[index] += inc; //it is an increment, so the similarity increases
            }
        }
    }
}

void LHRR::sortTmpList(int qimg) {
    std::vector<std::pair<int, float>> rkTmp;

    for (int i = 0; i < tmpList[qimg].size(); i++) {
        int img = tmpList[qimg][i];
        long int index = ((long int) n)*qimg + img;
        rkTmp.push_back(std::make_pair(img, matrix[index]));
    }

    std::sort(rkTmp.begin(), rkTmp.end(), [](const std::pair<int, float> &x, const std::pair<int, float> &y) { //sorting
        return x.second > y.second;
    });

    for (int i = 0; i < tmpList[qimg].size(); i++) {
        tmpList[qimg][i] = rkTmp[i].first;
    }
}

void LHRR::computeHyperEdgesSimilarities() {
    std::cout << "\n\t + Similarity between Hyper Edges ...\n";
    #pragma omp parallel for
    for (int qimg = 0; qimg < n; qimg++) {

        int j = 0;
        float tmpArr[n];
        for (int i = 0; i < n; i++) {
            tmpArr[i] = 0;
        }
        for (const int& img : imgList[qimg]) {
            tmpArr[img] = valList[qimg][j];
            j++;
        }


        for (int j = 1; j < l; j++) {
            int img = rkLists[l*qimg + j];

            double simValue = 0;
            int o = 0;
            for (const int& e : imgList[img]) {
                simValue += tmpArr[e]*valList[img][o];
                o++;
            }

            long int index = ((long int) n)*qimg + img;
            double curValue = matrix[index];
            double newValue = curValue * simValue;

            matrix[index] = newValue;
        }
    }
}

void LHRR::computeReciprocalHyperEdgesSimilarities() {
    std::cout << "\n\t + Similarity between Reverse and Hyper Edges ...\n";
    #pragma omp parallel for
    for (int qimg = 0; qimg < n; qimg++) {
        int j = 0;
        float tmpArr[n];
        for (int i = 0; i < n; i++) {
            tmpArr[i] = 0;
        }
        for (const int& img : imgListRev[qimg]) {
            tmpArr[img] = valListRev[qimg][j];
            j++;
        }

        for (int j = 1; j < l; j++) {
            int img = rkLists[l*qimg + j];

            double simValue = 0;
            int o = 0;
            for (const int& e : imgList[img]) {
                simValue += tmpArr[e]*valList[img][o];
                o++;
            }

            long int index = ((long int) n)*qimg + img;
            double curValue = matrix[index];
            double newValue = curValue * simValue;
            matrix[index] = newValue;
        }
    }
}

double LHRR::multiplyMaps(std::vector<std::pair<int, double>>& he1, std::vector<std::pair<int, double>>& he2) {
    double acumValue = 0;
    for (auto entry1 : he1) {
        int img = entry1.first;
        int pos = searchPairByKey(img, he2);
        if (pos != -1) {
            double v1 = entry1.second;
            double v2 = he2[pos].second;
            acumValue += v1 * v2;
        }
    }
    return acumValue;
}

void LHRR::computeDBBySimilarities() {
    std::cout << "\n\t + Computing New DB by Similarities ...\n";

    for (int i = 0; i < n; i++) {
        sortTmpList(i);
    }

    execSortRankedLists();

    for (int i = 0; i < n; i++) {
        joinRks(i);
    }
}

void LHRR::resetDB(int value) {
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            long int index = ((long int) n)*i + rkLists[l*i + j];
            matrix[index] = value;
        }
    }
}

void LHRR::loadHyperEdges() {
    std::cout << "\n\t + Creating HyperEdges ...\n";
    for (int i = 0; i < n; i++) {
        createHyperEdge(i);
    }
}

void LHRR::createHyperEdge(int img) {
    includeHyperEdgeValue(img, img, 1);
    for (int o = 0; o < k; o++) {
        int imgo = rkLists[l*img + o];
        int poso = o + 1;
        for (int j = 0; j < k; j++) {
            int imgj = rkLists[l*imgo + j];
            int posj = j + 1;
            double w = weightPosition(poso) * weightPosition(posj);
            includeHyperEdgeValue(img, imgj, w);
        }
    }

    //sort
    std::sort(hyperEdges[img].begin(), hyperEdges[img].end(), [](const std::pair<int, double> &x, const std::pair<int, double> &y) { //sorting
        return x.second > y.second;
    });

    //compute confidence
    for (int i = 0; i < k; i++) {
        confid[img] += hyperEdges[img][i].second;
    }
}

void LHRR::loadRevHyperEdges() {
    std::cout << "\n\t + Creating Reverse HyperEdges ...\n";
    for (int img = 0; img < n; img++) {
        for (auto entry : hyperEdges[img]) {
            int curImg =  entry.first;
            double curValue = entry.second;
            //add elements
            imgListRev[curImg].push_back(img);
            valListRev[curImg].push_back(curValue);
        }
    }
}

double LHRR::weightPosition(int pos) {
    double logValue = log(pos) / log(k);
    return (1.0 - logValue);
}

int LHRR::searchPairByKey(int key, std::vector<std::pair<int, double>>& hyperEdge) {
    int i = 0;
    for (auto it : hyperEdge) {
        if (it.first == key) {
            return i;
        }
        i++;
    }
    return -1;
}

void LHRR::includeHyperEdgeValue(int i, int j, double value) {
    double curValue = 0;
    int pos = searchPairByKey(j,hyperEdges[i]);
    if (pos != -1) {
        curValue = hyperEdges[i][pos].second;
        hyperEdges[i][pos].second = curValue + value;
    } else {
        hyperEdges[i].push_back(std::make_pair(j,value));
    }

    //add elements
    imgList[i].push_back(j);
    valList[i].push_back(value);
}

void LHRR::execFillPosMatrix() {
    std::cout << "\n\t Fill Distance Matrix with Positions ... \n";
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        kernelFillPosMatrix(i);
    }
}

void LHRR::execSortRankedLists() {
    std::cout << "\n\t Sort Ranked Lists ... \n";
    #pragma omp parallel for
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists(i);
    }
}

void LHRR::kernelFillPosMatrix(int rk) {
    long int lrk = ((long int) l)*rk;
    long int nrk = ((long int) n)*rk;
    long int img;

    for (int pos = 0; pos < l; pos++) {
        img = rkLists[lrk + pos];
        matrix[nrk + img] += l - pos;
        matrix[n*img + rk] += l - pos;
    }
}

void LHRR::kernelSortRankedLists(int rk) {
    long int LcurRL = ((long int) l)*rk;
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

    for (j = 2; j < l; j++) {
        keyA = a[j];
        keyR = rkLists[LcurRL + j];
        i = j - 1;
        while (i > 0 && (a[i] < keyA)) {
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


void LHRR::joinRks(int qimg) {
    if (tmpList[qimg].size() == 0) {
        return;
    }


    float a1[l];
    long int LcurRL = ((long int) l)*qimg;
    long int cNcurRL = ((long int) n)*qimg;

    float a2[tmpList[qimg].size()];

    long int index;
    for (int j = 0; j < l; j++) {
        index = cNcurRL + rkLists[LcurRL + j];
        a1[j] = matrix[index];
    }

    std::vector<int> rkTmp;
    rkTmp.clear();

    for (int j = 0; j < tmpList[qimg].size(); j++) {
        index = cNcurRL + tmpList[qimg][j];
        a2[j] = matrix[index];
        rkTmp.push_back(tmpList[qimg][j]);
    }

    int j = 0;
    for (int i = 1; i < l; i++) {
        if (a2[j] > a1[i]) {
            for (int o = l-1; o > i; o--) {
                rkLists[LcurRL + o] = rkLists[LcurRL + o - 1];
            }
            rkLists[LcurRL + i] = tmpList[qimg][j];
            j++;
            if (j >= tmpList[qimg].size()) {
                break;
            }
        }
    }
}

void LHRR::sortAll(int qimg) {
    std::vector<std::pair<int, float>> rkTmp;

    for (int i = 0; i < l; i++) {
        int img = rkLists[l*qimg + i];
        long int index = ((long int) n)*qimg + img;
        rkTmp.push_back(std::make_pair(img, matrix[index]));
    }

    for (int i = 0; i < tmpList[qimg].size(); i++) {
        int img = tmpList[qimg][i];
        long int index = ((long int) n)*qimg + img;
        rkTmp.push_back(std::make_pair(img, matrix[index]));
    }

    std::sort(rkTmp.begin(), rkTmp.end(), [](const std::pair<int, float> &x, const std::pair<int, float> &y) { //sorting
        return x.second > y.second;
    });

    for (int i = 0; i < l; i++) {
        rkLists[l*qimg + i] = rkTmp[i].first;
    }
}
