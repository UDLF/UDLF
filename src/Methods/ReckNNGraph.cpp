/* <ReckNNGraph.cpp>
 *
 * ReckNNGraph method implementation file
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

#include <iostream>

#include "ReckNNGraph.hpp"

/* Constructor */
ReckNNGraph::ReckNNGraph() {

}

void ReckNNGraph::loadParameters() {
    Exec exec = Exec::getInstance();
    exec.getConfigVariable(l,       "PARAM_RECKNNGRAPH_L");
    exec.getConfigVariable(k,       "PARAM_RECKNNGRAPH_K");
    exec.getConfigVariable(epsilon, "PARAM_RECKNNGRAPH_EPSILON");
}

void ReckNNGraph::checkParameters() {
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

void ReckNNGraph::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*l);
    initSparseMatrix(matrix);

    std::cout << "Initialized successfully!" << std::endl;
}

void ReckNNGraph::initDataStructuresFusion() {
    initDataStructuresUdl();
}

void ReckNNGraph::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
    }
}

void ReckNNGraph::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void ReckNNGraph::runUdlMethod() {
    std::cout << "\n * Reciprocal kNN Graph: Starting (...) \n";

    initializeVariables();
    do {
        std::cout << "\n\n************************************************************\n";
        std::cout << "             >> ITERATION [" << curIterationReckNN << "]: ReckNN Graph!";
        std::cout << "\n************************************************************\n";
        initializeDB();
        iterationReciprocalkNNGraph();
        computeNewDB();
        k++;
        diffConvScore = convScore - prevConvScore;
        prevConvScore = convScore;
        std::cout << "\n\n   <o> Convergence Score Difference: " << diffConvScore;
        curIterationReckNN++;
    } while (diffConvScore > epsilon);
}

void ReckNNGraph::runFusionMethod() {
    std::cout << "\n * Reciprocal kNN Graph: Starting (...) \n";

    //read the first file
    gettimeofday(&startTimeToDecrement, NULL);
        readInputFile(fusionFiles[0]);
        prepareInput();
    totalTimeToDecrement = Time::addTime(startTimeToDecrement, totalTimeToDecrement);

    initializeVariables();
    do {
        std::cout << "\n\n************************************************************\n";
        std::cout << "             >> ITERATION [" << curIterationReckNN << "]: ReckNN Graph!";
        std::cout << "\n************************************************************\n";
        initializeDB();
        iterationReciprocalkNNGraph();
        computeNewDB();
        if (curIterationReckNN == 0) {
            runRankAggregation();
        }
        k++;
        diffConvScore = convScore - prevConvScore;
        prevConvScore = convScore;
        std::cout << "\n\n   <o> Convergence Score Difference: " << diffConvScore;
        curIterationReckNN++;
    } while (diffConvScore > epsilon);
}

void ReckNNGraph::runRankAggregation() {
    std::cout << "\n\n * Starting rank aggregation! \n";
    int cont = 0;
    for (std::string const& file : fusionFiles) {
        cont++;

        //skip the first file
        if (cont == 1) {
            continue;
        }

        //init tmpMatrix
        delete [] tmpMatrix;
        tmpMatrix = new float[n*n];

        //copy the content of the main matrix to the tmpMatrix
        copyDB();

        //read the current file
        gettimeofday(&startTimeToDecrement, NULL);
            readInputFile(file);
            prepareInput();
        totalTimeToDecrement = Time::addTime(startTimeToDecrement, totalTimeToDecrement);

        //Re-Ranking
        initializeDB();
        iterationReciprocalkNNGraph();
        computeNewDB();

        //Aggregate
        aggregate();
    }
}

void ReckNNGraph::copyDB() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float curDist = getDistance(i, j);
            if (curDist != 0) {
                setTempNewDistance(i, j, curDist);
            }
        }
    }
}

void ReckNNGraph::aggregate() {
    std::cout << "\n\n<AGGREGATE>\n\n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            float dist = getDistance(i, j);
            float tempDist = getTempDistance(i, j);
            if ((dist != 0) || (tempDist != 0)) {
                dist = checkAggregateDistance(dist);
                tempDist = checkAggregateDistance(tempDist);
                float newDist = dist * tempDist;
                setNewDistance(i, j, newDist);
            }
        }
    }
    updateModels();
}

float ReckNNGraph::checkAggregateDistance(float curDist) {
    if (curDist == 0) {
        return l;
    } else {
        return curDist;
    }
}

void ReckNNGraph::initializeVariables() {
    std::cout << "\t - Intializing general structures ...\n";
    convScore = 0;
    prevConvScore = 0;
    diffConvScore = 0;
    curIterationReckNN = 0;
    std::cout << "\t - Intializing general structures ... Done!\n";
}

void ReckNNGraph::initializeMatrix() {
    std::cout << "\t - Intializing Matrix A ...\n";
    unsigned int n2 = n*n;
    delete [] matrix;
    matrix = new float[n2];
    for (int l = 0; l < n2; l++) {
        if (matrix[l] != 0) {
            matrix[l] = 0;
        }
    }
    std::cout << "\t - Intializing Matrix A ... Done!\n";
}

void ReckNNGraph::initializePositions() {
    std::cout << "\t - Intializing Positions... \n";
    unsigned int n2 = n*n;
    delete [] posMatrix;
    posMatrix = new unsigned int[n2];
    for (int l = 0; l < n2; l++) {
        if (posMatrix[l] != 0) {
            posMatrix[l] = 0;
        }
    }
    for (int i = 0; i < n; i++) {
        unsigned int Ni = n*i;
        unsigned int Li = l*i;
        for (int j = 0; j < l; j++) {
            posMatrix[Ni + rkLists[Li + j]] = j + 1;
        }
    }
    std::cout << "\t - Intializing Positions ... Done!\n";
}

void ReckNNGraph::initializeDB() {
    std::cout << "\t - Intializing DB ... \n";
    initializeMatrix();
    initializePositions();
}

void ReckNNGraph::iterationReciprocalkNNGraph() {
    std::cout << "\n\n - Iteration Reciprocal kNN Graph ... \n\n";
    double acumScore = 0;
    for (int imgId = 0; imgId < n; imgId++) {
        for (knn = 1; knn < k; knn++) {
            double scoreRL = evaluateRankedList(imgId, 0);
            acumScore += scoreRL;
        }
    }
    acumScore = acumScore / ((double) (k * n));
    convScore = acumScore;
    std::cout << "   + Convergence Value: " << acumScore;
}

double ReckNNGraph::evaluateRankedList(int imgId, float scoreThreshold) {
    std::vector<unsigned int> set(knn + 1);
    for (int i = 0; i <= knn; i++) {
        set[i] = getImgAtPos(imgId, i);
    }
    int size = knn + 1;
    double score = evaluateSetByCliqueNonWeithted(set, size);
    computeIncrements(score, set, scoreThreshold, size);
    return score;
}

double ReckNNGraph::evaluateSetByCliqueNonWeithted(const std::vector<unsigned int>& set, int size) {
    double score = 0;
    double maxScore = pow(size, 2);
    unsigned int imgId1;
    for (int i = 0; i < size; i++) {
        imgId1 = set[i];
        for (int j = 0; j < size; j++) {
            int imgId2 = getImgAtPos(imgId1, j);
            for (int l = 0; l < size; l++) {
                if (set[l] == imgId2) {
                    score++;
                }
            }
        }
    }
    score = score / maxScore;
    return score;
}

void ReckNNGraph::computeIncrements(double score, const std::vector<unsigned int>& set, double scoreThreshold, int size) {
    if (score < scoreThreshold) {
        return;
    }
    float squaredScore = pow(score, 2);
    unsigned int imgId1, imgId2;
    for (int i = 0; i < size; i++) {
        imgId1 = set[i];
        for (int j = 0; j < size; j++) {
            imgId2 = set[j];
            float curDist = getDistance(imgId1, imgId2);
            float increment = squaredScore;
            curDist += increment;
            setNewDistance(imgId1, imgId2, curDist);
        }
    }
}

void ReckNNGraph::computeNewDB() {
    std::cout << "\n\n + Computing new distances ... ";
    unsigned int imgId1, imgId2, pos;
    for (int i = 0; i < n; i++) {
        imgId1 = i;
        for (int j = 0; j < n; j++) {
            imgId2 = j;
            float currentDist = getDistance(imgId1, imgId2);
            pos = getPosition(imgId1,imgId2);
            if (currentDist == 0) {
                if (pos < l) {
                    setNewDistance(imgId1, imgId2, pos);
                }  else {
                    setNewDistance(imgId1, imgId2, l);
                }
            } else {
                float recKNNPos = std::max(getPosition(imgId1, imgId2), getPosition(imgId2, imgId1));
                float posDistance = recKNNPos / ((float) (l));
                float newDist = posDistance / (1 + currentDist);
                setNewDistance(imgId1, imgId2, newDist);
            }
        }
        setNewDistance(imgId1, imgId1, 0);
    }
    updateModels();
}

void ReckNNGraph::updateModels() {
    std::cout << "\n\n [*] Updating Models ...";
    std::cout << "\n\t + ReckNN:    Sorting ranked lists ... \n";
    for (int i = 0; i < n; i++) {
        sortRankedList(i);
    }
}

void ReckNNGraph::sortRankedList(int curRL) {
    int NRL = n*curRL;
    int LRL = l*curRL;

    float a[n];
    unsigned int r[n];

    std::vector<unsigned int> processed(n);
    int sRK = l;

    r[0] = curRL;
    a[0] = 0;
    processed[curRL] = 1;
    for (int j = 1; j < l; j++) {
        r[j] = rkLists[LRL + j];
        float curDist = matrix[NRL + r[j]];
        if (curDist == 0) {
            curDist = j;
        }
        a[j] = curDist;
        processed[r[j]] = 1;
    }

    for (int j = 0; j < n; j++) {
        float curDist = matrix[NRL + j];
        if ((processed[j] != 1)&&(curDist != 0)) {
            r[sRK] = j;
            a[sRK] = curDist;
            sRK++;
        }
    }

    //---------------------- INSERTION SORT --------------------------
    int i, j, keyR;
    float keyA;

    for (j = 1; j < sRK; j++) {
        keyA = a[j];
        keyR = r[j];
        i = j - 1;
        while (i >= 0 && a[i] > keyA) {
            a[i + 1] = a[i];
            r[i + 1] = r[i];
            i--;
        }
        a[i + 1] = keyA;
        r[i + 1] = keyR;
    }
    //----------------------------------------------------------------

    //Setting query image at first position
    i = 0;
    while ((r[i] != curRL)&&(i < sRK)) {
        i++;
    }
    if ((i > 0)&&(i != sRK)) {
        int aux = r[0];
        r[0] = r[i];
        r[i] = aux;
    }

    //Redefining rkLists
    for (int j = 0; j < l; j++) {
        rkLists[LRL + j] = r[j];
    }
}

unsigned int ReckNNGraph::getImgAtPos(unsigned int qImg, unsigned int iPos) {
    return rkLists[l * qImg + iPos];
}

float ReckNNGraph::getDistance(unsigned int imgId1, unsigned int imgId2) {
    return matrix[n * imgId1 + imgId2];
}

void ReckNNGraph::setNewDistance(unsigned int imgId1, unsigned int imgId2, float newDist) {
    matrix[n * imgId1 + imgId2] = newDist;
}

float ReckNNGraph::getTempDistance(unsigned int imgId1, unsigned int imgId2) {
    return tmpMatrix[n * imgId1 + imgId2];
}

void ReckNNGraph::setTempNewDistance(unsigned int imgId1, unsigned int imgId2, float newDist) {
    tmpMatrix[n * imgId1 + imgId2] = newDist;
}

unsigned int ReckNNGraph::getPosition(unsigned int imgId1, unsigned int imgId2) {
    unsigned int pos = posMatrix[n*imgId1 + imgId2];
    if (pos == 0) {
        return l;
    } else {
        return pos;
    }
}
