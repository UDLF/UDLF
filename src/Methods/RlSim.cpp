/* <RlSim.cpp>
 *
 * RlSim method implementation file
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

#include <iostream>

#include "RlSim.hpp"

/* Constructor */
RlSim::RlSim() {

}

void RlSim::loadParameters() {
    Exec exec = Exec::getInstance();
    exec.getConfigVariable(topK,   "PARAM_RLSIM_TOPK");
    exec.getConfigVariable(cK,     "PARAM_RLSIM_CK");
    exec.getConfigVariable(t,      "PARAM_RLSIM_T");
    exec.getConfigVariable(metric, "PARAM_RLSIM_METRIC");
}

void RlSim::checkParameters() {
    if (cK > n) {
        std::cout << "cK can't be greater than N" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
}

void RlSim::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*n);
    initMatrix(matrix);
    initMatrix(tmpMatrix);

    std::cout << "Initialized successfully!" << std::endl;
}

void RlSim::initDataStructuresFusion() {
    initDataStructuresUdl();
}

void RlSim::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
    } else { //RK
        //not really necessary since the sorting is stable
        genDistMatrixFromRks();
    }
}

void RlSim::prepareInputFusion() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "SIM") {
            convertSimToDistMatrix();
        }
    } else { //RK
        genDistMatrixFromRks();
    }
}

void RlSim::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void RlSim::runUdlMethod() {
    for (int curIteration = 0; curIteration < t; curIteration++) {
        std::cout << "\n + RL-Sim:\tRunning iteration " << (curIteration + 1) << "  \n";

        resetTmpMatrix();
        execUpdateDistances();
        normalizeMinDistances();
        execSortRankedLists();

        topK++;
    }

    releaseTmpMatrix();
}

void RlSim::runFusionMethod() {
    prepareRankAggregationMatrix();
    genRksFromDistMatrix();

    runUdlMethod();
}

void RlSim::execUpdateDistances() {
    std::cout << "\n\t Update Distances ... \n";
    for (int i = 0; i < n; i++) {
        kernelUpdateDistances(i);
    }
}

void RlSim::execSortRankedLists() {
    std::cout << "\n\t Sort Ranked Lists ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists(i);
    }
}

void RlSim::kernelUpdateDistances(int i1) {
    int Ni1 = n*i1;

    for (int j = 1; j < cK + 1; j++) {
        int i2 = rkLists[Ni1 + j];

        //Compute similarity
        float simRL = 0;
        if (metric == "INTERSECTION") {
            simRL = intersection(i1, i2);
        } else if (metric == "RBO") {
            simRL = rbo(i1, i2);
        } else if (metric == "KENDALL_TAU") {
            simRL = kendallTau(i1, i2);
        } else if (metric == "KENDALL_TAU_W") {
            simRL = kendallTauW(i1, i2);
        } else if (metric == "JACCARD") {
            simRL = jaccard(i1, i2);
        } else if (metric == "JACCARD_K") {
            simRL = jaccardK(i1, i2);
        } else if (metric == "SPEARMAN") {
            simRL = spearman(i1, i2);
        } else if (metric == "GOODMAN") {
            simRL = goodman(i1, i2);
        } else {
            std::cout << " Metric " << metric << " not implemented! Aborting ...\n";
            exit(1);
        }

        //Attribute the similarity result as distance
        if (simRL > 0) {
            tmpMatrix[Ni1 + i2] = 1/(1 + simRL);          //first segment
        } else {
            tmpMatrix[Ni1 + i2] = matrix[Ni1 + i2] + 1;   //second segment
        }
    }

    //Third segment
    for (int j = cK + 1; j < n; j++) {
        int i2 = rkLists[Ni1 + j];
        tmpMatrix[Ni1 + i2] = matrix[Ni1 + i2] + 2; //sum the matrix is not really necessary since the sorting is stable
    }

    tmpMatrix[Ni1 + i1] = 0; //the distance of i1 from i1 is always 0
}

void RlSim::kernelSortRankedLists(int rk) {
    int l = rkLists.size()/n;
    int LcurRL = l*rk;
    int cNcurRL = n*rk;
    float a[l];

    for (int j = 0; j < l; j++) {
        a[j] = matrix[cNcurRL + rkLists[LcurRL + j]];
    }

    //---------------------- INSERTION SORT --------------------------
    int i, j, keyR;
    float keyA;

    for (j = 1; j < l; j++) {
        keyA = a[j];
        keyR = rkLists[LcurRL + j];
        i = j - 1;
        while (i >= 0 && a[i] > keyA) {
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

void RlSim::resetTmpMatrix() {
    initMatrix(tmpMatrix);
}

void RlSim::releaseTmpMatrix() {
    delete [] tmpMatrix;
}

void RlSim::normalizeMinDistances() {
    std::cout << "\n\t Normalize Distances ... \n";
    float normDist;
    unsigned int Ni, indexIJ, indexJI;
    for (int i = 0; i < n; i++) {
        Ni = n * i;
        for (int j = i; j < n; j++) {
            indexIJ = Ni + j;
            indexJI = n * j + i;
            normDist = std::min(tmpMatrix[indexIJ], tmpMatrix[indexJI]);
            matrix[indexIJ] = normDist;
            matrix[indexJI] = normDist;
        }
    }
}

void RlSim::prepareRankAggregationMatrix() {
    std::cout << "\n\t Preparing Aggregation Matrix ... \n";

    //Init the temporary matrix
    for (long int i = 0; i < n*n; i++) {
        tmpMatrix[i] = 1;
    }

    for (std::string const& file : fusionFiles) {
        //Read matrix from descriptor
        gettimeofday(&startTimeToDecrement, NULL);
            readInputFile(file);
            prepareInputFusion();
        totalTimeToDecrement = Time::addTime(startTimeToDecrement, totalTimeToDecrement);

        //Multiplying
        for (long int i = 0; i < n*n; i++) {
            tmpMatrix[i] = tmpMatrix[i] * (1 + matrix[i]);
        }
    }

    //Root
    float raExp = 1 / ((float) fusionFiles.size());
    for (long int i = 0; i < n*n; i++) {
        tmpMatrix[i] = pow(tmpMatrix[i], raExp);
    }

    //Redirect Matrix
    delete [] matrix;
    matrix = tmpMatrix;
    tmpMatrix = NULL;
}

// ---------------- METRICS ------------------------- //

float RlSim::intersection(int i1, int i2) {
    //If it is the same image, the rankings are identical
    if (i1 == i2) {
        return 1;
    }

    //Initialize sets
    std::vector<int> set1(topK);
    std::vector<int> set2(topK);
    set1[0] = i1;
    set2[0] = i2;
    int setSize = 1; //starting with the first image in it

    //Initialize Intersection Variables
    int acumInter = 0;
    int curInterSize = 0;
    int Ni1 = n*i1;
    int Ni2 = n*i2;

    //Initialize the kNN sets (rkk1 and rkk2)
    std::vector<int> rkk1(topK);
    std::vector<int> rkk2(topK);

    //Build tradicional kNNs
    for (int i = 1; i < topK; i++) {
        rkk1[i] = rkLists[Ni1 + i];
    }
    for (int i = 1; i < topK; i++) {
        rkk2[i] = rkLists[Ni2 + i];
    }

    //Compute Intersection
    for (int i = 1; i < topK; i++) {
        int iRL1 = rkk1[i];
        int iRL2 = rkk2[i];

        if (iRL1 == iRL2) {
            curInterSize++;
        } else {
            set1[setSize] = iRL1;
            set2[setSize] = iRL2;
            setSize++;
            for (int j = 0; j < setSize; j++) {
                if (set1[j] == iRL2) {
                    curInterSize++;
                    break;
                }
            }
            for (int j = 0; j < setSize; j++) {
                if (set2[j] == iRL1) {
                    curInterSize++;
                    break;
                }
            }
        }

        acumInter += curInterSize;
    }

    //The intersection value is computed in the interval [0,1]
    float interSim = ((float) acumInter) / ((float) (topK*(topK+1))/2);

    return interSim;
}

float RlSim::rbo(int i1, int i2) {
    //Initialize RBO Variables
    float p = 0.9;
    float inter = 0;
    float result = 0;
    int k = 1;

    //Compute RBO
    while (k <= topK) {
        inter = 0;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                if (rkLists[n*i1 + i] == rkLists[n*i2 + j]) {
                    inter += 1;
                }
            }
        }
        result += pow(p, k-1)*(inter/k);
        k++;
    }

    result = (1-p)*result;

    return result;
}

float RlSim::kendallTau(int i1, int i2) {
    //Initialize KendallTau Variables
    float ktau = 0;
    int maxSize = topK*2;
    int vectorIdxRkl1[maxSize];
    int vectorIdxRkl2[maxSize];
    int unionSet[maxSize];
    int unionSize = 0;

    //Compute the union of the kNNs
    for (int i = 0; i < topK; i++) {
        unionSet[i] = rkLists[n*i1 + i];
        unionSize++;
    }
    for (int i = 0; i < topK; i++) {
        int j = 0;

        while (rkLists[n*i2 + i] != unionSet[j] && j < unionSize) {
            j++;
        }

        if (j == unionSize) {
            unionSet[unionSize] = rkLists[n*i2 + i];
            unionSize++;
        }
    }

    //If the kNNs are totally different (0 elements in common), return 0 (no similarity)
    if (unionSize == maxSize) {
        return 0;
    }

    //Compute the indexes
    for (int i = 0; i < unionSize; i++) {
        int j = 0;
        while (rkLists[n*i1 + j] != unionSet[i]) {
            j++;
        }
        vectorIdxRkl1[i] = j;

        j = 0;
        while (rkLists[n*i2 + j] != unionSet[i]) {
            j++;
        }
        vectorIdxRkl2[i] = j;
    }

    //Count the number of discordant pairs
    for (int i = 0; i < unionSize; i++) {
        for (int j = i + 1; j < unionSize; j++) {
            int comp1 = vectorIdxRkl1[i] >= vectorIdxRkl1[j];
            int comp2 = vectorIdxRkl2[i] >= vectorIdxRkl2[j];

            if (comp1 != comp2) {
                ktau += 1;
            }
        }
    }

    //Normalize in the [0,1] interval
    ktau = ktau/(topK*(topK-1));

    //Convert distance to similarity
    ktau = 1/(1 + ktau);

    return ktau;
}

float RlSim::kendallTauW(int i1, int i2) {
    //Initialize KendallTauW Variables
    float p = 0.8;
    float ktau = 0;
    int maxSize = topK*2;
    int vectorIdxRkl1[maxSize];
    int vectorIdxRkl2[maxSize];
    int unionSet[maxSize];
    int unionSize = 0;

    //Compute the union of the kNNs
    for (int i = 0; i < topK; i++) {
        unionSet[i] = rkLists[n*i1 + i];
        unionSize++;
    }
    for (int i = 0; i < topK; i++) {
        int j = 0;

        while (rkLists[n*i2 + i] != unionSet[j] && j < unionSize) {
            j++;
        }

        if (j == unionSize) {
            unionSet[unionSize] = rkLists[n*i2 + i];
            unionSize++;
        }
    }

    //If the kNNs are totally different (0 elements in common), return 0 (no similarity)
    if (unionSize == maxSize) {
        return 0;
    }

    //Compute the indexes
    for (int i = 0; i < unionSize; i++) {
        int j = 0;
        while (rkLists[n*i1 + j] != unionSet[i]) {
            j++;
        }
        vectorIdxRkl1[i] = j;

        j = 0;
        while (rkLists[n*i2 + j] != unionSet[i]) {
            j++;
        }
        vectorIdxRkl2[i] = j;
    }

    //Count the number of discordant pairs applying weights
    for (int i = 0; i < unionSize; i++) {
        for (int j = i + 1; j < unionSize; j++) {
            bool comp1 = (vectorIdxRkl1[i] >= vectorIdxRkl1[j]);
            bool comp2 = (vectorIdxRkl2[i] >= vectorIdxRkl2[j]);
            if (comp1 != comp2) {
                int minPos1 = std::min(vectorIdxRkl1[i], vectorIdxRkl1[j]);
                int minPos2 = std::min(vectorIdxRkl2[i], vectorIdxRkl2[j]);
                float minPos = std::min(minPos1, minPos2);
                float diff1 = abs(vectorIdxRkl1[i] - vectorIdxRkl1[j]);
                float diff2 = abs(vectorIdxRkl2[i] - vectorIdxRkl2[j]);

                float sizeMax = topK;
                float minValue = 2;
                float multFactor = std::min(minValue, ((diff1 + diff2) / sizeMax));

                ktau += pow(p, minPos) * multFactor;
            }
        }
    }

    //Normalize in the [0,1] interval
    ktau = ktau/(2*topK*(topK-1));

    //Convert distance to similarity
    ktau = 1/(1 + ktau);

    return ktau;
}

float RlSim::jaccard(int i1, int i2) {
    //Initialize Jaccard Variables
    float jaccard = 0;
    float interSize = 0;
    float unionSize = 0;

    //Compute Intersection
    for (int i = 0; i < topK; i++) {
        for (int j = 0; j < topK; j++) {
            if (rkLists[n*i1 + i] == rkLists[n*i2 + j]) {
                interSize += 1;
                break;
            }
        }
    }

    //Compute Union
    unionSize = (2*topK) - interSize;

    //If the kNNs are totally different (0 elements in common), return 0 (no similarity)
    if (unionSize == 2*topK) {
        return 0;
    }

    //Calculate Jaccard Index
    jaccard = interSize/unionSize;

    return jaccard;
}

float RlSim::jaccardK(int i1, int i2) {
    //Initialize JaccardK Variables
    float jaccard = 0;
    float interSize = 0;
    float unionSize = 0;
    int k = 1;

    //Calculate JaccardK Index
    while (k <= topK) {
        interSize = 0;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                if (rkLists[n*i1 + i] == rkLists[n*i2 + j]) {
                    interSize += 1;
                    break;
                }
            }
        }
        unionSize = (2*k) - interSize;
        jaccard += (interSize/unionSize);
        k++;
    }

    jaccard = jaccard/topK;

    return jaccard;
}

float RlSim::spearman(int i1, int i2) {
    //Initialize Spearman Variables
    float spearman = 0;
    int maxSize = topK*2;
    int unionSet[maxSize];
    int vectorIdxRkl1[maxSize];
    int vectorIdxRkl2[maxSize];
    int unionSize = 0;

    //Compute the union of the kNNs
    for (int i = 0; i < topK; i++) {
        unionSet[i] = rkLists[n*i1 + i];
        unionSize++;
    }
    for (int i = 0; i < topK; i++) {
        int j = 0;

        while (rkLists[n*i2 + i] != unionSet[j] && j < unionSize) {
            j++;
        }

        if (j == unionSize) {
            unionSet[unionSize] = rkLists[n*i2 + i];
            unionSize++;
        }
    }

    //If the kNNs are totally different (0 elements in common), return 0 (no similarity)
    if (unionSize == maxSize) {
        return 0;
    }

    //Compute the indexes
    for (int i = 0; i < unionSize; i++) {
        int j = 0;
        while (rkLists[n*i1 + j] != unionSet[i]) {
            j++;
        }
        vectorIdxRkl1[i] = j;

        j = 0;
        while (rkLists[n*i2 + j] != unionSet[i]) {
            j++;
        }
        vectorIdxRkl2[i] = j;
    }

    //Calculate the distance between ranks
    for (int i = 0; i < unionSize; i++) {
        spearman += abs(vectorIdxRkl1[i] - vectorIdxRkl2[i]);
    }

    //Normalize in the [0,1] interval
    spearman = spearman/(2*topK*n);

    //Convert distance to similarity
    spearman = 1/(1 + spearman);

    return spearman;
}

float RlSim::goodman(int i1, int i2) {
    //Initialize Goodman Variables
    float goodman;
    int ns = 0;
    int nd = 0;
    int maxSize = topK*2;
    int unionSet[maxSize];
    int vectorIdxRkl1[maxSize];
    int vectorIdxRkl2[maxSize];
    int unionSize = 0;

    //Compute the union of the kNNs
    for (int i = 0; i < topK; i++) {
        unionSet[i] = rkLists[n*i1 + i];
        unionSize++;
    }
    for (int i = 0; i < topK; i++) {
        int j = 0;

        while (rkLists[n*i2 + i] != unionSet[j] && j < unionSize) {
            j++;
        }

        if (j == unionSize) {
            unionSet[unionSize] = rkLists[n*i2 + i];
            unionSize++;
        }
    }

    //If the kNNs are totally different (0 elements in common), return 0 (no similarity)
    if (unionSize == maxSize) {
        return 0;
    }

    //Compute the indexes
    for (int i = 0; i < unionSize; i++) {
        int j = 0;
        while (rkLists[n*i1 + j] != unionSet[i]) {
            j++;
        }
        vectorIdxRkl1[i] = j;

        j = 0;
        while (rkLists[n*i2 + j] != unionSet[i]) {
            j++;
        }
        vectorIdxRkl2[i] = j;
    }

    //Calculate Goodman
    for (int i = 0; i < unionSize; i++) {
        for (int j = i + 1; j < unionSize; j++) {
            bool comp1 = (vectorIdxRkl1[i] >= vectorIdxRkl1[j]);
            bool comp2 = (vectorIdxRkl2[i] >= vectorIdxRkl2[j]);

            if (comp1 != comp2) {
                nd += 1;
            } else {
                ns += 1;
            }
        }
    }

    //Normalize in the [0,1] interval and return similarity
    goodman = ((float) (ns - nd)) / ((float) (ns + nd));
    goodman = (goodman + 1) / 2;

    return goodman;
}

// -------------------------------------------------- //
