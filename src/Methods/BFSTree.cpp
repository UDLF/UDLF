/* <BFSTree.cpp>
 *
 * BFSTree implementation file
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

#include <iostream>

#include "BFSTree.hpp"

/* Constructor */
BFSTree::BFSTree() {

}

void BFSTree::loadParameters() {
    Exec exec = Exec::getInstance();

    //BFSTree Parameters (read from the configuration file)
    exec.getConfigVariable(l, "PARAM_BFSTREE_L");
    exec.getConfigVariable(k, "PARAM_BFSTREE_K");
    topK = k;
    exec.getConfigVariable(correlation_metric, "PARAM_BFSTREE_CORRELATION_METRIC");

}

void BFSTree::checkParameters() {
    if (l > n) {
        std::cout << "L can't be greater than N" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
}

void BFSTree::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*l);
    initSparseMatrix(matrix);
    //initSparseMatrix(W);

    std::cout << "Initialized successfully!" << std::endl;
}

void BFSTree::initDataStructuresFusion() {
    std::cerr << "WARNING: The Fusion Mode is not implemented for this method! Aborting ...\n";
    exit(1);
}

void BFSTree::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
        initSparseMatrix(matrix);
    }
}

void BFSTree::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void BFSTree::runUdlMethod() {
    std::cout << "\n Executing BFSTree!\n\n";

    computeRankNormalization(1);
    computeRankNormalization(2);
    computeRankCorrelation();
    initSparseMatrix(matrix);
    acumBFSTree();
    normalizeWProb();
    computeWDiffusion();
    normalizeWProb();
    execSortRankedLists();
}

void BFSTree::runFusionMethod() {
    std::cerr << "WARNING: The Fusion Mode is not implemented for this method! Aborting ...\n";
    exit(1);
}

void BFSTree::normalizeWProb() {
    for (int i = 0; i < n; i++) {
        double acum = 1.0;
        for (int j = 0; j < n; j++) {
            acum += getMatrixElem(matrix, j, i);
        }
        for (int j = 0; j < n; j++) {
            setMatrixElem(matrix, j, i,
                          getMatrixElem(matrix, j, i) / acum);
        }
    }
}

void BFSTree::computeWDiffusion() {
    float* rc = NULL;
    initSparseMatrix(rc);

    std::cout << " - Computing W Diffusion L...\n";

    for (int imgI = 0; imgI < n; imgI++) {
        for (int posJ = 0; posJ < l; posJ++) {
            int imgJ = getRKElem(imgI, posJ);
            for (int posM = 0; posM < l; posM++) {
                int imgM = getRKElem(imgI, posM);
                float value = getMatrixElem(matrix, imgI, imgM)
                              * getMatrixElem(matrix, imgJ, imgM);
                setMatrixElem(rc, imgI, imgJ,
                              getMatrixElem(rc, imgI, imgJ) + value);
            }
        }
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            setMatrixElem(matrix, i, j,
                          getMatrixElem(rc, i, j));
        }
    }

    delete [] rc;
}

////////// Get and set elems //////////////

inline int BFSTree::getRKElem(long int query, long int pos) {
    long int index = ((long int) (l*query + pos));
    return rkLists[index];
}

inline float BFSTree::getMatrixElem(float* m, long int i, long int j) {
    long int index = ((long int) n*i + j);
    return m[index];
}

inline void BFSTree::setMatrixElem(float* m, long int i, long int j, float value) {
    long int index = ((long int) n*i + j);
    m[index] = value;
}

inline float BFSTree::minZero(float x, float y) {
    if (x == 0) {
        return y;
    }
    return std::min(x, y);
}

inline void BFSTree::incSimW(int img1, int img2, float newSim) {
    setMatrixElem(matrix, img1, img2,
                  getMatrixElem(matrix, img1, img2) + newSim);
}

///////////////////////////////////////

////////// RL TREE ////////////////////

void BFSTree::acumBFSTree() {
    std::cout << " - Computing RL Tree Increments ... \n";
    for (int qimg = 0; qimg < n; qimg++) {
        computeIncTree(qimg);
    }
}

void BFSTree::computeIncTree(int qImg) {
    const int incType = 2;
    std::vector<int> l1;
    std::vector<int> l2;
    std::vector<double> l2d;
    
    l1.push_back(qImg);

    for (int i = 1; i < k; i++) {
        int imgI = getRKElem(qImg, i);
        l1.push_back(imgI);

        float ircor = imgCor[qImg][imgI];
        if (incType != 2) {
            incSimW(qImg, imgI, ircor);
        }

        for (int j = 1; j < k; j++) {
            int imgJ = getRKElem(imgI, j);
            l2.push_back(imgJ);

            float jcor = ircor * imgCor[imgI][imgJ];
            l2d.push_back(jcor);
            if (incType != 2) {
                incSimW(qImg, imgJ, jcor);
            }
        }
    }

    if (incType != 1) {
        for (int img1a : l1) {
            for (int img1b : l1) {
                    incSimW(img1a, img1b,
                    imgCor[qImg][img1a]*imgCor[qImg][img1b]);
            }
        }

        for (int i = 0; i < l2.size(); i++) {
            for (int j = 0; j < l2.size(); j++) {
                    incSimW(l2[i], l2[j],
                    l2d[i]*l2d[j]);
            }
        }

        for (int i = 0; i < l1.size(); i++) {
            for (int j = 0; j < l2.size(); j++) {
                    incSimW(l1[i], l2[j],
                    imgCor[qImg][l1[i]]*l2d[j]);
            }
        }

    }
}

///////////////////////////////////////

////////// NORMALIZATION //////////////

void BFSTree::computeRankNormalization(int type)
{
    std::cout << " - Computing Rank Normalization " << type << " ...\n";

    switch (type) {
        case 1: 
            for (int qImg = 0; qImg < n; qImg++) {
                for (int pos = 0; pos < l; pos++) {
                    int img = getRKElem(qImg, pos);
                    setMatrixElem(matrix, qImg, img,
                                  getMatrixElem(matrix, qImg, img) + (l - pos));
                    setMatrixElem(matrix, img, qImg,
                                  getMatrixElem(matrix, img, qImg) + (l - pos));
                }
            }
            break;
        case 2:
            for (int qImg = 0; qImg < n; qImg++) {
                for (int pos = 0; pos < l; pos++) {
                    int img = getRKElem(qImg, pos);
                    setMatrixElem(matrix, qImg, img,
                                  minZero(getMatrixElem(matrix, qImg, img), (l - pos)));
                    setMatrixElem(matrix, img, qImg,
                                  minZero(getMatrixElem(matrix, img, qImg), (l - pos)));
                }
            } 
            break;
    }

    execSortRankedLists();
}

void BFSTree::computeRankNormalization2(int type)
{
    std::cout << " - Computing Rank Normalization " << type << " ...\n";

    switch (type) {
        case 1: 
            for (int qImg = 0; qImg < n; qImg++) {
                for (int pos = 0; pos < l; pos++) {
                    int img = getRKElem(qImg, pos);
                    setMatrixElem(matrix, qImg, img,
                                  getMatrixElem(matrix, qImg, img) + (pos));
                    setMatrixElem(matrix, img, qImg,
                                  getMatrixElem(matrix, img, qImg) + (pos));
                }
            }
            break;
        case 2:
            for (int qImg = 0; qImg < n; qImg++) {
                for (int pos = 0; pos < l; pos++) {
                    int img = getRKElem(qImg, pos);
                    setMatrixElem(matrix, qImg, img,
                                  std::max(getMatrixElem(matrix, qImg, img), ((float)pos)));
                    setMatrixElem(matrix, img, qImg,
                                  std::max(getMatrixElem(matrix, img, qImg), ((float)pos)));
                }
            } 
            break;
    }

    execSortRankedLists2();
}

////////////////////////////

////////// SORTING /////////

void BFSTree::execSortRankedLists()
{
    std::cout << "\n\t Sort Ranked Lists ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists(i);
    }
}

void BFSTree::kernelSortRankedLists(int rk)
{
    int l = rkLists.size()/n;
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

void BFSTree::execSortRankedLists2()
{
    std::cout << "\n\t Sort Ranked Lists ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists2(i);
    }
}

void BFSTree::kernelSortRankedLists2(int rk)
{
    int l = rkLists.size()/n;
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

////////////////////////////


// ---------------- METRICS ------------------------- //

void BFSTree::computeRankCorrelation()
{
    std::cout << "\t- Computing Rank Correlation ...\n";

    imgCor.clear();
    imgCor.resize(n);

    for (int qImg = 0; qImg < n; qImg++) {
        imgCor[qImg][qImg] = 1.0;
        for (int pos = 0; pos < k; pos++) {
            int img = getRKElem(qImg, pos);
            float correlation = rankCorrelationMeasure(qImg, img);
            imgCor[qImg][img] = correlation;
        }
    }
}

float BFSTree::rankCorrelationMeasure(int i1, int i2) {
    //Compute similarity
    float simRL = 0;
    if (correlation_metric == "INTERSECTION") {
        simRL = intersection(i1, i2);
    } else if (correlation_metric == "RBO") {
        simRL = rbo(i1, i2);
    } else if (correlation_metric == "KENDALL_TAU") {
        simRL = kendallTau(i1, i2);
    } else if (correlation_metric == "KENDALL_TAU_W") {
        simRL = kendallTauW(i1, i2);
    } else if (correlation_metric == "JACCARD") {
        simRL = jaccard(i1, i2);
    } else if (correlation_metric == "JACCARD_K") {
        simRL = jaccardK(i1, i2);
    } else if (correlation_metric == "SPEARMAN") {
        simRL = spearman(i1, i2);
    } else if (correlation_metric == "GOODMAN") {
        simRL = goodman(i1, i2);
    } else {
        std::cout << " Metric "
                  << correlation_metric
                  << " not implemented! Aborting ...\n";
        exit(1);
    }

    return simRL;
}

float BFSTree::intersection(int i1, int i2) {
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

float BFSTree::rbo(int i1, int i2) {
    //Initialize RBO Variables
    float p = 0.7;
    float inter = 0;
    float result = 0;
    int k = 1;

    //Compute RBO
    while (k <= topK) {
        inter = 0;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                if (getRKElem(i1, i) == getRKElem(i2, j)) {
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

float BFSTree::kendallTau(int i1, int i2) {
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

float BFSTree::kendallTauW(int i1, int i2) {
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

float BFSTree::jaccard(int i1, int i2) {
    //Initialize Jaccard Variables
    float jaccard = 0;
    float interSize = 0;
    float unionSize = 0;

    //Compute Intersection
    for (int i = 0; i < topK; i++) {
        for (int j = 0; j < topK; j++) {
            if (getRKElem(i1, i) == getRKElem(i2, j)) {
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

float BFSTree::jaccardK(int i1, int i2) {
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

float BFSTree::spearman(int i1, int i2) {
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

float BFSTree::goodman(int i1, int i2) {
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
