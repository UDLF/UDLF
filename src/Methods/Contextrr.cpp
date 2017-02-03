/* <Contextrr.cpp>
 *
 * Contextrr method implementation file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 ******************************************************************************************
 *
 * ContextRR is presented in the paper:
 *   D. C. G. Pedronette and R. da S. Torres.
 *   "Exploiting contextual information for image re-ranking."
 *   Iberoamerican Congress on Pattern Recognition (CIARP’2010), pages 541–548, 2010.
 *   http://dl.acm.org/citation.cfm?id=1948207.1948291
 *
 ******************************************************************************************
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

#include "Contextrr.hpp"

/* Constructor */
Contextrr::Contextrr() {

}

void Contextrr::loadParameters() {
    Exec exec = Exec::getInstance();
    exec.getConfigVariable(l,                "PARAM_CONTEXTRR_L");
    exec.getConfigVariable(k,                "PARAM_CONTEXTRR_K");
    exec.getConfigVariable(t,                "PARAM_CONTEXTRR_T");
    exec.getConfigVariable(nByK,             "PARAM_CONTEXTRR_NBYK");
    exec.getConfigVariable(hasOptimizations, "PARAM_CONTEXTRR_OPTIMIZATIONS");
}

void Contextrr::checkParameters() {

}

void Contextrr::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    //default structures
    rkLists.resize(n*n);
    initMatrix(matrix);

    //set core values
    n2 = n*n;
    l2 = l*l;
    incSz = n*k*l*l*5;
    nullImg = (n*k*l*l*5) + 1;
    maxHip = sqrt(2*l2);

    //method specific structures
    initMatrix(matrixW);
    incX.resize(incSz);
    incY.resize(incSz);
    incValue.resize(incSz);

    std::cout << "Initialized successfully!" << std::endl;
}

void Contextrr::initDataStructuresFusion() {
    initDataStructuresUdl();
}

void Contextrr::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
    } else { //RK
        genDistMatrixFromRks();
    }
}

void Contextrr::prepareInputFusion() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "SIM") {
            convertSimToDistMatrix();
        }
        genRksFromDistMatrix();
    } else { //RK
        genDistMatrixFromRks();
    }
}

void Contextrr::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void Contextrr::runUdlMethod() {
    for (curIteration = 1; curIteration <= t; curIteration++) {
        std::cout << "\n\n\n + Processing Iteration T = " << curIteration << "\n";

        setMatrixW(1);

        contextImgProc();

        searchCurrentMaxDistance();

        if (hasOptimizations) {
            execComputeNewDistsNoMin();
        } else {
            execComputeNewDists();
        }

        if (!hasOptimizations) {
            normalizeMinDistances();
        }

        execSortRankedLists();
    }
}

void Contextrr::runFusionMethod() {
    for (curIteration = 1; curIteration <= t; curIteration++) {
        std::cout << "\n\n\n + Processing Iteration T = " << curIteration << "\n";

        setMatrixW(1);

        if (curIteration == 1) {
            for (std::string const& file : fusionFiles) {
                //read the current file
                gettimeofday(&startTimeToDecrement, NULL);
                    readInputFile(file);
                    prepareInputFusion(); //assure that the both the dist matrix and the ranked lists are filled
                totalTimeToDecrement = Time::addTime(startTimeToDecrement, totalTimeToDecrement);

                contextImgProc();
            }
        } else {
            contextImgProc();
        }

        searchCurrentMaxDistance();

        if (hasOptimizations) {
            execComputeNewDistsNoMin();
        } else {
            execComputeNewDists();
        }

        if (!hasOptimizations) {
            normalizeMinDistances();
        }

        execSortRankedLists();
    }
}

void Contextrr::setMatrixW(float value) {
    std::cout << "\n\t Setting MatrixW ... \n";
    for (long int i = 0; i < n2; i++) {
        matrixW[i] = value;
    }
}

void Contextrr::contextImgProc() {
    if (hasOptimizations) {
        execContextImgProcNoSyncW();
    } else {
        execContextImgProc();
        processMatrixWUpdates();
    }
}

void Contextrr::execContextImgProc() {
    std::cout << "\n\t Processing context images ... \n";
    for (int i = 0; i < n; i++) {
        for (int j = 0;j < k; j++) {
            kernelCtxImageProc(i, j);
        }
    }
}

void Contextrr::execContextImgProcNoSyncW() {
    std::cout << "\n\t Processing context images - NoSyncW ... \n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < k; j++) {
            kernelCtxImageProcNoSyncW(i, j);
        }
    }
}

void Contextrr::kernelCtxImageProc(int pN, int pK) {
    //Constants to facilitate the reference to some variables
    const int& N = n;
    const int& K = k;
    const int& L = l;

    //General Variables
    unsigned int ctxImg[l2];
    unsigned int imgX, imgY; //Indexes
    unsigned int threshold = 0;
    unsigned int imgK = rkLists[N * pN + (pK + 1)];
    unsigned int Li;

    //Reseting Increments
    unsigned int NULL_IMG = (N * K * L * L * 5) + 1;
    unsigned int startInc = pN * (K * (L * L * 5)) + pK * (L * L * 5);
    unsigned int endInc = startInc + (L * L * 5);
    for (int i = startInc; i < endInc; i++) {
        incX[i] = NULL_IMG;
    }

    //****************** Computing Threshold Value ********************
    float curValue = 0, maxValue = 0, avgValue = 0;
    unsigned int Nn = N*pN;
    unsigned int Nk = N*imgK;
    unsigned int count = 0;
    for (int j = 1; j < L; j++) {
        int i = 1;
        imgX = rkLists[Nn + i];
        imgY = rkLists[Nk + j];
        curValue = matrix[imgX * N + imgY];
        count++;

        if (maxValue < curValue) {
            maxValue = curValue;
        }

        avgValue += curValue;
    }

    avgValue = avgValue / ((float) count);
    threshold = (unsigned int) ((avgValue / maxValue) * 255);

    //****************  Thresholding  ********************
    for (int i = 0; i < L; i++) {
        Li = L*i;
        for (int j = 0; j < L; j++) {
            imgX = rkLists[Nn + i];
            imgY = rkLists[Nk + j];
            curValue = matrix[imgX * N + imgY];
            curValue = curValue / maxValue;
            ctxImg[Li + j] = (unsigned int) (curValue * 255);
            if (ctxImg[Li + j] > threshold) {
                ctxImg[Li + j] = white;
            } else {
                ctxImg[Li + j] = black;
            }
        }
    }

    //************************ Median Filter **********************************
    unsigned int halfh = (h - 1) / 2;
    unsigned int halfi = ((((h * h) - 1) / 2) + 1);
    unsigned int wBegin = (halfh);
    unsigned int wEnd = L - (halfh) - 1;
    unsigned int hBegin = (halfh);
    unsigned int hEnd = L - (halfh) - 1;
    unsigned int indexFilter = 0;
    unsigned int countFilter = 0;

    for (int i = wBegin; i <= wEnd; i++) {
        Li = L*i;
        for (int j = hBegin; j <= hEnd; j++) {
            indexFilter = 0;
            countFilter = 0;
            for (int m = (i - halfh); m <= (i + halfh); m++) {
                for (int l = (j - halfh); l <= (j + halfh); l++) {
                    countFilter = countFilter + ctxImg[m * L + l];
                    indexFilter++;
                }
            }
            if (countFilter >= halfi - 1) {
                ctxImg[Li + j] = white;
            } else {
                ctxImg[Li + j] = black;
            }
        }
    }

    //************************  Increments  **********************************
    float weigth = (K - pK);
    float distFromOrigin;
    float baseInc;
    float baseIncBy4;
    unsigned int countInc = pN * (K * (L * L * 5)) + pK * (L * L * 5);

    for (int i = 0; i < L; i++) {
        Li = L*i;
        for (int j = 0; j < L; j++) {
            if (ctxImg[Li + j] == black) {
                //Indexes
                imgX = rkLists[Nn + i];
                imgY = rkLists[Nk + j];
                distFromOrigin = sqrt(pow(i + 1, 2) + pow(j + 1, 2));
                baseInc = weigth * (maxHip / distFromOrigin);
                baseIncBy4 = baseInc / 4;

                //------ Increments --------
                // Increment X x Y
                incX[countInc] = imgX * N + imgY;
                //incY[countInc] = imgY;
                incValue[countInc] = baseInc;
                countInc++;

                // Increment N x X
                incX[countInc] = pN * N + imgX;
                //incY[countInc] = imgX;
                incValue[countInc] = baseIncBy4;
                countInc++;

                // Increment N x Y
                incX[countInc] = pN * N + imgY;
                //incY[countInc] = imgY;
                incValue[countInc] = baseIncBy4;
                countInc++;

                // Increment K x X
                incX[countInc] = imgK * N + imgX;
                //incY[countInc] = imgX;
                incValue[countInc] = baseIncBy4;
                countInc++;

                // Increment K x Y
                incX[countInc] = imgK * N + imgY;
                //incY[countInc] = imgY;
                incValue[countInc] = baseIncBy4;
                countInc++;
            }
        }
    }
}

void Contextrr::kernelCtxImageProcNoSyncW(int pN, int pK) {
    //Constants to facilitate the reference to some variables
    const int& N = n;
    const int& K = k;
    const int& L = l;

    //General Variables
    unsigned int ctxImg[l2];
    unsigned int imgX, imgY; //Indexes
    unsigned int threshold = 0;
    unsigned int imgK = rkLists[N * pN + (pK + 1)];
    unsigned int Li;

    //****************** Computing Threshold Value ********************
    float curValue = 0, maxValue = 0, avgValue = 0;
    unsigned int Nn = N*pN;
    unsigned int Nk = N*imgK;
    unsigned int count = 0;
    //for (int i=1;i<L;i++) {
    for (int j = 1; j < L; j++) {
        int i = 1;
        imgX = rkLists[Nn + i];
        imgY = rkLists[Nk + j];
        curValue = matrix[imgX * N + imgY];
        count++;

        if (maxValue < curValue) {
            maxValue = curValue;
        }
        avgValue += curValue;
    }
    avgValue = avgValue / ((float) count);
    threshold = (unsigned int) ((avgValue / maxValue) * 255);

    //****************  Thresholding  ********************
    for (int i = 0; i < L; i++) {
        Li = L*i;
        for (int j = 0; j < L; j++) {
            imgX = rkLists[Nn + i];
            imgY = rkLists[Nk + j];
            curValue = matrix[imgX * N + imgY];
            curValue = curValue / maxValue;
            ctxImg[Li + j] = (unsigned int) (curValue * 255);
            if (ctxImg[Li + j] > threshold) {
                ctxImg[Li + j] = white;
            } else {
                ctxImg[Li + j] = black;
            }
        }
    }

    //************************ Median Filter **********************************
    unsigned int halfh = (h - 1) / 2;
    unsigned int halfi = ((((h * h) - 1) / 2) + 1);
    unsigned int wBegin = (halfh);
    unsigned int wEnd = L - (halfh) - 1;
    unsigned int hBegin = (halfh);
    unsigned int hEnd = L - (halfh) - 1;
    unsigned int indexFilter = 0;
    unsigned int countFilter = 0;

    for (int i = wBegin; i <= wEnd; i++) {
        Li = L*i;
        for (int j = hBegin; j <= hEnd; j++) {
            indexFilter = 0;
            countFilter = 0;
            for (int m = (i - halfh); m <= (i + halfh); m++) {
                for (int l = (j - halfh); l <= (j + halfh); l++) {
                    countFilter = countFilter + ctxImg[m * L + l];
                    indexFilter++;
                }
            }
            if (countFilter >= halfi - 1) {
                ctxImg[Li + j] = white;
            } else {
                ctxImg[Li + j] = black;
            }
        }
    }

    //************************  Increments  **********************************
    float weigth = (K - pK);
    float distFromOrigin;
    float baseInc;
    float baseIncBy4;
    unsigned int countInc = pN * (K * (L * L * 5)) + pK * (L * L * 5);

    for (int i = 0; i < L; i++) {
        Li = L*i;
        for (int j = 0; j < L; j++) {
            if (ctxImg[Li + j] == black) {
                //Indexes
                imgX = rkLists[Nn + i];
                imgY = rkLists[Nk + j];
                distFromOrigin = sqrt(pow(i + 1, 2) + pow(j + 1, 2));
                baseInc = weigth * (maxHip / distFromOrigin);
                baseIncBy4 = baseInc / 4;

                // Increment X x Y
                matrixW[imgX * N + imgY] += baseInc;

                // Increment N x X
                matrixW[pN * N + imgX] += baseIncBy4;

                // Increment N x Y
                matrixW[pN * N + imgY] += baseIncBy4;

                // Increment K x X
                matrixW[imgK * N + imgX] += baseIncBy4;

                // Increment K x Y
                matrixW[imgK * N + imgY] += baseIncBy4;

            }
        }
    }
}

void Contextrr::processMatrixWUpdates() {
    unsigned int index;
    for (int i = 0; i < incSz; i++) {
        if (incX[i] != nullImg) {
            matrixW[incX[i]] += incValue[i];
        }
    }
}

void Contextrr::searchCurrentMaxDistance() {
    std::cout << "\n\t Search current max distance ... \n";
    if (curIteration == 1) {
        maxDist = 0;
        for (long int i = 0; i < n2; i++) {
            maxDist = std::max(matrix[i], maxDist);
        }
    } else {
        maxDist = 2;
    }
}

void Contextrr::execComputeNewDists() {
    std::cout << "\n\t Computing new distances ... \n";
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            kernelComputeNewDistances(i, j);
        }
    }
}

void Contextrr::execComputeNewDistsNoMin() {
    std::cout << "\n\t Computing new distances - No Min ... \n";
    for (int i = 0; i < n; i++) {
        kernelComputeNewDistancesNoMin(i);
    }
}

void Contextrr::kernelComputeNewDistancesNoMin(int posX) {
    //Retrieving Dim. IDs
    int NPosX = cN*posX;

    for (int posY = posX + 1; posY < cN; posY++) {
        //Computing new distance value
        unsigned int indexXY = NPosX + posY;
        unsigned int indexYX = cN * posY + posX;

        //X x Y
        if (matrixW[indexXY] == 1) {
            matrix[indexXY] = 1 + (matrix[indexXY] / maxDist);
        } else {
            matrix[indexXY] = (1 / matrixW[indexXY]) * 2;
        }
        //Y x X
        if (matrixW[indexYX] == 1) {
            matrix[indexYX] = 1 + (matrix[indexYX] / maxDist);
        } else {
            matrix[indexYX] = (1 / matrixW[indexYX]) * 2;
        }

        matrix[indexXY] = std::min(matrix[indexXY], matrix[indexYX]);
        matrix[indexYX] = matrix[indexXY];
    }
    matrix[NPosX + posX] = 0;
}

void Contextrr::kernelComputeNewDistances(int posX, int posY) {
    //Computing new distance value
    unsigned int index = n*posX + posY;

    if (posX == posY) {
        matrix[index] = 0;
    } else {
        if (matrixW[index] == 1) {
            matrix[index] = 1 + (matrix[index] / maxDist) ;
        } else {
            matrix[index] = (1/matrixW[index]) * 2;
        }
    }
}

void Contextrr::normalizeMinDistances() {
    float normDist;
    unsigned int Ni, indexIJ, indexJI;
    for (int i = 0; i < n; i++) {
        Ni = n * i;
        for (int j = 0; j < n; j++) {
            indexIJ = Ni + j;
            indexJI = n*j + i;
            normDist = std::min(matrix[indexIJ], matrix[indexJI]);
            matrix[indexIJ] = normDist;
            matrix[indexJI] = normDist;
        }
    }
}

void Contextrr::execSortRankedLists() {
    std::cout << "\n\t Sorting ranked lists ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists(i);
    }
}

void Contextrr::kernelSortRankedLists(int curRL) {
    int cNcurRL = cN*curRL;
    float a[cN];

    for (int j = 0; j < cN; j++) {
        a[j] = matrix[cNcurRL + rkLists[cNcurRL + j]];
    }

    //---------------------- INSERTION SORT --------------------------
    int i, j, keyR;
    float keyA;

    for (j = 1; j < cN; j++) {
        keyA = a[j];
        keyR = rkLists[cNcurRL + j];
        i = j - 1;
        while (i >= 0 && a[i] > keyA) {
            a[i + 1] = a[i];
            rkLists[cNcurRL + i + 1] = rkLists[cNcurRL + i];
            i--;
        }
        a[i + 1] = keyA;
        rkLists[cNcurRL + i + 1] = keyR;
    }
    //----------------------------------------------------------------

    //Setting query image at first position
    i = 0;
    //while ((rk[i]!=curRL)&&(i<cN)) {
    while ((rkLists[cNcurRL + i] != curRL)&&(i < cN)) {
        i++;
    }
    if (i > 0) {
        int aux = rkLists[cNcurRL + 0];
        rkLists[cNcurRL + 0] = rkLists[cNcurRL + i];
        rkLists[cNcurRL + i] = aux;
    }
}
