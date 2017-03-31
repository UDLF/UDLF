/* <Cprr.cpp>
 *
 * CPRR method implementation file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *******************************************************************************************************************
 *
 * Cartesian Product of Ranking References (CPRR) is presented in the paper:
 *   L. P. Valem and D. C. G. Pedronette.
 *   "Unsupervised similarity learning through cartesian product of ranking references for image retrieval tasks."
 *   2016 29th SIBGRAPI Conference on Graphics, Patterns and Images (SIBGRAPI).
 *   http://dx.doi.org/10.1109/SIBGRAPI.2016.042
 *   http://sibgrapi.sid.inpe.br/col/sid.inpe.br/sibgrapi/2016/07.22.13.56/doc/PaperSIBGRAPI-2016_vFinal.pdf
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

#include "Cprr.hpp"

/* Constructor */
Cprr::Cprr() {

}

void Cprr::loadParameters() {
    Exec exec = Exec::getInstance();
    exec.getConfigVariable(l, "PARAM_CPRR_L");
    exec.getConfigVariable(k, "PARAM_CPRR_K");
    exec.getConfigVariable(t, "PARAM_CPRR_T");
}

void Cprr::checkParameters() {
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

void Cprr::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*l);
    initSparseMatrix(matrix);

    imgRef.resize(n*l);
    posRef.resize(n*l);
    nRef.resize(n);

    std::cout << "Initialized successfully!" << std::endl;
}

void Cprr::initDataStructuresFusion() {
    std::cout << "Initializing data structures..." << std::endl;

    initDataStructuresUdl(); //init the main structures

    //structures used for aggregation (fusion) tasks
    initSparseMatrix(matrixAgg);
    rkListsAgg.resize(n);  //make sure it's empty

    std::cout << "Initialized successfully!" << std::endl;
}

void Cprr::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
        initSparseMatrix(matrix);
    }
}

void Cprr::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void Cprr::runUdlMethod() {
    std::cout << std::endl << "  -> Noniterative Steps" << std::endl;
    execFillPosMatrix();
    execSortRankedListsZero();

    int iteration = 1;
    while (iteration <= t) {
        std::cout << std::endl << "  -> Iteration " << iteration << std::endl;
        execCartProd();
        execReverseCartProd();
        execCleanRef();
        if (iteration != t) {
            execSortRankedListsZero();
        } else {
            execSortRankedLists();
        }
        iteration++;
    }
}

void Cprr::runFusionMethod() {
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
        runUdlMethod();

        //store values in a secondary matrix
        for (long int i = 0; i < n*n; i++) {
            matrixAgg[i] += matrix[i];
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
    int steps = 2; //this is a fixed value for the final iterations
    int iteration = 1;
    while (iteration <= steps) {
        std::cout << std::endl << "  -> Step " << iteration << std::endl;
        execCartProd();
        execReverseCartProd();
        execCleanRef();
        if (iteration != steps) {
            execSortRankedListsZero();
        } else {
            execSortRankedLists();
        }
        iteration++;
    }
}

void Cprr::execFillPosMatrix() {
    std::cout << "\n\t Fill Distance Matrix with Positions ... \n";
    for (int i = 0; i < n; i++) {
        kernelFillPosMatrix(i);
    }
}

void Cprr::execSortRankedLists() {
    std::cout << "\n\t Sort Ranked Lists ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists(i);
    }
}

void Cprr::execSortRankedListsZero() {
    std::cout << "\n\t Sort Ranked Lists (Zero) ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedListsZero(i);
    }
}

void Cprr::execSortRankedListsAgg(std::vector<std::vector<int>>& rk) {
    std::cout << "\n\t Sort Ranked Lists (Fusion) ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedListsAgg(i, rk);
    }
}

void Cprr::execCartProd() {
    std::cout << "\n\t Cartesian Product ... \n";
    for (int i = 0; i < n; i++) {
        kernelCartProd(i);
    }
}

void Cprr::execReverseCartProd() {
    std::cout << "\n\t Reverse Cartesian Product ... \n";
    for (int i = 0; i < n; i++) {
        kernelReverseCartProd(i);
    }
}

void Cprr::execCleanRef() {
    std::cout << "\n\t Clean Reference List ... \n";
    for (int i = 0; i < n; i++) {
        kernelCleanRef(i);
    }
}

void Cprr::kernelFillPosMatrix(int rk) {
    long int lrk = l*rk;
    long int nrk = n*rk;
    long int img;

    for (int pos = 0; pos < l; pos++) {
        img = rkLists[lrk + pos];
        matrix[nrk + img] += l - pos;
        matrix[n * img + rk] += l - pos;
    }
}

void Cprr::kernelSortRankedLists(int rk) {
    long int LcurRL = l*rk;
    long int cNcurRL = n*rk;
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

void Cprr::kernelSortRankedListsZero(int rk) {
    long int LcurRL = l*rk;
    long int cNcurRL = n*rk;
    float a[l];

    long int index;
    for (int j = 0; j < l; j++) {
        index = cNcurRL + rkLists[LcurRL + j];
        a[j] = matrix[index];
        matrix[index] = 0;
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

void Cprr::kernelSortRankedListsAgg(int curRL, std::vector<std::vector<int>>& rk) {
    int LcurRL = rk[curRL].size();
    int cNcurRL = n*curRL;
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

void Cprr::kernelCartProd(int rk) {
    int set[k];
    set[0] = rk;

    long int Li = l*rk;
    long int imgNei, index;

    for (long int j = 1; j < k; j++) {
        imgNei = rkLists[Li + j];
        set[j] = imgNei;
        index = l * imgNei + nRef[imgNei];
        imgRef[index] = rk;
        posRef[index] = k - j;
        nRef[imgNei]++;
    }

    long int idx, idy, Nidx;
    double score;

    for (long int j = 0; j < k; j++) {
        idx = set[j];
        Nidx = n*idx;
        for (long int l = 0; l < k; l++) {
            idy = set[l];
            score = ((k - (j + 1) + 1)*(k - (l + 1) + 1));
            matrix[Nidx + idy] += score;
            matrix[n * idy + idx] += score;
        }
    }

}

void Cprr::kernelReverseCartProd(int rk) {
    long int size = nRef[rk];
    long int Li = l*rk;
    long int idx, idy, Nidx, index, index2;
    double score;

    for (long int j = 0; j < size; j++) {
        index = Li + j;
        idx = imgRef[index];
        Nidx = n*idx;
        for (long int l = 0; l < size; l++) {
            index2 = Li + l;
            idy = imgRef[index2];
            score = (posRef[index])*(posRef[index2]);
            matrix[Nidx + idy] += score;
            matrix[n * idy + idx] += score;
        }
    }
}

void Cprr::kernelCleanRef(int rk) {
    nRef[rk] = 0;
}

bool Cprr::search(int img, const std::vector<int>& rk) {
    for (int i = 0; i < rk.size(); i++) {
        if (img == rk[i]) {
            return true;
        }
    }
    return false;
}
