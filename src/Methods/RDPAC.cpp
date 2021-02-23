/* <RDPAC.cpp>
 *
 * RDPAC implementation file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *******************************************************************************************************************
 *
 * Rank-based Diffusion Process with Assured Convergence (RDPAC) is presented in the paper:
 *    PEDRONETTE, D. C. G.; VALEM, L. P.; LATECKI, L. J. .
 *    Efficient Rank-based Diffusion Process with Assured Convergence.
 *    Journal of Imaging (ISSN 2313-433X).
 *    (Under minor revision)
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

#include "RDPAC.hpp"

/* Constructor */
RDPAC::RDPAC() {

}

void RDPAC::loadParameters() {
    Exec exec = Exec::getInstance();

    //RDPAC Parameters (read from the configuration file)
    exec.getConfigVariable(l, "PARAM_RDPAC_L");
    exec.getConfigVariable(l_mult, "PARAM_RDPAC_L_MULT");
    exec.getConfigVariable(p, "PARAM_RDPAC_P");
    exec.getConfigVariable(pl, "PARAM_RDPAC_PL");
    exec.getConfigVariable(k_start, "PARAM_RDPAC_K_START");
    exec.getConfigVariable(k_inc, "PARAM_RDPAC_K_INC");
    exec.getConfigVariable(k_end, "PARAM_RDPAC_K_END");
}

void RDPAC::checkParameters() {
    if (l > n) {
        std::cout << "L can't be greater than N" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
}

void RDPAC::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    if (l_mult*l > n) {
        std::cout << "WARNING: Lx" << l_mult;
        std::cout << " > N. Setting L = N/" << l_mult << "\n";
        l = ((int) n/l_mult);
    }

    rkLists.resize(l_mult*l*n);
    initSparseMatrix(matrix);
    W = matrix;
    initSparseMatrix(RWL);
    tmpLine = (float*) malloc(n*sizeof(float));

    std::cout << "Initialized successfully!" << std::endl;
}

void RDPAC::initDataStructuresFusion() {
    initDataStructuresUdl();

    initSparseMatrix(fuseW);
    setnz.clear();
    for (int i = 0; i < n; i++) {
        setnz[i] = std::vector<int>();
    }
}

void RDPAC::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
    }
}

void RDPAC::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

inline int RDPAC::getRKElem(long int query, long int pos) {
    long int index = ((long int) ((l_mult*l)*query + pos));
    return rkLists[index];
}

inline float RDPAC::getMatrixElem(float* m, long int i, long int j) {
    long int index = ((long int) n*i + j);
    return m[index];
}

inline void RDPAC::setMatrixElem(float* m, long int i, long int j, float value) {
    long int index = ((long int) n*i + j);
    m[index] = value;
}

void RDPAC::generateFullLMatrixW() {
    std::cout << "\t - Generate Full L Matrix W\n";
    for (int imgI = 0; imgI < n; imgI++) {
        setMatrixElem(W, imgI, imgI, 1.0);
        for (int j = 0; j < l; j++) {
            int imgJ = getRKElem(imgI, j);
            setMatrixElem(W, imgI, imgJ, pow(pl, (j+1)));
        }
    }
}

void RDPAC::copyMatrices(float* src, float* dst, long int lim) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < lim; j++) {
            int img = getRKElem(i, j);
            setMatrixElem(dst, i, img, getMatrixElem(src, i, img));
        }
    }
}

void RDPAC::sumTranspostL() {
    float* tmpMatrix = NULL;
    initSparseMatrix(tmpMatrix);

    for (int img1 = 0; img1 < n; img1++) {
        for (int pos = 0; pos < l_mult*l; pos++) {
            int img2 = getRKElem(img1, pos);
            float value = getMatrixElem(W, img1, img2) + getMatrixElem(W, img2, img1);
            setMatrixElem(tmpMatrix, img1, img2, value);
        }
    }

    copyMatrices(tmpMatrix, W, (long int) l_mult*l);
    delete [] tmpMatrix; // dispose temporary matrix
}


void RDPAC::copyingMatrixWToDBL() {
    std::cout << "\t - Copying Matrix W to DBL ...\n";
    // iterate over the ranked lists of all images
    // in the top-2l postions to fill the remaining
    // value in the distance matrix
    for (int i = 0; i < n; i++) {
        for (int pos = 0; pos < l_mult*l; pos++) {
            int imgJ = getRKElem(i, pos);
            float newDist =  getMatrixElem(W, i, imgJ);
            if (newDist < 0) {
                setMatrixElem(matrix, i, imgJ, n-pos);
            }
        }
    }
}

void RDPAC::reciprocalReferences() {
    std::cout << "\t - Performing Reciprocal References\n";

    generateFullLMatrixW();
    sumTranspostL();
    copyingMatrixWToDBL();
}

void RDPAC::execFillPosMatrix() {
    std::cout << "\n\t Fill Distance Matrix with Positions ... \n";
    for (int i = 0; i < n; i++) {
        kernelFillPosMatrix(i);
    }
    copyingMatrixWToDBL();
}

void RDPAC::kernelFillPosMatrix(int rk) {
    for (int pos = 0; pos < l; pos++) {
        long int img = getRKElem(rk, pos);
        float value = pow(pl, pos+1);
        setMatrixElem(matrix, rk, img, getMatrixElem(matrix, rk, img) + value);
        setMatrixElem(matrix, img, rk, getMatrixElem(matrix, img, rk) + value);
    }
}

void RDPAC::execSortRankedLists() {
    std::cout << "\n\t Sort Ranked Lists ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists(i);
    }
}

void RDPAC::kernelSortRankedLists(int rk) {
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

void RDPAC::generatekNNMatrixW(int k) {
    std::cout << "\t - Generate kNN Matrix W ...\n";

    cleanSimMatrix();

    for (int imgI = 0; imgI < n; imgI++) {
        setMatrixElem(W, imgI, imgI, 1.0);
        for (int j = 0; j < k; j++) {
            int imgJ = getRKElem(imgI, j);
            setMatrixElem(W, imgI, imgJ, pow(p, (j+1)));
        }
    }
}

void RDPAC::normalizeMatrix(float* m, int lim) {
    std::cout << "\t - Normalizing Matrix with lim=" << lim <<  " ...\n";

    float* vAcum = (float*) malloc(n*sizeof(float));

    for (int i = 0; i < n; i++) {
        vAcum[i] = 0;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < lim; j++) {
            int jidx = getRKElem(i, j);
            vAcum[jidx] += getMatrixElem(m, i, jidx);
        }
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < lim; j++) {
            int jidx = getRKElem(i, j);
            float value = getMatrixElem(m, i, jidx) / vAcum[jidx];
            setMatrixElem(m, i, jidx, value);
        }
    }

    delete [] vAcum; // dispose array of acummulated values
}

void RDPAC::multiplyByTranspostLStoreList(int k) {
    for (int i = 0; i < n; i++) {
        std::vector<int> listNonZero;
        for (int j = 0; j < l; j++) {
            int jidx = getRKElem(i, j);
            tmpLine[jidx] = 0;
            for (int l = 0; l < k; l++) {
                int lidx = getRKElem(jidx, l);
                tmpLine[jidx] += getMatrixElem(RWL, i, lidx) * getMatrixElem(W, jidx, lidx);
            }
            if (tmpLine[jidx] > 0) {
                listNonZero.push_back(jidx);
            }
        }
        nonZero[i] = listNonZero;
        for (int j = 0; j < l; j++) {
            int jidx = getRKElem(i, j);
            setMatrixElem(RWL, i, jidx, tmpLine[jidx]);
        }
    }
}

void RDPAC::multiplyByTranspostL(int k) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            int jidx = getRKElem(i, j);
            tmpLine[jidx] = 0;
            for (int l = 0; l < k; l++) {
                int lidx = getRKElem(jidx, l);
                tmpLine[jidx] += getMatrixElem(RWL, i, lidx) * getMatrixElem(W, jidx, lidx);
            }
        }
        for (int j = 0; j < l; j++) {
            int jidx = getRKElem(i, j);
            setMatrixElem(RWL, i, jidx, tmpLine[jidx]);
        }
    }
}

void RDPAC::sumIdentityAlpha(float alpha) {
    std::cout << "\t - Adding Identity Alpha ...\n";

    for (int i = 0; i < n; i++) {
        for (int jc = 0; jc < l_mult*l; jc++) {
            int j = getRKElem(i, jc);
            if (i!=j) {
                float value = getMatrixElem(RWL, i, j)*alpha;
                setMatrixElem(RWL, i, j, value);
            } else {
                float value = getMatrixElem(RWL, i, j)*alpha + (1.0 * (1.0-alpha));
                setMatrixElem(RWL, i, j, value);
            }
        }
    }
}

void RDPAC::multiplyL1(float* m1, float* m2, float* dst) {
    float* tmpColumn = tmpLine;
    for (int jc = 0; jc < l_mult*l; jc++) {
        for (int i = 0; i < n; i++) {
            int j = getRKElem(i, jc);
            tmpColumn[i] = 0;
            for (int l : nonZero[i]) {
                tmpColumn[i] += getMatrixElem(m1, i, l) * getMatrixElem(m2, l, j);
            }
        }
        for (int i = 0; i < n; i++) {
            int j = getRKElem(i, jc);
            setMatrixElem(dst, i, j, tmpColumn[i]);
        }
    }
}

void RDPAC::multiplyL2(float* m1, float* m2, float* dst) {
    float* tmpMatrix = NULL;
    initSparseMatrix(tmpMatrix);

    for (int i = 0; i < n; i++) {
        std::vector<int> listNonZero = nonZero[i];
        for (int jc = 0; jc < l_mult*l; jc++) { // change to 2l later
            int j = getRKElem(i, jc);
            float value = 0;
            for (int l : listNonZero) {
                value += getMatrixElem(m1, i, l) * getMatrixElem(m2, l, j);
            }
            setMatrixElem(tmpMatrix, i, j, value);
        }
    }

    copyMatrices(tmpMatrix, dst, (long int) l_mult*l);
    delete [] tmpMatrix;
}

void RDPAC::cleanSimMatrix() {
    delete [] matrix;
    matrix = NULL;
    initSparseMatrix(matrix);
    W = matrix;
}

void RDPAC::cleanRWLMatrix() {
    delete [] RWL;
    RWL = NULL;
    initSparseMatrix(RWL);
}

void RDPAC::performSingleExecution() {
    execFillPosMatrix();
    execSortRankedLists();

    for (int k_iter = k_start; k_iter <= k_end; k_iter += k_inc) {
        int curK = k_end;

        std::cout << "\n\t* Running iteration for K=" << curK << "\n";

        generatekNNMatrixW(curK);

        if (k_iter == k_start) { // first execution
            cleanRWLMatrix();
            copyMatrices(W, RWL, curK); // RWL = W
        }

        normalizeMatrix(RWL, l);
        normalizeMatrix(W, curK);

        if (k_iter == k_end) { // last execution (P = P*Wt)
            multiplyByTranspostLStoreList(curK);
        } else {
            multiplyByTranspostL(curK);
        }

        sumIdentityAlpha(0.95);

        std::cout << "\n\n";
    }
    std::cout << "\n\t* Running final steps\n";
    normalizeMatrix(RWL, l);

    multiplyL1(RWL, RWL, W); // W = RWL x RWL
    multiplyL1(RWL, W, W); // W = RWL x W
}

void RDPAC::runUdlMethod() {
    std::cout << "\n Executing RDPAC!\n\n";

    performSingleExecution();

    copyingMatrixWToDBL();
    execSortRankedLists();
}

void RDPAC::copyMatrixWtoFuseWL() {
    for (int i = 0; i < n; i++) {
        int imgi = i;
        for (int jc = 0; jc < l_mult*l; jc++) {
            int j = getRKElem(i, jc);
            float value = getMatrixElem(fuseW, i, j) + (1.0 + getMatrixElem(W, i, j));
            setMatrixElem(fuseW, i, j, value);
            int imgj = j;
            setnz[imgi].push_back(imgj);
        }
    }
}

void RDPAC::copyMatrixFuseWtoDBL() {
    std::cout << "\t\t L-Copying matrix fuseW to DB...\n";
    for (int i = 0; i < n; i++) {
        int imgi = i;
        for (int imgj : setnz[imgi]) {
            int j = imgj;
            float newSim = getMatrixElem(fuseW, i, j);
            setMatrixElem(matrix, imgi, imgj, newSim);
        }
    }
    execSortRankedLists();
}

void RDPAC::runFusionMethod() {
    std::cout << "\n * Running Fusion for RDPAC \n";

    for (std::string curFile : fusionFiles) {
        std::cout << "\n Processing file " << curFile << "\n";
        cleanSimMatrix();
        readInputFile(curFile);
        prepareInput();
        execSortRankedLists();
        performSingleExecution();
        generateFullLMatrixW();
        copyMatrixWtoFuseWL();
    }
    cleanSimMatrix();

    copyMatrixFuseWtoDBL();

    performSingleExecution();

    copyingMatrixWToDBL();
    execSortRankedLists();
}
