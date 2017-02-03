/* <RlRecom.cpp>
 *
 * RL-Recommendation method implementation file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 **************************************************************************************************
 *
 * RL-Recommendation is presented in the paper:
 *   L. P. Valem, D. C. G. Pedronette, R. d. S. Torres, E. Borin, and J. Almeida.
 *   "Effective, efficient, and scalable unsupervised distance learning in image retrieval tasks."
 *   ICMR, 2015
 *   http://dx.doi.org/10.1145/2671188.2749336
 *
 **************************************************************************************************
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
#include <algorithm>

#include "RlRecom.hpp"

/* Constructor */
RlRecom::RlRecom() {

}

void RlRecom::loadParameters() {
    Exec exec = Exec::getInstance();
    exec.getConfigVariable(l,       "PARAM_RLRECOM_L");
    exec.getConfigVariable(k,       "PARAM_RLRECOM_K");
    exec.getConfigVariable(lambda,  "PARAM_RLRECOM_LAMBDA");
    exec.getConfigVariable(epsilon, "PARAM_RLRECOM_EPSILON");
}

void RlRecom::checkParameters() {
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

void RlRecom::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*l);
    initSparseMatrix(matrix);

    cohesionVector.resize(n);

    std::cout << "Initialized successfully!" << std::endl;
}

void RlRecom::initDataStructuresFusion() {
    std::cerr << "WARNING: The Fusion Mode is not implemented for this method! Aborting ...\n";
    exit(1);
}

void RlRecom::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
        initSparseMatrix(matrix);
    }
}

void RlRecom::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void RlRecom::runUdlMethod() {
    int iteration = 1;
    float curCohesion = 0;
    float lastCohesion = curCohesion;

    execFillMatrix();
    execFillPosMatrix();

    do {

        std::cout << std::endl << "  -> Iteration " << iteration << std::endl;

        execCalcCohesion();

        lastCohesion = curCohesion;
        for (int i = 0; i < n; i++) {
            curCohesion += cohesionVector[i];
        }
        curCohesion = curCohesion/((float)(n));

        execPerformRecommendations();
        execSortRankedLists();

        k++;
        iteration++;

        std::cout << "\n\t Current Cohesion: " << curCohesion << std::endl;

    } while ( (curCohesion*epsilon) <= (curCohesion-lastCohesion) );

}

void RlRecom::runFusionMethod() {
    std::cerr << "WARNING: The Fusion Mode is not implemented for this method! Aborting ...\n";
    exit(1);
}

void RlRecom::execFillMatrix() {
    std::cout << "\n\t Fill Distance Matrix ... \n";
    for (int i = 0; i < n; i++) {
        kernelFillMatrix(i);
    }
}

void RlRecom::execFillPosMatrix() {
    std::cout << "\n\t Fill Distance Matrix with Positions ... \n";
    for (int i = 0; i < n; i++) {
        kernelFillPosMatrix(i);
    }
}

void RlRecom::execCalcCohesion() {
    std::cout << "\n\t Calculate Cohesions ... \n";
    for (int i = 0; i < n; i++) {
        kernelCalcCohesion(i);
    }
}

void RlRecom::execPerformRecommendations() {
    std::cout << "\n\t Perform Recommendations ... \n";
    for (int i = 0; i < n; i++) {
        kernelPerformRecommendations(i);
    }
}

void RlRecom::execSortRankedLists() {
    std::cout << "\n\t Sort Ranked Lists ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists(i);
    }
}

void RlRecom::kernelFillMatrix(int rk) {
    for (int pos = 0; pos < l; pos++) {
        int img = rkLists[l*rk + pos];
        matrix[n*rk + img] = 2*l;
        matrix[n*img + rk] = 2*l;
    }
}

void RlRecom::kernelFillPosMatrix(int rk) {
    for (int pos = 0; pos < l; pos++) {
        int img = rkLists[l*rk + pos];
        matrix[n*rk + img] = matrix[n*rk + img] + pos - l;
        matrix[n*img + rk] = matrix[n*img + rk] + pos - l;
    }
}

void RlRecom::kernelCalcCohesion(int rk) {
    if (k <= 1) {
        cohesionVector[rk] = 0;
        return;
    }

    float cohesion = 0;
    int j;
    for (int i = 1; i <= k; i++) {
        j = 1;
        while (j + 1 <= k) {
            int image = rkLists[l*rkLists[rk*l + i] + j];
            if (isIn(rk, image, k)) {
                cohesion += weightNormCohesion(j);
            }
            j++;
        }
    }

    cohesionVector[rk] = ( cohesion / getAcumWeightCohesion(k) );
}

void RlRecom::kernelPerformRecommendations(int rk) {
    float cohesion = cohesionVector[rk];
    for (int i = 0; i <= k; i++) {
        int posI = i + 1;
        for (int j = 0; j <= k; j++) {
            int posJ = j + 1;
            float wi = 1 - ((float) posI) / ((float) k + 1);
            float wj = 1 - ((float) posJ) / ((float) k + 1);
            float weight = (wi * wj) * (cohesion);
            float multBy = (1 - std::min(1.0, lambda * weight));
            float currentDistance = (getCurrentDistance(rkLists[rk * l + i], rkLists[rk * l + j]));
            float newDistance = currentDistance * multBy;
            float currentInvDist = (getCurrentDistance(rkLists[rk * l + j], rkLists[rk * l + i]));
            setNewDistance(rkLists[rk * l + i], rkLists[rk * l + j], std::min(newDistance, currentInvDist));
        }
    }
}

void RlRecom::kernelSortRankedLists(int rk) {
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

float RlRecom::weightNormCohesion(int position) {
    return 1 / ((float) position+1);
}

float RlRecom::getAcumWeightCohesion(int size) {
    float acumWeight = 0;
    for (int i = 2; i <= size; i++) {
        acumWeight += weightNormCohesion((i-1));
    }
    return (((float)size) * acumWeight);
}

bool RlRecom::isIn(int rk, int image, int size) {
    bool contains = false;
    int i = 1;
    while ( (!contains) and (i <= size) ) {
        if (rkLists[rk*l + i] == image) {
            contains = true;
        }
        i++;
    }

    return contains;
}

float RlRecom::getCurrentDistance(int i, int j) {
    return matrix[i*n + j];
}

void RlRecom::setNewDistance(int i, int j, float newDist) {
    matrix[i*n + j] = newDist;
}
