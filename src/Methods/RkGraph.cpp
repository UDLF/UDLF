/* <RkGraph.cpp>
 *
 * RkGraph method implementation file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *************************************************************************************************
 *
 * Ranked List Graph Distance is presented in the paper:
 *   D. C. G. Pedronette, J. Almeida, and R. da S. Torres.
 *   "A graph-based ranked-list model for unsupervised distance learning on shape retrieval."
 *   Pattern Recognition Letters, 83, Part 3:357 – 367, 2016. Efficient Shape Representation,
 *   Matching, Ranking, and its Applications.
 *   http://dx.doi.org/10.1016/j.patrec.2016.05.021
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

#include "RkGraph.hpp"

/* Constructor */
RkGraph::RkGraph() {

}

void RkGraph::loadParameters() {
    Exec exec = Exec::getInstance();
    exec.getConfigVariable(k, "PARAM_RKGRAPH_K");
    exec.getConfigVariable(t, "PARAM_RKGRAPH_T");
    exec.getConfigVariable(p, "PARAM_RKGRAPH_P");
    exec.getConfigVariable(l, "PARAM_RKGRAPH_L");
}

void RkGraph::checkParameters() {
    if (l > n) {
        std::cout << "L can't be greater than N" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
}

void RkGraph::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*l);
    initMatrix(matrix);

    std::cout << "Initialized successfully!" << std::endl;
}

void RkGraph::initDataStructuresFusion() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*l);
    initMatrix(matrix);
    initMatrix(matrixAgg);
    for (long int i = 0; i < n*n; i++) {
        matrixAgg[i] = 1;
    }

    std::cout << "Initialized successfully!" << std::endl;
}

void RkGraph::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
    }
}

void RkGraph::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void RkGraph::runUdlMethod() {
    std::cout << "\n # Executing RkGraph Algorithm #\n\n";
    computeMutualRankDists();
    execSortRankedLists();
    for (int ite = 0; ite < t; ite++) {
        std::cout << "\n\n============== ITERATION " << ite << " ===============\n";
        runIteration();
    }
}

void RkGraph::runFusionMethod() {
    std::cout << "\n # Executing RkGraph Algorithm (Rank Aggregation) #\n\n";

    for (std::string const& file : fusionFiles) {
        //Read File
        gettimeofday(&startTimeToDecrement, NULL);
            readInputFile(file);
            prepareInput();
        totalTimeToDecrement = Time::addTime(startTimeToDecrement, totalTimeToDecrement);

        //Perform Algorithm
        computeMutualRankDists();
        execSortRankedLists();
        runIteration();

        //Store values in aggregation matrix (multiplying)
        for (long int i = 0; i < n*n; i++) {
            matrixAgg[i] *= (1 + matrix[i]);
        }
    }
    delete [] matrix;   //the main values are not necessary anymore
    matrix = matrixAgg; //the current matrix is now the one that has the accumulated values

    //Update rankings
    execSortRankedLists();

    //Perform final iterations
    for (int ite = 0; ite < t; ite++) {
        std::cout << "\n\n============== ITERATION (RA): " << ite << " ===============\n";
        runIteration();
    }
}

void RkGraph::initStructure(double*& structure) {
    delete [] structure;
    structure = new double[n*n](); //initialize all values to zero, but allocate memory immediately
}

void RkGraph::computeMutualRankDists() {
    std::cout << "\n\t Rank Normalization ... \n";

    for (long int i = 0; i < n*n; i++) {
        matrix[i] = 2*l;
    }

    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            int imgj = rkLists[l*i + j];
            matrix[n*i + imgj] = getPosition(i, imgj) + getPosition(imgj, i);
        }
    }
}

int RkGraph::getPosition(int qImg, int img) {
    for (int i = 0; i < l; i++) {
        if (rkLists[l*img + i] == qImg) {
            return (i+1);
        }
    }
    return (l+1);
}

void RkGraph::initGraphStructures() {
    initStructure(smetric);
    initStructure(adj);
}

void RkGraph::runIteration() {
    initGraphStructures();

    std::cout << "\n\n (*) Computing Similarity Metric ...";
    for (int i = 0; i < n; i++) {
        fillSimMetricForImage(i);
    }

    std::cout << "\n\n (*) Incrementing Adjacencies ...";
    for (int i = 0; i < n; i++) {
        computeIncAdjForImage(i);
    }

    std::cout << "\n\n (*) Computing New Distances ...";
    computeNewDists();

    std::cout << "\n\n (*) Updating Models ...";
    execSortRankedLists();
}

void RkGraph::fillSimMetricForImage(int img) {
    for (int i = 0; i <= k; i++) {
        int imgi = rkLists[l*img + i];
        for (int j = 0; j <= k; j++) {
            int imgj = rkLists[l*img + j];
            setSimMetric(imgi, imgj);
        }
    }
}

void RkGraph::setSimMetric(int img1, int img2) {
    if (smetric[n*img1 + img2] == 0) {
        double sim = calcRBO(img1, img2);
        smetric[n*img1 + img2] = sim;
        smetric[n*img2 + img1] = sim;
    }
}

double RkGraph::calcRBO(int i1, int i2) {
    double inter = 0;
    double result = 0;
    int curK = 1;

    while (curK <= k) {
        inter = 0;
        for (int i = 0; i < curK; i++) {
            for (int j = 0; j < curK; j++) {
                if (rkLists[l*i1 + i] == rkLists[l*i2 + j]) {
                    inter += 1;
                }
            }
        }
        result += pow(p, curK - 1) * (inter / curK);
        curK++;
    }

    result = (1 - p) * result;

    return result;
}

void RkGraph::computeIncAdjForImage(int imgq) {
    for (int i = 0; i < k; i++) {
        int imgi = rkLists[l*imgq + i];
        double valueSim = smetric[n*imgq + imgi];
        valueSim = valueSim*valueSim*2;
        adj[n*imgq + imgi] += valueSim;
        adj[n*imgi + imgq] += valueSim;
    }

    for (int i = 0; i <= k; i++) {
        int imgi = rkLists[l*imgq + i];
        for (int j = 0; j <= k; j++) {
            int imgj = rkLists[l*imgq + j];
            double simqi = smetric[n*imgq + imgi];
            double simqj = smetric[n*imgq + imgj];
            double valueSim = (simqi*simqj)*2;
            adj[n*imgi + imgj] += valueSim;
        }
    }
}

void RkGraph::computeNewDists() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            matrix[n*i + j] = 1.0/(1 + adj[n*i + j]);
        }
    }
}

/*
void RkGraph::computeNewDists() {
    for (int i = 0; i < n; i++) {
        int qimg = i;
        for (int j = 0; j < l; j++) {
            int img = rkLists[l*qimg + j];
            double valueAdj = adj[n*qimg + img];
            double newDist;
            if (valueAdj > 0) {
                newDist = 1.0/valueAdj;
            } else {
                newDist = i+1;
            }
            matrix[n*qimg + img] = newDist;
        }
    }
}
*/

void RkGraph::execSortRankedLists() {
    std::cout << "\n\t Sort Ranked Lists ... \n";
    for (int i = 0; i < n; i++) {
        kernelSortRankedLists(i);
    }
}

void RkGraph::kernelSortRankedLists(int rk) {
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
