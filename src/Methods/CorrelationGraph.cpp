/* <CorrelationGraph.cpp>
 *
 * CorrelationGraph method implementation file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 *************************************************************************************************
 *
 * Correlation Graph Manifold Learning is presented in the paper:
 *   D. C. G. Pedronette and R. da S. Torres.
 *   "A correlation graph approach for unsupervised manifold learning in image retrieval tasks."
 *   Neurocomputing, 208:66 – 79, 2016. SI: BridgingSemantic.
 *   http://dx.doi.org/10.1016/j.neucom.2016.03.081
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
#include <algorithm>

#include "CorrelationGraph.hpp"

/* Constructor */
CorrelationGraph::CorrelationGraph() {

}

void CorrelationGraph::loadParameters() {
    Exec exec = Exec::getInstance();

    //K
    exec.getConfigVariable(curK, "PARAM_CORGRAPH_K");

    //Threshold
    exec.getConfigVariable(thStart, "PARAM_CORGRAPH_THRESHOLD_START");
    exec.getConfigVariable(thEnd,   "PARAM_CORGRAPH_THRESHOLD_END");
    exec.getConfigVariable(thInc,   "PARAM_CORGRAPH_THRESHOLD_INC");

    //Correlation Function
    exec.getConfigVariable(correlationFunc, "PARAM_CORGRAPH_CORRELATION");

    //Rk size
    exec.getConfigVariable(l, "PARAM_CORGRAPH_L");
}

void CorrelationGraph::checkParameters() {
    if (curK > l) {
        std::cout << "K can't be greater than L" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
    if (l > n) {
        std::cout << "L can't be greater than N" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
    if ((thStart < 0) || (thEnd <= 0) || (thInc <= 0)) {
        std::cout << "The thresholds can't be negative" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
}

void CorrelationGraph::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    //default structures
    rkLists.resize(n*l);
    initMatrix(matrix);

    //algorithm structures
    correlationList.clear();
    correlationList.resize(n);

    //graph
    graph.clear();
    graph.resize(n);

    //tarjan
    index = 0;
    indexes.clear();
    lowlink.clear();
    stack.clear();
    onStack.clear();
    scc.clear();

    std::cout << "Initialized successfully!" << std::endl;
}

void CorrelationGraph::initDataStructuresFusion() {
    initMatrix(tmpMatrix);

    initDataStructuresUdl();
}

void CorrelationGraph::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
            genDistMatrixFromRks();
        }
    } else { //RK
        genDistMatrixFromRks();
    }
}

void CorrelationGraph::prepareInputFusion() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "SIM") {
            convertSimToDistMatrix();
        }
    } else { //RK
        genDistMatrixFromRks();
    }
}

void CorrelationGraph::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void CorrelationGraph::runUdlMethod() {
    runCorrelationGraphReRanking();
    computeFinalDistances();
    execSortRankedLists();
}

void CorrelationGraph::runFusionMethod() {
    prepareRankAggregationMatrix();
    execSortRankedLists();

    runUdlMethod();
}

void CorrelationGraph::prepareRankAggregationMatrix() {
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

        normalizeProbDB();

        //Multiply and store values
        for (long int i = 0; i < n*n; i++) {
            tmpMatrix[i] = tmpMatrix[i] * (1 + matrix[i]);
        }
    }

    //Redirect Matrix
    delete [] matrix;
    matrix = tmpMatrix;
    tmpMatrix = NULL;
}

void CorrelationGraph::normalizeProbDB() {
    //apply expoent
    for (int i = 0; i < n; i++) {
        double acumDistK = 0;
        for (int j = 0; j < curK; j++) {
            acumDistK += matrix[n*i + j];
        }
        double acumDistL = acumDistK;
        for (int j = curK; j < l; j++) {
            acumDistL += matrix[n*i + j];
        }
        double qest = 1.0 - (acumDistK/acumDistL);
        for (int j = 0; j < n; j++) {
            matrix[n*i + j] = pow(matrix[n*i + j], qest);
        }
    }

    //set reciprocal distances as the max
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            double dist1 = matrix[n*i + j];
            double dist2 = matrix[n*j + i];
            double newDist = std::max(dist1, dist2);
            matrix[n*i + j] = newDist;
            matrix[n*j + i] = newDist;
        }
    }
}

void CorrelationGraph::runCorrelationGraphReRanking() {
    std::cout << "\n\n >  Starting CorrelationGraph Re-Ranking ... \n\n";

    std::cout << ">>:  Computing correlation graph for  K = " << curK << "\n";
    cleanCorrelationValues();
    computeCorrelationForKNN();
    initSparseMatrix(matrix); //clean matrix values
    curThreshold = thStart;
    while (curThreshold < thEnd) {
        std::cout << "\n\t>>:  Computing correlation graph for  K = " << curK << ", THRESHOLD = " << curThreshold << "\n";
        cleanGraphStructures();
        buildCorrelationGraph();
        computeIncrementsByGraph();
        buildSCC();
        computeIncrementsBySCC();
        cleanGraphStructures();
        curThreshold += thInc;
    }
    cleanCorrelationValues();
}

void CorrelationGraph::computeFinalDistances() {
    //Normalization
    for (int i = 0; i < n; i++) {
        double acum =0;
        for (int j = 0; j < n; j++) {
            acum += matrix[j*n + i];
        }
        for (int j = 0; j < n; j++) {
            matrix[j*n + i] = matrix[j*n + i]/acum;
        }
    }
    //Compute Distances from Similarities
    for (long int i = 0; i < n*n; i++) {
        matrix[i] = 1/(1+matrix[i]);
    }
}

void CorrelationGraph::cleanCorrelationValues() {
    for (int i = 0; i < n; i++) {
        correlationList[i].clear();
        correlationList[i].resize(l);
    }
}

void CorrelationGraph::computeCorrelationForKNN() {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            int img = rkLists[i*l + j];
            if (correlationFunc == "RBO") {
                correlationList[i][j] = rbo(i, img, curK);
            } else if (correlationFunc == "PEARSON") {
                correlationList[i][j] = pearson(i, img, curK);
            }
        }
    }
}

void CorrelationGraph::cleanGraphStructures() {
    for (int i = 0; i < n; i++) {
        graph[i].clear();
    }
}

void CorrelationGraph::buildCorrelationGraph() {
    std::cout << "\n \t # Building Correlation Graph... \n";
    initializeGraph();
    for (int i = 0; i < n; i++) {
        for (int j = 1; j < l; j++) {
            if (correlationList[i][j] >= curThreshold) {
                //define edge
                Edge edge;
                edge.destination = rkLists[i*l + j];
                edge.weight = correlationList[i][j];
                //add edge to node
                graph[i].push_back(edge);
            }
        }
    }
}

void CorrelationGraph::initializeGraph() {
    for (int i = 0; i < n; i++) {
        //define edge
        Edge edge;
        edge.destination = i;
        edge.weight = correlationList[i][0]; //max correlation value
        //add edge to node
        graph[i].push_back(edge);
    }
}

void CorrelationGraph::computeIncrementsByGraph() {
    std::cout << "\n\t # Computing Increments by Graph... \n";
    double incValue = curThreshold;
    for (int i = 0; i < n; i++) { //for each node of the graph
        for (Edge edge1 : graph[i]) { //for each edge
            incDistance(i, edge1.destination, incValue);
            for (Edge edge2 : graph[i]) { //for each edge
                incDistance(edge1.destination, edge2.destination, incValue);
            }
        }
    }
}

void CorrelationGraph::computeIncrementsBySCC() {
    std::cout << "\n\t # Computing Increments by SCCs... \n";
    double incValue = curThreshold;
    for (std::vector<int> currentSCC : scc) { //for each SCC
        for (int node1 : currentSCC) {
            for (int node2 : currentSCC) {
                incDistance(node1, node2, incValue);
            }
        }
    }
}

void CorrelationGraph::buildSCC() {
    tarjan();
}

void CorrelationGraph::initTarjanStructures() {
    //clean
    indexes.clear();
    lowlink.clear();
    stack.clear();
    onStack.clear();
    scc.clear();

    //resize
    indexes.resize(n);
    lowlink.resize(n);
    onStack.resize(n);

    //fill
    for (int i = 0; i < n; i++) {
        indexes[i] = -1;
        lowlink[i] = 0;
        onStack[i] = false;
    }

    //init index
    index = 0;
}

void CorrelationGraph::tarjan() {
    initTarjanStructures();

    //visit all graph nodes
    for (int i = 0; i < n; i++) {
        if (indexes[i] == -1) {
            strongConnect(i);
        }
    }
}

void CorrelationGraph::strongConnect(int v) {
    //set the depth index for v to the smallest unused index
    indexes[v] = index;
    lowlink[v] = index;
    index += 1;
    stack.push_back(v);
    onStack[v] = true;

    //consider successors of v
    for (Edge e : graph[v]) {
        int w = e.destination;
        if (indexes[w] == -1) {
            //successor w has not yet been visited; recurse on it
            strongConnect(w);
            lowlink[v] = std::min(lowlink[v], lowlink[w]);
        } else if (onStack[w]) {
            //successor w is in stack S and hence in the current SCC
            lowlink[v] = std::min(lowlink[v], lowlink[w]);
        }
    }

    //if v is a root node, pop the stack and generate an SCC
    if (lowlink[v] == indexes[v]) {
        std::vector<int> sc; //start a new strongly connected component
        int w;
        do {
            w = stack.back(); //get element at the top
            stack.pop_back(); //pop (remove) element
            onStack[w] = false;
            sc.push_back(w); //add w to current strongly connected component
        } while (w != v);
        scc.push_back(sc); //output the current strongly connected component
    }
}

float CorrelationGraph::rbo(int i1, int i2, int paramK) {
    float p = 0.95;
    float inter = 0;
    int k = 1;
    float result = 0;

    while (k <= paramK) {
        inter = 0;
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                if (rkLists[l*i1 + i] == rkLists[l*i2 + j]) {
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

float CorrelationGraph::pearson(int i1, int i2, int paramK) {
    std::set<int> rkUnion;

    for (int j = 0; j < paramK; j++) {
        rkUnion.insert(rkLists[l*i1 + j]);
        rkUnion.insert(rkLists[l*i2 + j]);
    }

    double avg1 = 0;
    double avg2 = 0;
    int avgCount = 0;
    for (int img : rkUnion) {
        double img1Dist = matrix[n*i1 + img];
        double img2Dist = matrix[n*i2 + img];
        avg1 += img1Dist;
        avg2 += img2Dist;
        avgCount++;
    }
    avg1 = avg1 / ((double) avgCount);
    avg2 = avg2 / ((double) avgCount);

    double dividendo = 0;
    double divisor1 = 0;
    double divisor2 = 0;
    double divisor = 0;
    for (int img : rkUnion) {
        double img1Dist = matrix[n*i1 + img];
        double img2Dist = matrix[n*i2 + img];
        dividendo += ((img1Dist - avg1) * (img2Dist - avg2));
        divisor1 += pow((img1Dist - avg1), 2);
        divisor2 += pow((img2Dist - avg2), 2);
    }

    divisor = sqrt(divisor1) * sqrt(divisor2);
    double pearsonCorrelation = dividendo / divisor;
    double pearsonDist = (pearsonCorrelation + 1) / 2;

    return pearsonDist;
}

void CorrelationGraph::incDistance(int i, int j, double inc) {
    matrix[n*i + j] += inc;
}

void CorrelationGraph::execSortRankedLists() {
    std::cout << "\n\t Sort Ranked Lists ... \n";
    for (int i = 0; i < n; i++) {
        reRankingImageHeapSort(i);
    }
}

void CorrelationGraph::reRankingImageHeapSort(unsigned int rk) {
    std::vector<float> distances(n);
    std::vector<int> curRk(n);
    for (int j = 0; j < n; j++) {
        curRk[j] = j;
        distances[j] = matrix[n*rk + j];
    }
    heapsort(distances, curRk, n);
    int l = rkLists.size()/n;
    for (int j = 0; j < l; j++) {
        rkLists[l*rk + j] = curRk[j];
    }
}

void CorrelationGraph::heapsort(std::vector<float>& distances, std::vector<int>& curRk, int n) {
    buildheap(distances, curRk, n);
    while (n > 1) {
        n--;
        exchange(distances, curRk, 0, n);
        downheap(distances, curRk, n, 0);
    }
}

void CorrelationGraph::exchange(std::vector<float>& distances, std::vector<int>& curRk, int i, int j) {
    //Distances
    float t = distances[i];
    distances[i] = distances[j];
    distances[j] = t;
    //Ranked Lists
    int trk = curRk[i];
    curRk[i] = curRk[j];
    curRk[j] = trk;
}

void CorrelationGraph::downheap(std::vector<float>& distances, std::vector<int>& curRk, int n, int v) {
    int w = 2 * v + 1; //first descendant of v
    while (w < n) {
        if (w + 1 < n) {
            if (distances[w + 1] > distances[w]) {
                w++;
            }
        }
        if (distances[v] >= distances[w]) {
            return;
        }
        exchange(distances, curRk, v, w);
        v = w;
        w = 2 * v + 1;
    }
}

void CorrelationGraph::buildheap(std::vector<float>& distances, std::vector<int>& curRk, int n) {
    for (int v = n / 2 - 1; v >= 0; v--) {
        downheap(distances, curRk, n, v);
    }
}
