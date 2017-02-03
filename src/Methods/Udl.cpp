/* <Udl.cpp>
 *
 * Unsupervised Distance Learning implementation file
 *
 * @Authors: Lucas Pascotti Valem <lucasvalem@rc.unesp.br>
 *           Daniel Carlos Guimarães Pedronette <daniel@rc.unesp.br>
 *
 ***********************************************************************************
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
#include <fstream>

#include "Udl.hpp"

/* Constructor */
Udl::Udl() {
    Exec exec = Exec::getInstance();

    exec.getConfigVariable(udlTask, "UDL_TASK");

    if (udlTask == "FUSION") {
        exec.getInputFusionFiles(fusionFiles);
    }

    exec.getConfigVariable(inputFile,       "INPUT_FILE");
    exec.getConfigVariable(listsFile,       "INPUT_FILE_LIST");
    exec.getConfigVariable(classFile,       "INPUT_FILE_CLASSES");
    exec.getConfigVariable(inputFileFormat, "INPUT_FILE_FORMAT");
    exec.getConfigVariable(inputRkFormat,   "INPUT_RK_FORMAT");
    exec.getConfigVariable(inputMatrixType, "INPUT_MATRIX_TYPE");
    exec.getConfigVariable(rkFirstSorting,  "MATRIX_TO_RK_SORTING");

    exec.getConfigVariable(hasOutput,             "OUTPUT_FILE");
    exec.getConfigVariable(outputFileFormat,      "OUTPUT_FILE_FORMAT");
    exec.getConfigVariable(outputFile,            "OUTPUT_FILE_PATH");
    exec.getConfigVariable(outputRkFormat,        "OUTPUT_RK_FORMAT");
    exec.getConfigVariable(outputMatrixType,      "OUTPUT_MATRIX_TYPE");
    exec.getConfigVariable(showRkHtmlBeforeAfter, "OUTPUT_HTML_RK_BEFORE_AFTER");

    exec.getConfigVariable(efficiencyEval,      "EFFICIENCY_EVAL");
    exec.getConfigVariable(effectivenessEval,   "EFFECTIVENESS_EVAL");
    exec.getConfigVariable(computePrecisions,   "EFFECTIVENESS_COMPUTE_PRECISIONS");
    exec.getConfigVariable(computeMap,          "EFFECTIVENESS_COMPUTE_MAP");
    exec.getConfigVariable(computeRecall,       "EFFECTIVENESS_COMPUTE_RECALL");
    exec.getConfigVariable(recallAt,            "EFFECTIVENESS_RECALL_AT");
    exec.getConfigVariable(precisionsToCompute, "EFFECTIVENESS_PRECISIONS_TO_COMPUTE");

    exec.getConfigVariable(n, "SIZE_DATASET");
}

/* Runs the method implemented in the subclass and export the results if necessary */
void Udl::run() {
    loadParameters();
    checkParameters();

    if (udlTask == "UDL") {
        initDataStructuresUdl();
    } else { //FUSION
        initDataStructuresFusion();
        showRkHtmlBeforeAfter = false; //there's no original ranked list in this mode
    }

    readImagesList();

    Effectiveness effectiveness = Effectiveness(n, rkLists, imgList);
    if (effectivenessEval) {
        effectiveness.readClassesFile(classFile);
    }

    if (udlTask == "UDL") {
        readInputFile(inputFile);

        //check if ranked lists will be exported to html with before and after
        if ( (outputFileFormat == "RK") && (outputRkFormat == "HTML" || outputRkFormat == "ALL") && (showRkHtmlBeforeAfter == true) ) {
            std::cout << "\n Storing original ranked lists for HTML export ...\n\n\n";
            int l = rkLists.size()/n; //get the l value
            rkListsBefore.resize(n*l);
            if (inputFileFormat == "MATRIX") {
                if (inputMatrixType == "DIST") {
                    genRksFromDistMatrix();
                } else { //SIM
                    genRksFromSimMatrix();
                }
                rkListsBefore = rkLists; //store the original ranked list
            }
        }

        prepareInput();
    } //if it's fusion, we will do this inside the runMethod function

    if ((udlTask == "UDL") && (effectivenessEval)) { //this pre-evaluation is not necessary for fusion
        if (computePrecisions) {
            effectiveness.fillPrecisionsMap(precisionsBefore, precisionsToCompute);
        }
        if (computeMap) {
            mapBefore = effectiveness.computeMAPMeasure();
        }
        if (computeRecall) {
            recallBefore = effectiveness.computeRecall(recallAt);
        }
        std::cout << "\n";
    }

    if (udlTask == "UDL") {
        gettimeofday(&startTimeElapsed, NULL);
            runUdlMethod();
        totalTimeElapsed = Time::addTime(startTimeElapsed, totalTimeElapsed);
    } else { //FUSION
        gettimeofday(&startTimeElapsed, NULL);
            runFusionMethod();
        totalTimeElapsed = Time::addTime(startTimeElapsed, totalTimeElapsed);
        totalTimeElapsed -= totalTimeToDecrement;
    }

    if (hasOutput) {
        prepareOutput();
        writeOutput(effectiveness);
    }

    if (effectivenessEval) {
        if (computePrecisions) {
            effectiveness.fillPrecisionsMap(precisionsAfter, precisionsToCompute);
        }
        if (computeMap) {
            mapAfter = effectiveness.computeMAPMeasure();
        }
        if (computeRecall) {
            recallAfter = effectiveness.computeRecall(recallAt);
        }
        std::cout << "\n";
    }

    generateExecutionLog();
    std::cout << "\n***********************************************************************\n\n";
    TxtFile::printFile("log.txt");

    releaseDataStructures();
}

/* Generates a log a file that contains the evaluation results requested by the user */
void Udl::generateExecutionLog() {
    std::cout << "\n Writing execution log file (log.txt) ...\n";

    //open file
    std::ofstream file;
    file.open("log.txt");

    file << std::fixed << std::setprecision(4); //decimal precision for recall@, map, etc.

    //title
    file << " # UNSUPERVISED DISTANCE LEARNING METHODS FRAMEWORK #\n";

    //main information
    file << "\n\n - GENERAL INFORMATION -\n";
    file << " --------------------------------------";
    std::string methodName;
    Exec::getInstance().getConfigVariable(methodName, "UDL_METHOD");
    file << "\n Task:             " << udlTask;
    file << "\n Method:           " << methodName;
    file << "\n Dataset Size:     " << n;
    file << "\n Image List File:  " << listsFile;
    if (effectivenessEval) {
        file << "\n Image Class File: " << classFile;
    }
    if (udlTask == "UDL") {
        file << "\n Input File:       " << inputFile;
    } else { //FUSION
        file << "\n Input Files:";
        int i = 1;
        for (std::string const& filename : fusionFiles) {
            file << "\n                   (" << i << ") " << filename;
            i++;
        }
    }
    file << "\n Input Format:     " << inputFileFormat;
    if (inputFileFormat == "MATRIX") {
        file << " " << inputMatrixType;
    } else {
        file << " " << inputRkFormat;
    }
    if (hasOutput) {
        file << "\n Output File:      " << outputFile;
        file << "\n Output Format:    " << outputFileFormat;
        if (outputFileFormat == "MATRIX") {
            file << " " << outputMatrixType;
        } else {
            file << " " << outputRkFormat;
        }
    }
    file << "\n --------------------------------------\n";

    //parameters
    file << "\n\n - METHOD PARAMETERS -\n";
    file << " --------------------------------------";
    file << Exec::getInstance().getMethodParameterVariables();
    file << "\n --------------------------------------\n";

    //evaluation info (efficiency and effectiveness)
    if (effectivenessEval || efficiencyEval) {
        file << "\n\n - EVALUATION RESULTS -\n";
        file << " --------------------------------------";

        if (efficiencyEval) {
            file << "\n * Efficiency:";
            file << "\n\t Total Time of the Algorithm Execution: " << totalTimeElapsed << " s\n";
        }
        if ( (effectivenessEval) && (computePrecisions || computeRecall || computeMap) ) {
            file << "\n * Effectiveness:";
            if (udlTask == "UDL") { //not necessary for fusion
                file << "\n Before: ";
                if (computePrecisions) {
                    for (auto const& elem : precisionsBefore) {
                        file << "\n\t P@" << elem.first << "\t\t" << elem.second;
                    }
                }
                if (computeRecall) {
                    file << "\n\t Recall@" << recallAt << "\t" << recallBefore;
                }
                if (computeMap) {
                    file << "\n\t MAP\t\t" << mapBefore;
                }
                file << "\n After: ";
            }
            if (computePrecisions) {
                for (auto const& elem : precisionsAfter) {
                    file << "\n\t P@" << elem.first << "\t\t" << elem.second;
                }
            }
            if (computeRecall) {
                file << "\n\t Recall@" << recallAt << "\t" << recallAfter;
            }
            if (computeMap) {
                file << "\n\t MAP\t\t" << mapAfter;
            }
            if (udlTask == "UDL") { //not necessary for fusion
                file << "\n Relative Gains: ";
                if (computePrecisions) {
                    for (auto const& elem : precisionsAfter) {
                        int value = elem.first;
                        float after = elem.second;
                        float before = precisionsBefore[value];
                        file << "\n\t P@" << value << "\t\t" << (after-before)/before*100 << "\%";
                    }
                }
                if (computeRecall) {
                    file << "\n\t Recall@" << recallAt << "\t" << ((recallAfter-recallBefore)/recallBefore)*100 << "\%";
                }
                if (computeMap) {
                    file << "\n\t MAP\t\t" << ((mapAfter-mapBefore)/mapBefore)*100 << "\%";
                }
            }
        }

        file << "\n --------------------------------------\n";
    }

    file << "\n Log generated at " << Time::getCurrentTime(); //print current time at the end of the log

    file.close();
}

/* Initializes the distance matrix allocating memory immediately */
void Udl::initMatrix(float*& matrix) {
    std::cout << "\t - Intializing Distance Matrix ..." << std::endl;

    delete [] matrix;
    matrix = new float[n*n]();  //initialize all values to zero, but allocate memory immediately

    std::cout << "\t - Matrix Successfully Initialized ..." << std::endl;
}

/* Initializes a sparse sparse distance matrix */
void Udl::initSparseMatrix(float*& matrix) {
    std::cout << "\t - Intializing Sparse Distance Matrix ..." << std::endl;

    delete [] matrix;
    matrix = new float[n*n];

    for (long int l = 0; l < n*n; l++) {
        if (matrix[l] != 0) {
            matrix[l] = 0;
        }
    }

    std::cout << "\t - Matrix Successfully Initialized ..." << std::endl;
}

/* Releases data structures that aren't automatically released */
void Udl::releaseDataStructures() {
    delete [] matrix;
}

/* Reads the images list, line by line */
void Udl::readImagesList() {
    std::cout << "\n Starting reading images list ... [" << listsFile << "]  \n";

    imgList.clear();
    std::ifstream inFile;
    inFile.open(listsFile.c_str());
    if (!inFile) {
        std::cerr << " Unable to open image list file [" << listsFile << "].";
        exit(1); //terminate with error
    }
    std::string line;
    int i = 0;
    while (inFile >> line) {
        imgList.push_back(line);
        i++;
    }
    inFile.close();

    if (imgList.size() != n) {
        std::cout << " Your image list file is invalid for a dataset of " << n <<
        " images! There are " << imgList.size() << " images in it! "<< std::endl;
        exit(1); //terminate with error
    }

    std::cout << " Done! \n";
}

/* Reads the input file, choosing the right action (matrix, rk, ...) */
void Udl::readInputFile(std::string inputFile) {
    if (inputFileFormat == "RK") {
        if (inputRkFormat == "NUM") {
            readRkListsNumeric(inputFile);
        } else { //STR
            readRkListsStr(inputFile);
        }
    } else { //MATRIX
        readDistMatrix(inputFile);
    }
}

/* Reads a numeric ranked lists file */
void Udl::readRkListsNumeric(std::string inputFile) {
    std::cout << "\n Starting reading rkLists ... [" << inputFile.c_str() << "]  \n";

    std::ifstream inFile;
    inFile.open(inputFile.c_str());

    int l = rkLists.size()/n; //get the l value

    if (!inFile) {
        std::cerr << " Unable to open rkLists file [" << inputFile << "].";
        exit(1); //terminate with error
    }

    std::string line;
    int i, j;
    int curValue;
    i = 0;
    while (i < n) {
        j = 0;
        while (j < l) {
            inFile >> line;
            curValue = atoi(line.c_str());
            rkLists[l * i + j] = curValue;
            j++;
        }
        char next;
        while (inFile.get(next)) {
            if (next == '\n') {
                break;
            }
        }
        i++;
    }

    inFile.close();

    std::cout << " Done! \n\n";
}

/* Reads a ranked lists file that contains the images by its names */
void Udl::readRkListsStr(std::string inputFile) {
    std::cout << "\n Starting reading rkLists ... [" << inputFile.c_str() << "]  \n";

    std::ifstream inFile;
    inFile.open(inputFile.c_str());

    int l = rkLists.size()/n; //get the l value

    if (!inFile) {
        std::cerr << " Unable to open rkLists file [" << inputFile << "].";
        exit(1); //terminate with error
    }

    std::string line;
    int i, j;
    int curValue;
    i = 0;
    while (i < n) {
        j = 0;
        while (j < l) {
            inFile >> line;
            curValue = getImageNumber(line.c_str());
            if (curValue == -1) {
                std::cerr << " Unable to open rkLists file [" << inputFile << "].\n";
                std::cerr << " Invalid image string: " << line.c_str() << "\n";
                exit(1); //terminate with error
            }
            rkLists[l * i + j] = curValue;
            j++;
        }
        char next;
        while (inFile.get(next)) {
            if (next == '\n') {
                break;
            }
        }
        i++;
    }

    inFile.close();

    std::cout << " Done! \n\n";
}

/* Reads the distance matrix */
void Udl::readDistMatrix(std::string inputFile) {
    std::cout << "\n Starting reading matrix A ... [" << inputFile.c_str() << "]  \n";
    std::ifstream inFile;
    inFile.open(inputFile.c_str());
    if (!inFile) {
        std::cerr << "Unable to open distances file [" << inputFile.c_str() << "].";
        exit(1); //terminate with error
    }
    std::string line;
    int i, j, Ni;
    double curValue;
    i = 0;
    while (i < n) {
        j = 0;
        Ni = n*i;
        while ((j < n) && (inFile >> line)) {
            curValue = atof(line.c_str());
            matrix[Ni + j] = curValue;
            j++;
        }
        i++;
    }
    inFile.close();
    std::cout << " Done! \n\n";
}

/* Writes the output file, choosing the right action (matrix, rk, ...) */
void Udl::writeOutput(Effectiveness& effectiveness) {
    std::cout << "\n Exporting to output file(s) ... " << std::endl;

    if (outputFileFormat == "RK") {
        if (outputRkFormat == "NUM") {
            exportRkListsNumeric(outputFile);
        } else if (outputRkFormat == "STR") {
            exportRkListsStr(outputFile);
        } else if (outputRkFormat == "HTML") {
            exportRkListsHtml(effectiveness, outputFile);
        } else { //ALL
            exportRkListsNumeric(outputFile + "_num");
            exportRkListsStr(outputFile + "_str");
            exportRkListsHtml(effectiveness, outputFile);
        }
    } else { //MATRIX
        exportDistMatrix(outputFile);
    }

    std::cout << "\n\n Exported successfully!" << std::endl;
}

/* Exports the ranked list as a file of numbers (positions) */
void Udl::exportRkListsNumeric(std::string path) {
    path += ".txt";
    std::cout << "\n\t [Numeric RK Lists] Writing to output file: " << path;

    std::ofstream file;
    file.open(path);

    int l = rkLists.size()/n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            file << rkLists[l*i + j] << " ";
        }
        file << "\n";
    }

    file.close();
}

/* Exports the ranked list as a file of strings (generally, images file names) */
void Udl::exportRkListsStr(std::string path) {
    path += ".txt";
    std::cout << "\n\t [String RK Lists] Writing to output file: " << path;

    std::ofstream file;
    file.open(path);

    int l = rkLists.size()/n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            file << imgList[rkLists[l*i + j]] << " ";
        }
        file << "\n";
    }

    file.close();
}

/* Exports the ranked list to html files for a better visualization of the results */
void Udl::exportRkListsHtml(Effectiveness& effectiveness, std::string path) {
    //get config variables
    int rkPerFile, imagesPerRk;
    std::string imagesPath;
    bool htmlColors;
    Exec exec = Exec::getInstance();
    exec.getConfigVariable(rkPerFile,   "OUTPUT_HTML_RK_PER_FILE");
    exec.getConfigVariable(imagesPerRk, "OUTPUT_HTML_RK_SIZE");
    exec.getConfigVariable(imagesPath,  "INPUT_IMAGES_PATH");
    exec.getConfigVariable(htmlColors,  "OUTPUT_HTML_RK_COLORS");

    //disable colors if effectiveness evaluation is disabled
    if (!effectivenessEval) {
        htmlColors = false;
    }

    //verify config variables
    int l = rkLists.size()/n;
    if (rkPerFile == 0) {
        rkPerFile = n;
    }
    if ((imagesPerRk == 0) || (imagesPerRk > l)) {
        imagesPerRk = l;
    }

    int number = ceil((float) n/(float) rkPerFile); //number of files to generate
    int image = 0; //current image to generate the ranked list

    //print information in the terminal
    std::cout << "\n\t [HTML RK Lists] Writing html files ... ";
    std::cout << "\n\t\t Config: " << rkPerFile << " ranked lists per file, " << imagesPerRk << " images each";
    if (number == 1) {
        std::cout << "\n\t\t Writing html file...";
    } else {
        std::cout << "\n\t\t Writing html files. There are " << number << " files in total... please wait.";
    }

    //write each one of the html files
    for (int k = 0; k < number; k++) {
        //create file to write
        std::ostringstream filename;
        filename << path << std::setfill('0') << std::setw(Type::numDigits(number)) << k+1 << ".html";
        std::ofstream file;
        file.open(filename.str());

        //file header
        file << "<html>\n";
        file << "<meta charset=\"UTF-8\">\n";
        file << "<title>Ranked Lists</title>\n";
        file << "<body bgcolor=\"FFFFFF\">\n";
        file << "<center>";
        file << "<b>Ranked Lists</b>\n";
        file << "<br><b>File " << k+1 << "/" << number << "</b>\n";
        file << "<br><br>\n";

        //check if rkPerFile exceeds the size of the dataset
        int lim = (image+rkPerFile);
        if (lim > n) {
            lim = n;
        }

        //write ranked lists
        for (int i = image; i < lim; i++) {
            //initiate variables
            int imgI;
            std::string baseClass, curClass;

            //table header
            file << "<table border=\"3\"><tr><td>\n";

            //export ranked list (before - original)
            if (showRkHtmlBeforeAfter) {
                imgI = rkListsBefore[l*i];

                if (htmlColors) {
                    baseClass = effectiveness.getClass(imgI);
                }

                file << "<table cellspacing=\"2\" cellpadding=\"2\" border=\"0\" bordercolor=\"black\"  BGCOLOR=\"#000000\"><tr>\n";
                file << std::setfill('0') << std::setw(Type::numDigits(n)) << i+1 << "\n";
                if (htmlColors) {
                    file << "<td bgcolor=\"#00FF00\"><img src=\"file://" << imagesPath << imgList[imgI] << "\" width=\"80\" height=\"80\"></td>\n";
                } else {
                    file << "<td><img src=\"file://" << imagesPath << imgList[imgI] << "\" width=\"80\" height=\"80\"></td>\n";
                }
                for (int j = 1; j < imagesPerRk; j++) {
                    int imgJ = rkListsBefore[l*i + j];
                    if (htmlColors) {
                        curClass = effectiveness.getClass(imgJ);
                    }
                    if (curClass == baseClass) {
                        file << "<td><img src=\"file://" << imagesPath << imgList[imgJ] << "\" width=\"80\" height=\"80\"></td>\n";
                    } else {
                        if (htmlColors) {
                            file << "<td bgcolor=\"#FF0000\"><img src=\"file://" << imagesPath << imgList[imgJ] << "\" width=\"80\" height=\"80\"></td>\n";
                        } else {
                            file << "<td><img src=\"file://" << imagesPath << imgList[imgJ] << "\" width=\"80\" height=\"80\"></td>\n";
                        }
                    }
                }
                file << "</tr></table>\n";
            }

            //export ranked list (after)
            imgI = rkLists[l*i];
            if (htmlColors) {
                baseClass = effectiveness.getClass(imgI);
            }
            file << "<table cellspacing=\"2\" cellpadding=\"2\" border=\"0\" bordercolor=\"black\"  BGCOLOR=\"#000000\"><tr>\n";
            file << "<br>";
            if (htmlColors) {
                file << "<td bgcolor=\"#00FF00\"><img src=\"file://" << imagesPath << imgList[imgI] << "\" width=\"80\" height=\"80\"></td>\n";
            } else {
                file << "<td><img src=\"file://" << imagesPath << imgList[imgI] << "\" width=\"80\" height=\"80\"></td>\n";
            }
            for (int j = 1; j < imagesPerRk; j++) {
                int imgJ = rkLists[l*i + j];
                if (htmlColors) {
                    curClass = effectiveness.getClass(imgJ);
                }
                if (curClass == baseClass) {
                    file << "<td><img src=\"file://" << imagesPath << imgList[imgJ] << "\" width=\"80\" height=\"80\"></td>\n";
                } else {
                    if (htmlColors) {
                        file << "<td bgcolor=\"#FF0000\"><img src=\"file://" << imagesPath << imgList[imgJ] << "\" width=\"80\" height=\"80\"></td>\n";
                    } else {
                        file << "<td><img src=\"file://" << imagesPath << imgList[imgJ] << "\" width=\"80\" height=\"80\"></td>\n";
                    }
                }
            }
            file << "</tr></table>\n";
            file << "</td></tr></table>\n\n";

            file << "<br>\n";
        }
        file << "</center>";
        file << "</body>\n";
        file << "</html>\n";

        file.close();

        image += rkPerFile;
    }
}

/* Exports the results to a distance matrix */
void Udl::exportDistMatrix(std::string path) {
    path += ".txt";
    std::cout << "\n [Distance Matrix] Writing to output file: " << path << std::endl;

    //export distance matrix to file
    std::ofstream file;
    file.open(path);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            file << matrix[n*i + j] << " ";
        }
        file << "\n";
    }

    file.close();
}

/* Reinitializes and fills the distance matrix based in the ranked lists */
void Udl::genDistMatrixFromRks() {
    //reinitialize the matrix
    initMatrix(matrix);

    //given the ranked lists, fill the matrix with the positions
    int l = rkLists.size()/n;
    for (long int i = 0; i < n*n; i++) {
        matrix[i] = l;
    }
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            int img = rkLists[l*i + j];
            matrix[n*i + img] += j - l;
        }
    }
}

/* Reinitializes and fills the similarity matrix based in the ranked lists */
void Udl::genSimMatrixFromRks() {
    //reinitialize the matrix
    initMatrix(matrix);

    //given the ranked lists, fill the matrix with the positions
    int l = rkLists.size()/n;
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < l; j++) {
            int img = rkLists[l*i + j];
            matrix[n*i + img] = l - j;
        }
    }
}

/* Given a string, returns the index of the image, returns -1 if not found */
int Udl::getImageNumber(std::string image) {
    for (int i = 0; i < imgList.size(); i++) {
        if (imgList[i] == image) {
            return i;
        }
    }

    return -1;
}

/* Fills the ranked lists based in the distance matrix */
void Udl::genRksFromDistMatrix() {
    if (rkFirstSorting == "HEAP") {
        mainHeapSort("DIST");
    } else { //INSERTION
        insertionSortDist();
    }
}

/* Fills the ranked lists based in the similarity matrix */
void Udl::genRksFromSimMatrix() {
    if (rkFirstSorting == "HEAP") {
        mainHeapSort("SIM");
    } else { //INSERTION
        insertionSortSim();
    }
}

/* Converts similarity matrix to distance matrix */
void Udl::convertSimToDistMatrix() {
    for (long int i = 0; i < n*n; i++) {
        matrix[i] = 1/(1 + matrix[i]);
    }
}

/* Converts distance matrix to similarity matrix */
void Udl::convertDistToSimMatrix() {
    for (long int i = 0; i < n*n; i++) {
        matrix[i] = 1/(1 + matrix[i]);
    }
}

// ---------------- SORTING METHODS ------------------- //
// Sorting methods to obtain ranked lists from distance/similarity matrices


// ---------------- INSERTION SORT ------------------- //

void Udl::insertionSortDist() {
    int l = rkLists.size()/n;

    std::vector<int> curRk(n); //temporary ranked list for current iteration

    for (int rk = 0; rk < n; rk++) {
        curRk.clear(); //clean current ranked list

        int cNcurRL = n*rk;
        double a[n];

        //-------  Creating auxiliar structures  --------
        for (int j = 0; j < n; j++) {
            curRk[j] = j;
        }

        for (int j = 0; j < n; j++) {
            a[j] = matrix[cNcurRL + curRk[j]];
        }

        //---------------------- INSERTION SORT --------------------------
        int i, j, keyR;
        double keyA;

        for (j = 1; j < n; j++) {
            keyA = a[j];
            keyR = curRk[j];
            i = j - 1;
            while (i >= 0 && a[i] > keyA) {
                a[i + 1] = a[i];
                curRk[i + 1] = curRk[i];
                i--;
            }
            a[i + 1] = keyA;
            curRk[i + 1] = keyR;
        }
        //----------------------------------------------------------------

        //Setting query image at first position
        i = 0;
        while ((curRk[i] != rk) && (i < n)) {
            i++;
        }
        if (i > 0) {
            int aux = curRk[0];
            curRk[0] = curRk[i];
            curRk[i] = aux;
        }

        //Copy current ranked list to ranked lists vector
        for (int j = 0; j < l; j++) {
            rkLists[rk*l + j] = curRk[j];
        }
    }
}

void Udl::insertionSortSim() {
    int l = rkLists.size()/n;

    std::vector<int> curRk(n); //temporary ranked list for current iteration

    for (int rk = 0; rk < n; rk++) {
        curRk.clear(); //clean current ranked list

        int cNcurRL = n*rk;
        double a[n];

        //-------  Creating auxiliar structures  --------
        for (int j = 0; j < n; j++) {
            curRk[j] = j;
        }

        for (int j = 0; j < n; j++) {
            a[j] = matrix[cNcurRL + curRk[j]];
        }

        //---------------------- INSERTION SORT --------------------------
        int i, j, keyR;
        double keyA;

        for (j = 1; j < n; j++) {
            keyA = a[j];
            keyR = curRk[j];
            i = j - 1;
            while (i >= 0 && a[i] < keyA) {
                a[i + 1] = a[i];
                curRk[i + 1] = curRk[i];
                i--;
            }
            a[i + 1] = keyA;
            curRk[i + 1] = keyR;
        }
        //----------------------------------------------------------------

        //Setting query image at first position
        i = 0;
        while ((curRk[i] != rk) && (i < n)) {
            i++;
        }
        if (i > 0) {
            int aux = curRk[0];
            curRk[0] = curRk[i];
            curRk[i] = aux;
        }

        //Copy current ranked list to ranked lists vector
        for (int j = 0; j < l; j++) {
            rkLists[rk*l + j] = curRk[j];
        }
    }
}


// ---------------- HEAP SORT ------------------- //

void Udl::mainHeapSort(std::string type) {
    for (int rk = 0; rk < n; rk++) {
        std::vector<float> distances(n);
        std::vector<int> curRk(n);
        for (int j = 0; j < n; j++) {
            curRk[j] = j;
            distances[j] = matrix[n*rk + j];
        }
        heapsort(distances, curRk, n, type);
        int l = rkLists.size()/n;
        for (int j = 0; j < l; j++) {
            rkLists[l*rk + j] = curRk[j];
        }
    }
}

void Udl::heapsort(std::vector<float>& distances, std::vector<int>& curRk, int n, std::string type) {
    buildheap(distances, curRk, n, type);
    while (n > 1) {
        n--;
        exchange(distances, curRk, 0, n);
        if (type == "SIM") {
            downheapSim(distances, curRk, n, 0);
        } else { //DIST
            downheapDist(distances, curRk, n, 0);
        }
    }
}

void Udl::exchange(std::vector<float>& distances, std::vector<int>& curRk, int i, int j) {
    //Distances
    float t = distances[i];
    distances[i] = distances[j];
    distances[j] = t;
    //Ranked Lists
    int trk = curRk[i];
    curRk[i] = curRk[j];
    curRk[j] = trk;
}

void Udl::downheapDist(std::vector<float>& distances, std::vector<int>& curRk, int n, int v) {
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

void Udl::downheapSim(std::vector<float>& distances, std::vector<int>& curRk, int n, int v) {
    int w = 2 * v + 1; //first descendant of v
    while (w < n) {
        if (w + 1 < n) {
            if (distances[w + 1] < distances[w]) {
                w++;
            }
        }
        if (distances[v] <= distances[w]) {
            return;
        }
        exchange(distances, curRk, v, w);
        v = w;
        w = 2 * v + 1;
    }
}

void Udl::buildheap(std::vector<float>& distances, std::vector<int>& curRk, int n, std::string type) {
    for (int v = n / 2 - 1; v >= 0; v--) {
        if (type == "SIM") {
            downheapSim(distances, curRk, n, v);
        } else { //DIST
            downheapDist(distances, curRk, n, v);
        }
    }
}
