/* <Udl.hpp>
 *
 * Unsupervised Distance Learning class header file
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

#ifndef UDL_HPP
#define UDL_HPP

#include <math.h>

#include <memory>
#include <map>
#include <vector>
#include <array>
#include <sstream>
#include <iomanip>

#include "Core/Exec.hpp"
#include "Evaluation/Effectiveness.hpp"
#include "Utils/Type.hpp"
#include "Utils/Time.hpp"
#include "Utils/TxtFile.hpp"

/* Abstract class which contains functions and attributes that are common to all UDL methods */
class Udl {
    public:
            Udl();
            void run();

    protected:
            //general functions that ALL METHODS must implement
            virtual void runUdlMethod() = 0;
            virtual void runFusionMethod() = 0;
            virtual void initDataStructuresUdl() = 0;
            virtual void initDataStructuresFusion() = 0;
            virtual void loadParameters() = 0;
            virtual void checkParameters() = 0;
            virtual void prepareInput() = 0;
            virtual void prepareOutput() = 0;

            //init/release structures
            void initMatrix(float*& matrix);
            void initSparseMatrix(float*& matrix);
            void releaseDataStructures();

            //reading files
            void readInputFile(std::string inputFile);

            //convert matrices to ranked lists
            void genDistMatrixFromRks();
            void genSimMatrixFromRks();
            int getImageNumber(std::string image);

            //convert ranked lists to matrices
            void genRksFromDistMatrix();
            void genRksFromSimMatrix();

            //convert matrix type
            void convertSimToDistMatrix();
            void convertDistToSimMatrix();

            //general structures for all methods
            int n; //dataset size
            std::vector<std::string> imgList;
            std::vector<int> rkLists;
            float* matrix = NULL;

            //time evaluation
            timeval startTimeToDecrement;
            float totalTimeToDecrement = 0;

            //values read from config
            std::string udlTask;
            std::string inputFile;
            std::vector<std::string> fusionFiles;
            std::string listsFile;
            std::string classFile;
            std::string inputFileFormat;
            std::string inputRkFormat;
            std::string inputMatrixType;
            std::string outputFile;
            std::string outputFileFormat;
            std::string outputRkFormat;
            std::string outputMatrixType;
            bool hasOutput;
            bool effectivenessEval;
            bool efficiencyEval;
            bool computePrecisions;
            bool computeRecall;
            bool computeMap;
            int recallAt;
            std::string precisionsToCompute;

    private:
            //reading files
            void readImagesList();
            void readDistMatrix(std::string inputFile);
            void readRkListsNumeric(std::string inputFile);
            void readRkListsStr(std::string inputFile);

            //writing files
            void writeOutput(Effectiveness& effectiveness);
            void exportDistMatrix(std::string path);
            void exportRkListsNumeric(std::string path);
            void exportRkListsStr(std::string path);
            void exportRkListsHtml(Effectiveness& effectiveness, std::string path);

            //generate log
            void generateExecutionLog();

            //sorting methods to obtain ranked lists from distance/similarity matrices
            //insertion
            void insertionSortDist();
            void insertionSortSim();
            //heapsort
            void mainHeapSort(std::string type);
            void heapsort(std::vector<float>& distances, std::vector<int>& curRk, int n, std::string type);
            void exchange(std::vector<float>& distances, std::vector<int>& curRk, int i, int j);
            void downheapDist(std::vector<float>& distances, std::vector<int>& curRk, int n, int v);
            void downheapSim(std::vector<float>& distances, std::vector<int>& curRk, int n, int v);
            void buildheap(std::vector<float>& distances, std::vector<int>& curRk, int n, std::string type);

            //rk html export
            bool showRkHtmlBeforeAfter;
            std::vector<int> rkListsBefore; //store the original input ranked lists if necessary

            //rk convertion
            std::string rkFirstSorting;

            //time evaluation
            timeval startTimeElapsed;
            float totalTimeElapsed = 0;

            //effectiveness evalution variables
            std::map<int, float> precisionsBefore;
            std::map<int, float> precisionsAfter;
            float mapBefore;
            float mapAfter;
            float recallBefore;
            float recallAfter;

};

#endif // UDL_HPP
