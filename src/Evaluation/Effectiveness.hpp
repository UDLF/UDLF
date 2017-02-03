/* <Effectiveness.hpp>
 *
 * Effectiveness class header file
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

#ifndef EFFECTIVENESS_HPP
#define EFFECTIVENESS_HPP

#include <iostream>
#include <memory>
#include <map>
#include <vector>
#include <sstream>
#include <fstream>

class Effectiveness {
    public:
            Effectiveness(int& n_in,
                          std::vector<int>& rkLists_in,
                          std::vector<std::string>& imgList_in);
            //Recall
            float computeRecall(int recallAt);
            //MAP
            float computeMAPMeasure();
            //Precision
            void fillPrecisionsMap(std::map<int, float>& precisions, std::string precisionsToCompute);
            //Files
            void readClassesFile(std::string classFile);
            //Given an image, return its class name
            std::string getClass(int x);

    private:
            //Recall
            void computeShowResults(int recallAt);
            void computeAvgArray();
            float computeFinalRecall();
            //Precision
            float computePrecision(int k);
            //MAP
            float computeMAP(int paramNK);
            float computeAveragePrecision(int qId, int d, int offset);
            void startMAPByClass();
            void endMAPByClass();
            //Generic
            int getClassSize(std::string classname);

            //Variables
            const int& n;
            const std::vector<int>& rkLists;
            const std::vector<std::string>& imgList;
            std::map<std::string, std::string> classes;
            std::map<std::string, unsigned int> classesSize;
            std::map<std::string, float> mapByClass;

};

#endif // EFFECTIVENESS_HPP
