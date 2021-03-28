/* <Exec.cpp>
 *
 * Exec class implementation file
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

#include <fstream>
#include <iostream>
#include <string>
#include <sstream>
#include <algorithm>

#include "Exec.hpp"
#include "Validation.hpp"
#include "Utils/Type.hpp"
#include "Methods/Udl.hpp"
#include "Methods/None.hpp"
#include "Methods/RDPAC.hpp"
#include "Methods/BFSTree.hpp"
#include "Methods/LHRR.hpp"
#include "Methods/Cprr.hpp"
#include "Methods/RlRecom.hpp"
#include "Methods/RlSim.hpp"
#include "Methods/Contextrr.hpp"
#include "Methods/ReckNNGraph.hpp"
#include "Methods/RkGraph.hpp"
#include "Methods/CorrelationGraph.hpp"

/* Constructor */
Exec::Exec() {

}

/* Allows to call the class public methods in any class (singleton class) */
Exec& Exec::getInstance() {
    static Exec self;
    return self;
}

/* Inserts a new element in the map and verify if it is duplicated */
void Exec::insertElementInMap(std::string name, std::string val) {
    auto it = m_variables.find(name);
    if (it != m_variables.end()) {
        std::cerr << name << " is being declared again!" << std::endl;
    }

    m_variables[name] = val;
}

/* Parses config file and store variable values in a map */
bool Exec::parseFile(std::string filename) {
    std::cout << "\nInput configuration file: " << filename << std::endl;

    std::ifstream file;
    file.open(filename);
    if (!file.is_open()) {
        std::cerr << "Configuration file not found!\n";
        return false;
    }

    int lineNum = 1;
    std::string line;
    while (std::getline(file, line)) { //read line
        line.erase(std::remove(line.begin(),line.end(),' '),line.end());  //remove spaces
        line.erase(std::remove(line.begin(),line.end(),'\t'),line.end()); //remove tabs
        line = line.substr(0, line.find("#", 0)); //remove comments

        if (line != "") {
            int pos = line.find("=");
            if (pos != -1) {
                std::string variable = line.substr(0, pos);
                std::string value = line.substr(pos+1, line.size());
                insertElementInMap(variable, value);  //save value to variable in the map
                //std::cout << variable << " set to " << value << std::endl;
            } else {
                std::cout << "Couldn't parse line " << lineNum << "! No attribution found!" << std::endl;
            }
        }

        lineNum++;
    }

    std::cout << std::endl;

    Validation& validation = Validation::getInstance();
    if (!validation.validate(m_variables)) {
        std::cerr << "Validation failed!" << std::endl;
        return false;
    }

    std::cout << std::endl;

    return true;
}

/* Gives the value of a variable of the config file as string */
void Exec::getConfigVariable(std::string& dest, std::string variable) {
    dest = "";

    auto it = m_variables.find(variable);
    if (it != m_variables.end()) {
        dest = it->second;
    } else {
        std::cerr << "Variable " << variable << " couldn't be found!" << std::endl;
    }
}

/* Gives the value of a variable of the config file as long int */
void Exec::getConfigVariable(long int& dest, std::string variable) {
    dest = 0;

    auto it = m_variables.find(variable);
    if (it != m_variables.end()) {
        if (!Type::isInteger(it->second)) {
            std::cerr << variable << " (" << it->second << ") can't be converted to long int! Returning 0!" << std::endl;
            return;
        }
        dest = stoi(it->second);
    } else {
        std::cerr << "Variable " << variable << " couldn't be found!" << std::endl;
    }
}

/* Gives the value of a variable of the config file as int */
void Exec::getConfigVariable(int& dest, std::string variable) {
    dest = 0;

    auto it = m_variables.find(variable);
    if (it != m_variables.end()) {
        if (!Type::isInteger(it->second)) {
            std::cerr << variable << " (" << it->second << ") can't be converted to int! Returning 0!" << std::endl;
            return;
        }
        dest = stoi(it->second);
    } else {
        std::cerr << "Variable " << variable << " couldn't be found!" << std::endl;
    }
}

/* Gives the value of a variable of the config file as double */
void Exec::getConfigVariable(double& dest, std::string variable) {
    dest = 0;

    auto it = m_variables.find(variable);
    if (it != m_variables.end()) {
        if (!Type::isNumeric(it->second)) {
            std::cerr << variable << " (" << it->second << ") can't be converted to double! Returning 0!" << std::endl;
            return;
        }
        dest = stod(it->second);
    } else {
        std::cerr << "Variable " << variable << " couldn't be found!" << std::endl;
    }
}

/* Gives the value of a variable of the config file as boolean */
void Exec::getConfigVariable(bool& dest, std::string variable) {
    dest = false;

    auto it = m_variables.find(variable);
    if (it != m_variables.end()) {
        if (!Type::isBoolean(it->second)) {
            std::cerr << variable << " (" << it->second << ") can't be converted to boolean! Returning false!" << std::endl;
            return;
        }
        dest = Type::convertToBoolean(it->second);
    } else {
        std::cerr << "Variable " << variable << " couldn't be found!" << std::endl;
    }
}

/* Fills a std::vector with the filenames of the input fusion files */
void Exec::getInputFusionFiles(std::vector<std::string>& fusionFiles) {
    fusionFiles.clear(); //just to be extra careful

    int numFusionFiles;
    getConfigVariable(numFusionFiles, "NUM_INPUT_FUSION_FILES");
    if (numFusionFiles == 0) {
        std::cerr << "WARNING: The Fusion mode requires at least one input file! Aborting ...\n";
        exit(1);
    }
    for (int i = 1; i <= numFusionFiles; i++) {
        std::string varName = "INPUT_FILES_FUSION_" + std::to_string(i);
        std::string varValue;
        getConfigVariable(varValue, varName);
        fusionFiles.push_back(varValue);
    }
}

/* Returns a string with the names and values of the parameter variables of the choosen method */
std::string Exec::getMethodParameterVariables() {
    std::stringstream ss;

    std::string delim;
    getConfigVariable(delim, "UDL_METHOD");
    delim = "PARAM_" + delim; //parameter variables MUST start with "PARAM_METHODNAME_"

    for (auto const& it : m_variables) {
        if (it.first.substr(0, delim.length()) == delim) {
            ss << "\n " << it.first << " = " << it.second;
        }
    }

    return ss.str();
}

/* Runs the application based in the info given in the file */
void Exec::run() {
    std::cout << "***********************************************************************\n\n";

    std::cout << "Running Method: ";

    std::string method;
    getConfigVariable(method, "UDL_METHOD");
    std::cout << method << "\n\n";

    if (method == "NONE") {
        None none;
        none.run();
    } else if (method == "RDPAC") {
        RDPAC rdpac;
        rdpac.run();
    } else if (method == "BFSTREE") {
        BFSTree bfstree;
        bfstree.run();
    } else if (method == "LHRR") {
        LHRR lhrr;
        lhrr.run();
    } else if (method == "CPRR") {
        Cprr cprr;
        cprr.run();
    } else if (method == "RLRECOM") {
        RlRecom rlRecom;
        rlRecom.run();
    } else if (method == "RLSIM") {
        RlSim rlSim;
        rlSim.run();
    } else if (method == "CONTEXTRR") {
        Contextrr contextrr;
        contextrr.run();
    } else if (method == "RECKNNGRAPH") {
        ReckNNGraph reckNNGraph;
        reckNNGraph.run();
    } else if (method == "RKGRAPH") {
        RkGraph rkGraph;
        rkGraph.run();
    } else if (method == "CORGRAPH") {
        CorrelationGraph corGraph;
        corGraph.run();
    }

    std::cout << "\nExecution Completed!" << std::endl;
}
