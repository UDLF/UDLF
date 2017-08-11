/* <Validation.cpp>
 *
 * Validation class implementation file
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
#include <algorithm>
#include <sstream>

#include "Validation.hpp"
#include "Utils/Type.hpp"
#include "Utils/TxtFile.hpp"

/* Constructor */
Validation::Validation() {

}

/* Allows to call the class public methods in any class (singleton class) */
Validation& Validation::getInstance() {
    static Validation self;
    return self;
}

/* Based on the string, return corresponding int */
int Validation::getTypeOf(std::string str) {
    if (str == "STR") {
        return 0;
    }
    if (str == "INT") {
        return 1;
    }
    if (str == "UINT") {
        return 2;
    }
    if (str == "DBL") {
        return 3;
    }
    if (str == "BOL") {
        return 4;
    }
    return -1;
}

/* Based on the int, return corresponding string */
std::string Validation::getTypeOf(int i) {
    switch (i) {
        case 0:
            return "STR";
        case 1:
            return "INT";
        case 2:
            return "UINT";
        case 3:
            return "DBL";
        case 4:
            return "BOL";
        default:
            return "";
    }
}

/* Verifies if the variable has the correct type */
bool Validation::checkType(std::string str, int type) {
    switch (type) {
        case 0: //STR
            return true;
        case 1: //INT
            return Type::isInteger(str);
        case 2: //UINT
            return Type::isIntegerPos(str);
        case 3: //DBL
            return Type::isNumeric(str);
        case 4: //BOL
            return Type::isBoolean(str);
        default: //Unknown type
            return false;
    }
}

/* Inserts a new element in the map and verify if it is duplicated */
bool Validation::insertElementInMap(std::string name, struct varInfo var) {
    auto it = m_validationMap.find(name);
    if (it != m_validationMap.end()) {
        std::cerr << name << ": there is already a variable with this name in the validation map!" << std::endl;
        return false;
    }

    m_validationMap[name] = var;
    m_variableSequence.push_back(name);
    return true;
}

/* Checks if teh variable has a value that matches one of the expected values */
bool Validation::checkExpected(std::string value, std::vector<std::string> expected) {
    for (std::string str : expected) {
        if (value == str) {
            return true;
        }
    }
    return false;
}

/* Parses a line of config file, adding a new element to the validation map */
bool Validation::parse(std::string line) {
    std::string name;
    struct varInfo var;

    var.vectorVariables.clear();
    std::string varname;
    int posOpenBracket = line.find("{", 0);
    if (posOpenBracket != -1) {
        line = line.substr(posOpenBracket+1, line.length());
        int posCloseBracket = line.find("}", 0);
        varname = line.substr(0, posCloseBracket);
        line = line.substr(posCloseBracket+1, line.length());
    }

    int posDots = line.find(":", 0);
    int posEq   = line.find("=", 0);
    int posPar  = line.find("(", 0);

    if ( posDots == -1 || (posPar != -1 && posEq != -1) ) {
        return false;
    }

    name = line.substr(0, posDots);

    if (posPar != -1) {
        int posPar2  = line.find(")", 0);
        if (posPar2 == -1) {
            return false;
        }

        var.type = getTypeOf(line.substr(posDots+1, posPar-(posDots+1)));
        var.userMustDefine = false;
        var.hasExpectedVal = true;

        std::istringstream f(line.substr(posPar+1, posPar2-(posPar+1)));
        std::string s;
        while (std::getline(f, s, ';')) {
            var.expectedVals.push_back(s);
        }
        var.defaultVal = var.expectedVals[0];
    } else if (posEq != -1) {
        var.type = getTypeOf(line.substr(posDots+1, posEq-(posDots+1)));
        var.userMustDefine = false;
        var.hasExpectedVal = false;
        var.defaultVal     = line.substr(posEq+1, line.length());
    } else {
        var.type = getTypeOf(line.substr(posDots+1, line.length()));
        var.userMustDefine = true;
    }

    if (var.type == -1) {
        return false;
    }

    if (posOpenBracket == -1) {
        if (!insertElementInMap(name, var)) {
            return false;
        }
    } else {
        auto it = m_validationMap.find(varname);
        if (it == m_validationMap.end()) {
            std::cerr << "ERROR: Can't define " << name << "! The variable " << varname << " hasn't been declared yet!\n";
            return false;
        }
        std::string type = getTypeOf(it->second.type);
        if ( (type != "INT") && (type != "UINT") ) {
            std::cerr << "ERROR: Can't define " << name << "! The variable " << varname << " is not an integer!\n";
            return false;
        }
        it->second.vectorVariables[name] = var;
        return true;
    }

    return true;
}

/* Given a config file, parses the validation array and fill the validation map */
bool Validation::fillValidationMap(std::string confStr) {
    m_validationMap.clear();
    m_variableSequence.clear();

    std::istringstream file(confStr);

    int lineNum = 1;
    std::string line;
    while (std::getline(file, line)) {
        line.erase(std::remove(line.begin(),line.end(),' '),line.end());  //remove spaces
        line.erase(std::remove(line.begin(),line.end(),'\t'),line.end()); //remove tabs
        line = line.substr(0, line.find("#", 0)); //remove comments
        if (line != "") {
            if (!parse(line)) {
                std::cerr << "Internal config file is corrupted! Problem in line " << lineNum << std::endl;
                return false;
            }
        }
        lineNum++;
    }

    return true;
}

/* From the filled validation map, validates the values of the parameters set by the user */
bool Validation::applyValidation(std::map<std::string, std::string>& variables, std::map<std::string, struct varInfo>& validationMap, std::vector<std::string>& variableSequence) {
    for (std::string name : variableSequence) {
        struct varInfo var = validationMap[name];

        auto it = variables.find(name);
        std::string& value = variables[name];
        if (var.userMustDefine) {
            if (it == variables.end()) {
                std::cerr << " *ERROR: You haven't declared " << name << "! You must declare it! Aborting..." << std::endl;
                return false;
            } else {
                std::cout << " " << name << " = " << value << std::endl;
            }
        } else {
            if (it == variables.end()) {
                variables[name] = var.defaultVal;
                std::cerr << " *WARNING: You haven't declared " << name << "!\n *SETTING: " << name << " = " << var.defaultVal << std::endl;
            } else if (var.hasExpectedVal) {
                if (!checkExpected(value, var.expectedVals)) {
                    std::cerr << " *WARNING: You've declared " << name << " = " << value << "!";
                    std::cerr << " Expected values are: [";
                    for (std::string str : var.expectedVals) {
                        std::cout << str << ",";
                    }
                    std::cerr << "\b]\n *SETTING: " << name << " = " << var.defaultVal << std::endl;
                    variables[name] = var.defaultVal;
                } else {
                    std::cout << " " << name << " = " << value << std::endl;
                }
            } else {
                std::cout << " " << name << " = " << value << std::endl;
            }
        }

        if (!checkType(value, var.type)) {
            std::cerr << " *WARNING: Type of " << name << " is " << getTypeOf(var.type) << "! " << value << " is invalid!" << std::endl;
            if (var.userMustDefine) {
                std::cerr << "This parameter needs to be specified correctly! Aborting..." << std::endl;
                return false;
            } else {
                variables[name] = var.defaultVal;
                std::cerr << " *SETTING: " << name << " = " << var.defaultVal << std::endl;
                if (!checkType(var.defaultVal, var.type)) {
                    std::cerr << " WARNING: Type of default value doesn't match with the var type! Aborting..." << std::endl;
                    return false;
                }
            }
        }

        if (!var.vectorVariables.empty()) {
            int value = std::stoi(variables[name]);
            std::map<std::string, struct varInfo> tmp;
            std::vector<std::string> tmpSeq;
            tmp.clear();
            for (auto it : var.vectorVariables) {
                std::string valueStr = it.second.defaultVal;
                for (int i = 1; i <= value; i++) {
                    std::stringstream ss;
                    ss << i;
                    std::string name = it.first;
                    TxtFile::replace(name, "*", ss.str());
                    if (!it.second.userMustDefine) {
                        std::string curValue = valueStr;
                        TxtFile::replace(curValue, "*", ss.str());
                        it.second.defaultVal = curValue;
                    }
                    tmp[name] = it.second;
                    tmpSeq.push_back(name);
                }
            }
            applyValidation(variables, tmp, tmpSeq);
        }

    }

    return true;
}

/* Validates the variables map of the Exec class and try to fix issues. If can't correct the issues, returns false */
bool Validation::validate(std::map<std::string, std::string>& variables) {
    std::cout << "Validating parameters...\n" << std::endl;

    /* MAIN VALIDATION */
    std::cout << "General parameters: " << std::endl;

    if (!fillValidationMap(mainConfStr)) {
        std::cerr << "Couldn't fill validation map! Check your config/*.conf files and try to recompile your executable!" << std::endl;
        return false;
    }

    if (!applyValidation(variables, m_validationMap, m_variableSequence)) {
        return false;
    }


    /* METHOD VALIDATION */
    std::cout << "\nMethod specific parameters: " << std::endl;

    auto it = variables.find(methodParam); //the method must be defined there at this point
    std::string method = it->second;

    auto it2 = methodsConfStr.find(method);
    std::string methodConfStr = it2->second;
    if (!fillValidationMap(methodConfStr)) {
        std::cerr << "Couldn't fill validation map! Check your config/*.conf files and try to recompile your executable!" << std::endl;
        return false;
    }

    if (!applyValidation(variables, m_validationMap, m_variableSequence)) {
        return false;
    }

    return true;
}
