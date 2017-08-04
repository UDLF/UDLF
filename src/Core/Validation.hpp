/* <Validation.hpp>
 *
 * Validation class header file
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

#ifndef VALIDATION_HPP
#define VALIDATION_HPP

#include <memory>
#include <map>
#include <vector>
#include <array>

#include "Conf.cpp" //specifies the internal configuration files

class Validation {
    public:
            Validation();
            static Validation& getInstance();
            bool validate(std::map<std::string, std::string>& variables);

    private:
            struct varInfo {
                int type;                                                //0 - STR; 1 - INT; 2 - DBL; 3 - BOL
                bool userMustDefine;                                     //If true and the user hasn't defined, abort
                std::string defaultVal;                                  //Default value
                bool hasExpectedVal;                                     //False = can be anything of the type; True = need to be one of the specfied
                std::vector<std::string> expectedVals;                   //Expected values
                std::map<std::string, struct varInfo> vectorVariables;
            };

            bool fillValidationMap(std::string confStr);
            bool applyValidation(std::map<std::string, std::string>& variables,
                                 std::map<std::string, struct varInfo>& validationMap,
                                 std::vector<std::string>& variableSequence);
            bool parse(std::string str);
            int getTypeOf(std::string str);
            std::string getTypeOf(int i);
            void printVarInfo(varInfo var);
            bool insertElementInMap(std::string name, struct varInfo var);
            bool checkType(std::string str, int type);
            bool checkExpected(std::string value, std::vector<std::string> expected);

            std::map<std::string, struct varInfo> m_validationMap;
            std::vector<std::string> m_variableSequence; //keep the sequence when inserting values in the validationMap

            /* Internal configuration files */
            //these files are the base of validation and are embedded into the binary at compile time
            //main config file
            const std::string mainConfStr = confGeneral;
            //file for each method, string ID, name of the file
            const std::map<std::string, std::string> methodsConfStr = { {"CPRR", confCprr},
                                                                        {"RLRECOM", confRlrecom},
                                                                        {"RLSIM", confRlsim},
                                                                        {"CONTEXTRR", confContextrr},
                                                                        {"RECKNNGRAPH", confRecknngraph},
                                                                        {"RKGRAPH", confRkgraph},
                                                                        {"CORGRAPH", confCorgraph},
                                                                        {"NONE", confNone}
                                                                      };
            //in case you want to change the parameter that defines the method to be used, make sure to change it here too
            const std::string methodParam = "UDL_METHOD";
};

#endif // VALIDATION_HPP
