/* <Exec.hpp>
 *
 * Exec class header file
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

#ifndef EXEC_HPP
#define EXEC_HPP

#include <memory>
#include <vector>
#include <map>

class Exec {
    public:
            Exec();
            static Exec& getInstance();
            bool parseFile(std::string filename);
            void insertElementInMap(std::string name, std::string val);
            void getConfigVariable(std::string& dest, std::string variable);
            void getConfigVariable(long int& dest, std::string variable);
            void getConfigVariable(int& dest, std::string variable);
            void getConfigVariable(double& dest, std::string variable);
            void getConfigVariable(bool& dest, std::string variable);
            void getInputFusionFiles(std::vector<std::string>& fusionFiles);
            std::string getMethodParameterVariables();
            void run();

    private:
            std::map<std::string, std::string> m_variables; //map to store config variables and their respective values

};

#endif // EXEC_HPP
