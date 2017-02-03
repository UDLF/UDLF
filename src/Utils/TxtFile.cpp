/* <TxtFile.cpp>
 *
 * TxtFile class implementantion file
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

#include "TxtFile.hpp"

/* Constructor */
TxtFile::TxtFile() {

}

/* Given a text file, prints it in the screen (like the cat command) */
void TxtFile::printFile(std::string path) {
    std::ifstream inFile;
    inFile.open(path.c_str());
    if (!inFile) {
        std::cout << "Unable to open file [" << path.c_str() << "].";
        exit(1); // terminate with error
    }
    std::string line;
    while (std::getline(inFile, line)) {
        std::cout << line.c_str() << std::endl;
    }
    inFile.close();
}

/* Given a string, replaces a substring for a given string */
bool TxtFile::replace(std::string& str, const std::string& from, const std::string& to) {
    size_t start_pos = str.find(from);
    if (start_pos == std::string::npos) {
        return false;
    }
    str.replace(start_pos, from.length(), to);
    return true;
}
