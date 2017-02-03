/* <None.cpp>
 *
 * Just an empty method implementation
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

#include "None.hpp"

/* Constructor */
None::None() {

}

void None::loadParameters() {
    Exec exec = Exec::getInstance();
    exec.getConfigVariable(l, "PARAM_NONE_L");
}

void None::checkParameters() {
    if (l > n) {
        std::cout << "L can't be greater than N" << std::endl;
        std::cout << "Aborting..." << std::endl;
        exit(1);
    }
}

void None::initDataStructuresUdl() {
    std::cout << "Initializing data structures..." << std::endl;

    rkLists.resize(n*l);
    initMatrix(matrix);

    std::cout << "Initialized successfully!" << std::endl;
}

void None::initDataStructuresFusion() {
    std::cerr << "WARNING: The Fusion Mode is not implemented for this method! Aborting ...\n";
    exit(1);
}

void None::prepareInput() {
    if (inputFileFormat == "MATRIX") {
        if (inputMatrixType == "DIST") {
            genRksFromDistMatrix();
        } else { //SIM
            genRksFromSimMatrix();
        }
    }
}

void None::prepareOutput() {
    if (outputFileFormat == "MATRIX") {
        if (outputMatrixType == "DIST") {
            genDistMatrixFromRks();
        } else { //SIM
            genSimMatrixFromRks();
        }
    }
}

void None::runUdlMethod() {
    std::cout << "\n Executing NONE!\n\n";
}

void None::runFusionMethod() {
    std::cerr << "WARNING: The Fusion Mode is not implemented for this method! Aborting ...\n";
    exit(1);
}
