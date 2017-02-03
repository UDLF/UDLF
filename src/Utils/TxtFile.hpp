/* <TxtFile.hpp>
 *
 * TxtFile class header file
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

#ifndef TXTFILE_HPP
#define TXTFILE_HPP

#include <iostream>
#include <sstream>
#include <fstream>

class TxtFile {
    public:
            TxtFile();
            static void printFile(std::string path);
            static bool replace(std::string& str, const std::string& from, const std::string& to);

    private:

};

#endif // TXTFILE_HPP
