/* <Type.cpp>
 *
 * Type class implementantion file
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
#include <locale>

#include "Type.hpp"

/* Constructor */
Type::Type() {

}

/* Given a string returns a boolean (need to put this in a utils class) */
bool Type::convertToBoolean(std::string str) {
    std::locale loc;
    std::string s = "";
    for (std::string::size_type i=0; i<str.length(); ++i)
        s += std::tolower(str[i],loc);

    if (s == "true" || s == "1")
        return true;

    return false;
}

/* Tests if a string can be converted to an integer */
bool Type::isInteger(const std::string& str) {
    if (str.empty() || ((!isdigit(str[0])) && (str[0] != '-') && (str[0] != '+')))
        return false;

    char* p ;
    strtol(str.c_str(), &p, 10) ;

    return (*p == 0);
}

/* Tests if a string can be converted to an unsigned integer */
bool Type::isIntegerPos(const std::string& str) {
    if (str.empty() || ((!isdigit(str[0])) && (str[0] != '+')))
        return false;

    char* p ;
    strtol(str.c_str(), &p, 10) ;

    return (*p == 0);
}

/* Tests if a string can be converted to an double (numeric) */
bool Type::isNumeric(const std::string& str) {
    std::istringstream ss(str);
    double dbl;
    ss >> dbl;      // try to read the number
    ss >> std::ws;  // eat whitespace after number

    if (!ss.fail() && ss.eof()) {
        return true;  // is-a-number
    } else {
        return false; // not-a-number
    }
}

/* Tests if a string can be converted to a boolean */
bool Type::isBoolean(const std::string& str) {
    std::locale loc;
    std::string s = "";
    for (std::string::size_type i=0; i<str.length(); ++i)
        s += std::tolower(str[i],loc);

    if ( (s == "true" || s == "false") || (s == "1" || s == "0") )
        return true;

    return false;
}

/* Returns the number of digits of a positive integer number */
int Type::numDigits(int number) {
    int digits = 0;
    if (number < 0) {
        digits = 1;
    }
    while (number) {
        number /= 10;
        digits++;
    }
    return digits;
}
