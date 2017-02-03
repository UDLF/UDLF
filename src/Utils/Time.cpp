/* <Time.cpp>
 *
 * Time class implementantion file
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

#include "Time.hpp"

/* Constructor */
Time::Time() {

}

float Time::addTime(timeval startTime, float currentTime) {
    timeval endTime;
    long seconds, nseconds;
    float duration, newTime;

    gettimeofday(&endTime, NULL);

    seconds  = endTime.tv_sec  - startTime.tv_sec;
    nseconds = endTime.tv_usec - startTime.tv_usec;

    duration = seconds + nseconds/1000000.0;
    newTime = currentTime + duration;

    return newTime;
}

std::string Time::getCurrentTime() {
    time_t t = time(0); // get the current time
    struct tm* now = localtime(&t); // convert to a tm struct

    std::stringstream time;
    time << (now->tm_year + 1900) << '/'  << (now->tm_mon + 1) << '/' <<  now->tm_mday
         << " " << now->tm_hour << ":" << now->tm_min << ":" << now->tm_sec << "\n";

    return time.str();
}
