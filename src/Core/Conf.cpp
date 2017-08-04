/* <Conf.cpp>
 *
 * Load config/*.conf (internal configuration) files into char buffers 
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

const char confGeneral[] =
    #include "../config/general.conf"
;

const char confContextrr[] =
    #include "../config/contextrr.conf"
;

const char confCorgraph[] =
    #include "../config/corgraph.conf"
;

const char confCprr[] =
    #include "../config/cprr.conf"
;

const char confNone[] =
    #include "../config/none.conf"
;

const char confRecknngraph[] =
    #include "../config/recknngraph.conf"
;

const char confRkgraph[] =
    #include "../config/rkgraph.conf"
;

const char confRlrecom[] =
    #include "../config/rlrecom.conf"
;

const char confRlsim[] =
    #include "../config/rlsim.conf"
;
