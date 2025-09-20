/*  Lattice Boltzmann sample, written in C++, using the OpenLB
 *  library
 *
 *  Copyright (C) 2025 Mingliang Zhong, Stephan Simonis
 *  E-mail contact: info@openlb.net
 *  The most recent release of OpenLB can be downloaded at
 *  <http://www.openlb.net/>
 *
 *  This program is free software; you can redistribute it and/or
 *  modify it under the terms of the GNU General Public License
 *  as published by the Free Software Foundation; either version 2
 *  of the License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public
 *  License along with this program; if not, write to the Free
 *  Software Foundation, Inc., 51 Franklin Street, Fifth Floor,
 *  Boston, MA  02110-1301, USA.
 */

#ifndef FILES_IO_H
#define FILES_IO_H

#include "filesIO.h"
#include <cstdlib> // For std::system

bool directoryExists(const std::string& path)
{
  struct stat statbuf;
  return (stat(path.c_str(), &statbuf) == 0 && S_ISDIR(statbuf.st_mode));
}

bool createDirectory(const std::string& path)
{
  return (mkdir(path.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) == 0);
}

bool deleteDirectory(const std::string& path)
{
  std::string command = "rm -rf " + path;
  return (std::system(command.c_str()) == 0);
}

bool fileExists(const std::string& name)
{
  std::ifstream f(name.c_str());
  return f.good();
}

#endif // FILES_IO_hh