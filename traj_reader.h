// Copyright (C) 2014 Alexey Anikeenko
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#ifndef CV_TRAJREADER_H
#define CV_TRAJREADER_H

#include <boost/array.hpp>
#include <vector>
#include <iostream>
#include <string>

// Interface class for iterating over MD trajectory
class TrajReader
{
public:
    typedef boost::array<double, 3> Array_double_3;
    typedef std::vector<Array_double_3> Atoms_container;

public:
    virtual bool open_traj_iter(const char *fname, int first, int last, int step, std::ostream &out_log) = 0;
    virtual bool next_timestep(Array_double_3 &box, Atoms_container &atoms, std::string &comment, std::ostream &out_log) = 0;
    virtual int get_current_pos() const = 0;  // current position in trajectory (number of timesteps from the begining)
    virtual ~TrajReader() {}
};

#endif /* CV_TRAJREADER_H */
