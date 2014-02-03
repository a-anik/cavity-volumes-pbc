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

#ifndef CV_NICE_TRAJREADER_H
#define CV_NICE_TRAJREADER_H

#include <boost/array.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "traj_reader.h"

// Iterates over text trajectory files. 
// File format:
//   #timestep1 comment
//   N_atoms
//   box_X
//   box_Y
//   box_Z
//   x1 y1 z1
//   ......
//   xN yN zN
//   #timestep2 comment
//   ......
// Implements TrajReader.
class NICE_TrajReader : public TrajReader
{
public:
    // Opens trajectory file for iteration
    bool open_traj_iter(const char *fname, int first, int last, int step, std::ostream &out_log)
    {
        // Initialize iteration parameters
        ts_first = first;  // first timestep to process
        ts_last = last;    // last timestep to process
        ts_step = step;    // iteration step
        ts_cnt = 0;        // current timestep

        in_dat.open(fname);
        if (!in_dat) {
            out_log << "unable to open trajectory file : " << fname << std::endl;
            return false;
        }
        return true;
    }

    // Iterates to the next requested timestep in trajectory and returns geometry and info througt arguments.
    // Returns: true if coordiates have actually been read, false on EOF and errors.
    bool next_timestep(Array_double_3 &box, Atoms_container &atoms, std::string &comment, std::ostream &out_log)
    {
        // read timesteps from file
        while (ts_cnt < ts_last) {
            std::string line;
            size_t nA = 0;

            getline(in_dat, comment);                               // comment line
            getline(in_dat, line);                                  
            static_cast<std::istringstream>(line) >> nA;
            getline(in_dat, line);
            static_cast<std::istringstream>(line) >> box[0];
            getline(in_dat, line);
            static_cast<std::istringstream>(line) >> box[1];
            getline(in_dat, line);
            static_cast<std::istringstream>(line) >> box[2];

            atoms.clear();  // clear atoms container before storing into it
            for (size_t i = 0; in_dat.good() && i < nA; i++) {
                Array_double_3 x;
                getline(in_dat, line);
                std::istringstream is(line);
                is >> x[0] >> x[1] >> x[2];
                atoms.push_back(x);
            }

            if (!in_dat) {
                out_log << "Trajectory read failed." << std::endl;
                return false;
            }

            ++ts_cnt;

            if (atoms.size() != nA) {
                out_log << "Number of coordinates in timestep is not nA :" << atoms.size() << " vs. " << nA << std::endl;
                return false;
            }

            bool out_conf = (ts_cnt >= ts_first && ((ts_cnt - ts_first) % ts_step == 0));
            if (out_conf)
                return true;
        }

        return false;
    }

    // Returns number of timesteps we have read so far (including skipped)
    int get_current_pos() const { return ts_cnt; }; 

private:
    std::fstream in_dat;
    int ts_first, ts_last, ts_step;  // trajectory iteration parameters
    int ts_cnt;  // current timestep
};

#endif /* CV_NICE_TRAJREADER_H */
