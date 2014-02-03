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

#ifndef CV_GRO_TRAJREADER_H
#define CV_GRO_TRAJREADER_H

#include <boost/array.hpp>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

#include "Gromacs.h"
#include "traj_reader.h"

// Iterates over GROMACS trajectory files (.xtc, .trr, .g96, .gro).
// Implements TrajReader.
// Uses Gromacs.h from VMD molfile plugin.
class Gromacs_TrajReader : public TrajReader
{
public:
    Gromacs_TrajReader() : mf(NULL), mdh_natoms(0) {};

    // Opens trajectory file for iteration
    bool open_traj_iter(const char *xtc_fname, int first, int last, int step, std::ostream &out_inf)
    {
        // Initialize iteration parameters
        ts_first = first;  // first timestep to process
        ts_last = last;    // last timestep to process
        ts_step = step;    // iteration step
        ts_cnt = 0;        // current timestep

        // Open Gromacs file with coordinates
        mf = mdio_open(xtc_fname, 0, MDIO_READ);
        if (!mf) {
            out_inf << "unable to open trajectory file : " << xtc_fname << " " << mdio_errmsg(mdio_errno()) << std::endl;
            return false;
        }

        // Read Gromacs header
        md_header mdh;
        if (mdio_header(mf, &mdh) < 0) {
            out_inf <<  "Cannot read header from " << xtc_fname;
            out_inf << " " << mdio_errmsg(mdio_errno()) << std::endl;
            return false;
        }

        // For .g96 files we don't have number of atoms in header so get it by counting
        if (mf->fmt == MDFMT_G96)
            mdh.natoms = g96_countatoms(mf);

        mdh_natoms = mdh.natoms;  // save number of atoms for future timestep reads
        return true;
    }


    // Iterates to the next requested timestep in trajectory and returns geometry and info througt arguments.
    // Returns: true if coordiates have actually been read, false on EOF and errors.
    bool next_timestep(Array_double_3 &box, Atoms_container &atoms, std::string &ts_info, std::ostream &out_inf)
    {
        md_ts mdts;  // molfile's timestep structure
        mdts.natoms = mdh_natoms;

        // read timesteps from file
        while (ts_cnt < ts_last && !(mdio_timestep(mf, &mdts) < 0)) {
            ++ts_cnt;

            // flag if current position has form first+i*step 
            bool out_conf = (ts_cnt >= ts_first && ((ts_cnt - ts_first) % ts_step == 0));  

            bool good_conf = true;  // flag if timestep contains geometry and was read without errors
            if (!mdts.box) {
                out_inf << "Timestep does not contain box dimensions!" << std::endl;
                good_conf = false;
            }

            // Output only requested configurations
            if (good_conf && out_conf) {
                // Format and store some information about current timestep
                std::ostringstream os;
                os << "step " << mdts.step << "  time " << std::setprecision(6) << mdts.time;
                ts_info = os.str();

                // Store the geometry
                box[0] = mdts.box->A;
                box[1] = mdts.box->B;
                box[2] = mdts.box->C;

                atoms.clear();  // clear atoms container before storing into it
                for (int i = 0; i < mdts.natoms; ++i) {
                    const float *a = mdts.pos + 3*i;
                    Array_double_3 x = {{a[0], a[1], a[2]}};
                    atoms.push_back(x);
                }
            }

            mdio_tsfree(&mdts);

            if (!good_conf)
                return false;
            if (out_conf)
                return true;
        }

        return false;
    }

    // Returns number of timesteps we have read so far (including skipped)
    int get_current_pos() const { return ts_cnt; }; 

    ~Gromacs_TrajReader() { if (mf != NULL) mdio_close(mf); }

private:
    md_file *mf;
    int mdh_natoms;  // number of atoms
    int ts_first, ts_last, ts_step;  // trajectory iteration parameters
    int ts_cnt;  // current timestep
};


#endif /* CV_GRO_TRAJREADER_H */
