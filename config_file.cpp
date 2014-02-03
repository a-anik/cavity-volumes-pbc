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

#include <fstream>
#include <iostream>
#include <iomanip>
#include <iterator>
#include <utility>
#include <vector>
#include <limits>
#include <boost/iostreams/device/null.hpp>
#include <boost/iostreams/tee.hpp>

#include "cfglinestream.h"
#include "config_file.h"
#include "traj_reader_gromacs.h"
#include "traj_reader_nice.h"

using boost::iostreams::tee_filter;


// Reads run control file, radii, opens input and output files
bool CavConfig::init(const char *run_control_file)
{
    std::string trj_fname, r_fname;
    int ts_first, ts_last, ts_step;
    int do_verbose, file_format;

    // Read run control file
    Cfglinestream cfg(run_control_file);

    cfg >> trj_fname;
    cfg >> ts_first >> ts_last >> ts_step;
    cfg >> r_fname;
    cfg >> file_format;
    cfg >> shell_width;
    cfg >> r_scale;
    cfg >> r_probe;

    cfg >> out_inf_ofstream;
    cfg >> do_verbose;
    cfg >> out_asf;
    cfg >> out_vol;
    cfg >> out_evf;
    cfg >> out_stl;
    cfg >> nSubdiv;

    // Initialize log streams (tee to file and screen)
    if (out_inf_ofstream.good())
        out_inf.push(tee_filter<std::ofstream>(out_inf_ofstream));

    if (do_verbose)
        out_inf.push(std::cout);
    else
        out_inf.push(boost::iostreams::null_sink());

    // Sanitize input numbers for cycling over configurations
    if (ts_first <= 0) ts_first = 1;
    if (ts_step <= 0) ts_step = 1;
    if (ts_last < ts_first)
        ts_last = std::numeric_limits<int>::max();

    out_inf << "Trajectory file : " << trj_fname << std::endl;
    out_inf << "traj first/last/step : " << ts_first << " " << ts_last << " " << ts_step << std::endl;
    out_inf << "file format : " << file_format << std::endl;
    out_inf << "Radii file : " << r_fname << std::endl;

    // Read radii file. Format:
    //   #header line
    //   N
    //   r_1
    //   ...
    //   r_N
    std::ifstream in_rad(r_fname.c_str(), std::ios::in);
    if (!in_rad) {
        out_inf << "unable to open file with radii : " << r_fname << std::endl;
        return false;
    }
    in_rad.ignore(std::numeric_limits<int>::max(), '\n');
    int nA_rad;
    in_rad >> nA_rad;
    for (int i = 0; i < nA_rad; i++) {
        double r;
        if (!(in_rad >> r)) {
            out_inf << "radii read error!" << std::endl;
            return false;
        }
        radii.push_back(r);
    }

    out_inf << "r[0] ... r[" << radii.size() << "] : " << radii.front() << " ... " << radii.back() << std::endl;
    out_inf << "Shell width : " << shell_width << std::endl;
    out_inf << "R scale factor : " << r_scale << std::endl;
    out_inf << "Rprobe : " << r_probe << std::endl;
    out_inf << "Subdivision depth : " << nSubdiv << std::endl;

    // Select and initialize trajectory reader
    switch (file_format) {
        case 0: trj = new Gromacs_TrajReader(); break;
        default: trj = new NICE_TrajReader(); break;
    }

    if (!trj->open_traj_iter(trj_fname.c_str(), ts_first, ts_last, ts_step, out_inf))
        return false;

    // Check output files and set precision
    if (!out_asf || !out_vol || !out_evf || !out_stl)
        return false;
    out_asf.precision(std::numeric_limits<double>::digits10);
    out_vol.precision(std::numeric_limits<double>::digits10);
    out_evf.precision(std::numeric_limits<double>::digits10);
    out_inf.precision(std::numeric_limits<double>::digits10);

    // Resize storage for per-atom contributions to void surface
    atom_confs_surf.resize(radii.size(), 0.0L);

    return true;
}

// Retruns total number of timestep reads from trajectory
int CavConfig::traj_ts_cnt()
{
    int ts_cnt = 0;
    if (trj != NULL)
        ts_cnt = trj->get_current_pos();

    return ts_cnt;
}

// Reads box, atom coordinates, info from the next timestep iteration.
bool CavConfig::next_timestep()
{
    bool good_conf = trj->next_timestep(box, atoms, ts_info, out_inf);

    // Consistency checks: number of atoms and box presence
    if (good_conf && atoms.size() != radii.size()) {
        out_inf << "Number of atoms in timestep is not equal to number of radii : " << atoms.size();
        out_inf << " vs. " << radii.size() << std::endl;
        good_conf = false;
    }

    return good_conf;
}


CavConfig::~CavConfig()
{
    if (trj != NULL)
        delete trj;
}
