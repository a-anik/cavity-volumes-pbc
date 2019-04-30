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

#ifndef CV_CONFIG_FIN_H
#define CV_CONFIG_FIN_H

#include <boost/array.hpp>
#include <boost/iostreams/filtering_stream.hpp>
#include <vector>
#include <string>
#include <iostream>
#include <fstream>

#include "traj_reader.h"

// Class for storing run control parameters and current atomic configuration.
class CavConfig
{
public:
    typedef boost::array<double, 3> Array_double_3;
    typedef std::vector<Array_double_3> Atoms_container;
    typedef std::vector<double> Radii_container;

public:
    double r_scale;          // input radii scale factor
    double r_probe;          // probe atom radius
    int do_pbc;              // do pbc unwrapping of coordinates

    std::vector<long double> atom_confs_surf;         // per atom contributions to void surface
    boost::iostreams::filtering_ostream out_inf;      // information logging tee stream
    std::ofstream out_asf;      // atom surface
    std::ofstream out_vol;      // raw cavity volumes
    std::ofstream out_stl;      // 3D surface file
    int nSubdiv;                // 3D surface subdivision parameter
    double stl_scale;           // output coordinates scale factor

    // current timestep data
    Array_double_3 box;
    Radii_container radii;
    Atoms_container atoms;
    std::string ts_info;


public:
    bool init(const char *run_control_file);
    bool next_timestep();       // fills next timestep data
    int traj_ts_cnt();          // number of timesteps we have read from trajectory (including skipped)

    CavConfig() : trj(NULL) {};
    ~CavConfig();

private:
    std::ofstream out_inf_ofstream;     // file for run information output
    TrajReader *trj;
};

#endif /* CV_CONFIG_FIN_H */
