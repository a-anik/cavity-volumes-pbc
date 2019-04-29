// Copyright (C) 2014,2019 Alexey Anikeenko
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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Timer.h>
#include <boost/array.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <cassert>
#include <vector>
#include <set>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <limits>
#include <iterator>
#include <algorithm>
#include <functional>
#include <utility>

#include "Regular_triangulation_cell_base_with_id_and_wcc_3.h"
#include "sastry_subsimplex_volume.h"
#include "skin_surface.h"
#include "config_file_fin.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef struct { int atom_id; }                             Point_info;
typedef K::FT                                               Weight;
typedef K::Point_3                                          Point;
typedef K::Weighted_point_3                                 Weighted_point;
typedef std::pair<Weighted_point, Point_info>               Weighted_point_with_info;
typedef std::vector<Weighted_point_with_info>               Wpi_container;

typedef CGAL::Regular_triangulation_vertex_base_3<K>                    Vb0;
typedef CGAL::Triangulation_vertex_base_with_info_3<Point_info, K, Vb0> Vb;
typedef CGAL::Regular_triangulation_cell_base_with_id_and_wcc_3<K>      Cb;  // no hidden points, int id, cached weighted circumcenter
typedef CGAL::Triangulation_data_structure_3<Vb, Cb>        Tds;
typedef CGAL::Regular_triangulation_3<K, Tds>               Rt;

typedef Rt::Finite_cells_iterator                           Finite_cells_iterator;
typedef Rt::Cell_handle                                     Cell_handle;
typedef CavConfig::Array_double_3                           Array_double_3;
typedef boost::array<double, 4>                             Array_double_4;  // for per atom surface in cell


struct Voids_result
{
    struct Void_info {
        long double volume;      // volume of each cavity
        long double surface;     // area of each cavity
        std::set<int> atoms;    // set of atoms for each cavity
        Void_info () : volume(0), surface(0), atoms() {}
    };
    std::vector <Void_info> voids;        // measurement results for discovered voids (num_voids size)
    std::vector<long double> atom_surf;   // for each atom : surface exposed to cavities (natoms values)
    std::ostream &out_log;                // information stream (debugging Voronoi calculation info)
    Voids_result(int natoms, std::ostream &strm) : atom_surf(natoms), out_log(strm) {}
};

bool void_greater(const Voids_result::Void_info& v1, const Voids_result::Void_info& v2) { return v1.volume > v2.volume; }

// Calculates void volume and spheres area (only for own four spheres) of regular triangulation cell.
// Doesn't consider extraneous atoms from neighbor cells! True volume can be obtained only for cluster of cells comprising one connected void.
// S. Sastry, D.S. Corti, P.G. Debenedetti, F.H. Stillinger, "Statistical geometry of particle packings. I. Algorithm for exact determination of connectivity,
// volume, and surface areas of void space in monodisperse and polydisperse sphere packings", Phys. Rev. E, V.56, N5, p. 5524, 1997. doi:10.1103/PhysRevE.56.5524
double cell_void_volume(Cell_handle c, double &out_surf, Array_double_4 &out_atom_surf)
{
    const Point *pts[4] = { &((c->vertex(0)->point()).point()),
                            &((c->vertex(1)->point()).point()),
                            &((c->vertex(2)->point()).point()),
                            &((c->vertex(3)->point()).point()) };
    const Point &V = c->weighted_circumcenter().point();   // Voronoi vertex V

    double volume = 0.0, surface = 0.0;
    Array_double_4 per_atom_surface;
    std::fill(per_atom_surface.begin(), per_atom_surface.end(), 0.0);

    // General idea: consider 24 subsimplexes with three orthogonal edges 
    // by constructing normals from weighted circumcenter to the faces and then to the edges
    for (int i = 0; i < 4; i++) {                             // iteration over cell facets
        const Point *p_vertex_i = pts[i];
        pts[i] = &V;     // replace vertex(i) in array  by weighted circumcenter of cell
        CGAL::Sign sV = orientation(*pts[0], *pts[1], *pts[2], *pts[3]);  // do V and vertex(i) on the same side of facet?

        const int a[3] = { (i+1)%4, (i+2)%4, (i+3)%4 };    // indicies of atoms of face opposite to cell vertex i
        const Point E = K::Plane_3(*pts[a[0]], *pts[a[1]], *pts[a[2]]).projection(V);  // projection of Voronoi vertex to facet
        const double z0 = sqrt((V-E).squared_length());    // length of normal to facet

        for (int j = 0; j < 3; j++) {                      // iteration over facet edges
            const int e[2] = { a[(j+1)%3], a[(j+2)%3] };   // indices of two edge atoms
            const Weighted_point WA[2] = { c->vertex(e[0])->point(), c->vertex(e[1])->point() };
            const Point A[2] = { WA[0].point(), WA[1].point() };
            const CGAL::Sign sE = coplanar_orientation(A[0], A[1], *pts[a[j]], E);

            const Point B = K::Line_3(A[0], A[1]).projection(E);  // projection of E to edge
            const double y0 = sqrt((E-B).squared_length());

            for (int k = 0; k < 2; k++) {                         // iteration over edge ends (two atoms)
                CGAL::Sign sB = CGAL::sign((B - A[k])*(A[(k+1)%2] - A[k]));  // negative if point B is in other direction than opposite A[(k+1)%2] looking from A[k]
                CGAL::Sign sF = sV * sE * sB;                                // final subsimplex sign
                double x0 = sqrt((B-A[k]).squared_length());
                double area;
                volume += sF * subsimplex_void_volume(x0, y0, z0, WA[k].weight(), &area);  // sum signed volume contribution from subsimplexes
                surface += sF * area;
                per_atom_surface[e[k]] += sF * area;
            }
        }

        pts[i] = p_vertex_i;  // place vertex(i) back to array
    }

    out_atom_surf = per_atom_surface;
    out_surf = surface;
    return volume;
}



// Constructs regular triangulation of weighted points.
// Builds graph of cells with voids. Finds connected components of this graph (voids as clusters of cells).
// For each component calculates total void volume and area.
// Excludes roughs on the surface (components, connected to the infinite cell).
bool regular_triangulation_voids(const Wpi_container& points, Voids_result &res)
{
    typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Cell_handle> CellGraph;
    typedef typename boost::graph_traits<CellGraph>::vertex_descriptor GVertex;
    typedef typename boost::graph_traits<CellGraph>::vertex_iterator GVertex_iterator;

    Rt T;            // regular triangulation
    CellGraph G;     // Voronoi subgraph where cells_sqr_r > 0 and edges_sqr_r > 0 (cells and edges with voids)
    GVertex_iterator vi, vi_end;


    res.out_log << "Number of input points for RT : " <<  points.size() << std::endl;
    T.insert(points.begin(), points.end());    // insert all points in a row (this is faster than one insert() at a time).
    T.infinite_vertex()->info().atom_id = -1;  // set special atom_id for dummy infinite triangulation vertex

    assert(T.dimension() == 3);
    /*res.out_log << "Is valid : " << T.is_valid() << std::endl;*/
    res.out_log << "Number of vertices in RT : " << T.number_of_vertices() << std::endl;
    res.out_log << "Number of cells : " << T.number_of_cells() << std::endl;
    res.out_log << "Number of finite cells : " << T.number_of_finite_cells() << std::endl;
    /*res.out_log << "Inf. vert. atom id : " << T.infinite_vertex()->info().atom_id << std::endl;*/

    // Add finite cells having sqr_r > 0 as graph verticies
    Finite_cells_iterator cit, cit_end = T.finite_cells_end();
    for (cit = T.finite_cells_begin(); cit != cit_end; cit++) {
        if (cit->weighted_circumcenter().weight() > 0)       // cached weighted circumcenter computation
            cit->id() = add_vertex(Cell_handle(cit), G);     // store corresponding graph vertex descriptor in cell's id()
    }

    // store one of infinite cells in graph: it will represent all of them. It should be the last vertex in graph with maximal vertex_descriptor.
    const GVertex inf_graph_v = add_vertex(T.infinite_cell(), G);
    res.out_log << "inf_graph_v descriptor : " << inf_graph_v << std::endl;

    // Add graph edges connecting cells with void
    for (boost::tie(vi, vi_end) = vertices(G); vi != vi_end; ++vi) {
        const GVertex v1 = *vi;           // current graph vertex

        if (v1 != inf_graph_v) { // check edges only from finite cells
            const Cell_handle c1 = G[v1];    // current cell

            for (int i = 0; i < 4; i++) {       // for each of four neighbor cells
                const Cell_handle c2 = c1->neighbor(i);
                bool ends_in_void = false;   // do both corresponding Voronoi verticies lie in void space?
                bool on_same_side = false;   // do both ends lie on the same side of corresponding cell facet?
                const Weighted_point &p1 = c1->vertex((i+1)&3)->point();     // weighted points of common facet between c and c2 cells
                const Weighted_point &p2 = c1->vertex((i+2)&3)->point();
                const Weighted_point &p3 = c1->vertex((i+3)&3)->point();
                const Point &wcc = c1->weighted_circumcenter().point();
                GVertex v2;  // neighbor vertex in graph to construct edge to

                if (c2->id() != -1) {    // neighbor cell is finite and in graph: id()=vertex_descriptor
                    v2 = c2->id();
                    ends_in_void = true;
                    on_same_side = (orientation(p1.point(), p2.point(), p3.point(), wcc) ==
                                    orientation(p1.point(), p2.point(), p3.point(), c2->weighted_circumcenter().point()));
                } else if (T.is_infinite(c2))   {
                    v2 = inf_graph_v;    // use our infinite graph vertex as destination if neighbor cell is infinite
                    ends_in_void = true;
                    on_same_side = (orientation(p1.point(), p2.point(), p3.point(), wcc) !=
                                    orientation(p1.point(), p2.point(), p3.point(), c1->vertex(i)->point().point()));  // same side with infinite vertex = different sides with cell vertex
                }

                typename K::Compute_squared_radius_smallest_orthogonal_sphere_3 r_mouth;  // facet bottleneck squared radius (facet closed if < 0)
                if (ends_in_void && (v2 > v1) && (on_same_side || r_mouth(p1, p2, p3) > 0))
                    add_edge(v1, v2, G);  // use ordering v2>v1 to rule out most of the parallel edges: rely on v_inf>v1 for all finite vertices
            }
        }
    }

    res.out_log << "Number of vertices in G : " << num_vertices(G) << std::endl;
    res.out_log << "Number of edges in G : " << num_edges(G) << std::endl;

    // Find connected components on this Voronoi subnetwork
    std::vector<int> component(num_vertices(G));
    int num_components = boost::connected_components(G, &component[0]);
    const int inf_comp_id = component[inf_graph_v];  // component of infinite vertex
    res.out_log << "Number of connected components : " << num_components << std::endl;
    res.out_log << "Component id of infinite vertex : " << inf_comp_id << std::endl;

    // Reserve storage for results
    res.voids.resize(num_components);
    std::fill(res.atom_surf.begin(), res.atom_surf.end(), 0.0L);

    // Calculate total volume of each finite component
    for (boost::tie(vi, vi_end) = vertices(G); vi != vi_end; ++vi) {
        const int comp_id = component[*vi];
        if (comp_id != inf_comp_id) {   // skip infinite component
            Cell_handle cell = G[*vi];  // current cell
            double Vc, Sc;              // cell volume and surface
            Array_double_4 Sa;          // per atom surface in cell
            Vc = cell_void_volume(cell, Sc, Sa);     // void volume of cell and its surface area

            res.voids[comp_id].volume += Vc;  // add to volume of component
            res.voids[comp_id].surface += Sc;
            for (int k = 0; k < 4; k++) {  // process four cell atoms
                int atom_id = cell->vertex(k)->info().atom_id;
                res.atom_surf[atom_id] += Sa[k];     // atom exposed surface
                res.voids[comp_id].atoms.insert(atom_id);  // add to set of cavity atoms
            }
        }
    }

    // Remove skipped infinite component with zero values
    res.voids.erase(res.voids.begin() + inf_comp_id);

    // Sort cavities by volume
    std::sort(res.voids.begin(), res.voids.end(), void_greater);

    return true;
}


// Creates weighted points from atoms
void build_input(const CavConfig::Atoms_container &in_atoms, const CavConfig::Radii_container &in_radii,
                 double r_scale, double r_probe, Wpi_container &out_points)
{
    for (std::size_t atom_id = 0; atom_id < in_atoms.size(); atom_id++) {
        const Array_double_3 &x = in_atoms[atom_id];
        const Point p(x[0], x[1], x[2]);
        Point_info info = { static_cast<int>(atom_id) };
        Weight weight = CGAL::square(in_radii[atom_id] * r_scale + r_probe);

        out_points.push_back(std::make_pair(Weighted_point(p, weight), info));
    }
}


// Scales coordinates for atoms forming cavity (for drawing that cavity)
void scale_cavity_atoms(const std::set<int> &cav_atoms, const CavConfig::Atoms_container &in_atoms,
                        const CavConfig::Radii_container &in_radii, const Array_double_3 &box,
                        double r_scale, double r_probe,
                        std::vector<Weighted_point> &out_points)
{
    const double ANGS_PER_NM = 1;  // since we output raw graphics for VMD (angst) from Gromacs file (nm)
    for (std::set<int>::const_iterator it = cav_atoms.begin(); it != cav_atoms.end(); ++it) {
        int atom_id = *it;

        Array_double_3 x = in_atoms[atom_id];
        Weight weight = CGAL::square((in_radii[atom_id] * r_scale + r_probe)*ANGS_PER_NM);
#ifdef CORRECT_PBC
        const Array_double_3 &x0 = in_atoms[*cav_atoms.begin()]; // cluster around first atom
        for (int k = 0; k < 3; k++) {
            double dx = x[k] - x0[k];
            x[k] = x[k] - floor(dx / box[k] + 0.5) * box[k];
        }
#endif

        const Point corrected_p(x[0]*ANGS_PER_NM, x[1]*ANGS_PER_NM, x[2]*ANGS_PER_NM);
        out_points.push_back(Weighted_point(corrected_p, weight));
    }
}


bool process_conf(CavConfig &cfg)
{
    Wpi_container points;
    Voids_result res(cfg.atoms.size(), cfg.out_inf);  // reserve space for atoms surf,  pass main log stream

    build_input(cfg.atoms, cfg.radii, cfg.r_scale, cfg.r_probe, points);
    if (!regular_triangulation_voids(points, res)) {
        cfg.out_inf << "algorithm fatal error!" << std::endl;
        return false;
    }

    Surf::write_header(cfg.out_stl);
    cfg.out_inf << std::setprecision(12);
    cfg.out_inf << "Void num :     Volume              Surface " << std::endl;

    for (std::size_t i = 0; i < res.voids.size(); i++) {
        cfg.out_inf << std::setw(8) << i+1 << " : " << std::left << std::setw(17) << res.voids[i].volume << "  ";
        cfg.out_inf << res.voids[i].surface << std::right << std::endl;

        if (cfg.out_vol.is_open()) cfg.out_vol << res.voids[i].volume << std::endl;

        if (cfg.out_stl.is_open()) {  // draw surface of void
            std::vector<Weighted_point> wp;
            scale_cavity_atoms(res.voids[i].atoms, cfg.atoms, cfg.radii, cfg.box, cfg.r_scale, cfg.r_probe, wp);
            Surf::write_cavity_surface(cfg.out_stl, wp.begin(), wp.end(), cfg.nSubdiv);
        }
    }

    Surf::write_footer(cfg.out_stl);

    cfg.out_inf << "Voids count : " << res.voids.size() << std::endl;

    long double total_void_volume = 0.0L;
    long double total_void_surface = 0.0L;
    long double total_void_surface_from_atoms = std::accumulate(res.atom_surf.begin(), res.atom_surf.end(), 0.0L);  // for check
    for (const auto &v: res.voids) {        
        total_void_volume += v.volume;
        total_void_surface += v.surface;
    }

    // accumulate per atom surfaces
    std::transform(res.atom_surf.begin(), res.atom_surf.end(), cfg.atom_confs_surf.begin(), cfg.atom_confs_surf.begin(), std::plus<long double>());

    cfg.out_inf << std::setprecision(std::numeric_limits<double>::digits10);
    cfg.out_inf << "voids total surface : " << total_void_surface << std::endl;
    cfg.out_inf << "voids surf atom sum : " << total_void_surface_from_atoms << std::endl;
    cfg.out_inf << "voids total volume  : " << total_void_volume << std::endl;
    cfg.out_inf << std::endl;

    return true;
}


int main(int argc, const char *argv[])
{
    CavConfig cfg;
    const char *run_control_file = "cavity_volumes_fin.inp";

    if (argc > 1)
        run_control_file = argv[1];

    if (!cfg.init(run_control_file)) {
        std::cerr << "Error while initializing from run control file : " << run_control_file << std::endl;
        return 2;
    }

    cfg.out_inf << "Run control file : " << run_control_file << std::endl;

    CGAL::Timer t;
    int processed_cnt = 0;

    cfg.out_inf << std::endl;
    t.start();
    while (cfg.next_timestep()) {
        const Array_double_3 &a = cfg.atoms.back();
        cfg.out_inf << "Number of input atoms : " << cfg.atoms.size() << std::endl;
        cfg.out_inf << "MD info : " << cfg.ts_info << std::endl;
        cfg.out_inf << "Box : [ " << cfg.box[0] << ", " << cfg.box[1] << ", " << cfg.box[2] << " ]" << std::endl;
        cfg.out_inf << "Last atom : " << a[0] << " " << a[1] << " " << a[2] << std::endl;
        if (!process_conf(cfg)) {
            cfg.out_inf << "process_conf() error. Exiting..." << std::endl;
            return 1;
        }
        processed_cnt++;
    }
    t.stop();

    // save accumulated per-atom surfaces
    if (cfg.out_asf.is_open()) {
        long double total_surf = std::accumulate(cfg.atom_confs_surf.begin(), cfg.atom_confs_surf.end(), 0.0L);
        cfg.out_asf << total_surf << std::endl;  // first line is total exposed surface of all atoms in all confs
        cfg.out_asf << cfg.atom_confs_surf.size() << std::endl;  // second line is the number of atoms
        for (size_t i = 0; i < cfg.atom_confs_surf.size(); i++)
            cfg.out_asf << cfg.atom_confs_surf[i] << std::endl;
    }

    cfg.out_inf << "Processed " << processed_cnt << " (of " << cfg.traj_ts_cnt() << ") configurations." << std::endl;
    cfg.out_inf << "Time: " << t.time() << " sec." << std::endl;

    return 0;
}
