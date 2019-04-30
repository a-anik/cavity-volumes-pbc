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

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>
#include <CGAL/Periodic_3_triangulation_traits_3.h>
#include <CGAL/Periodic_3_offset_3.h>
#include <CGAL/Random.h>
#include <CGAL/Timer.h>
#include <CGAL/point_generators_3.h>
#include <boost/array.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/connected_components.hpp>

#include <cassert>
#include <vector>
#include <set>
#include <map>
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
#include "config_file.h"

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;

typedef CGAL::Periodic_3_offset_3                           Offset;
typedef struct { int atom_id;  Offset off; }                Point_info;
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

typedef CGAL::Periodic_3_triangulation_traits_3<K>          PTraits;  // used for offset -> points construction
typedef PTraits::Iso_cuboid_3                               Iso_cuboid;

typedef Rt::Finite_cells_iterator                           Finite_cells_iterator;
typedef Rt::Cell_handle                                     Cell_handle;
typedef boost::array<int, 4>                                Array_int_4;
typedef boost::array<Offset, 3>                             Array_off_3;
typedef std::pair<Array_int_4, Array_off_3>                 Cell_key;  // sorted set of cell atoms and three relative periodic offsets
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
    double domain_cells_volume;           // sum of primary domain tetrahedral cells (for checking)
    std::vector<long double> atom_surf;   // for each atom : surface exposed to cavities (natoms values)
    std::ostream &out_log;                // information stream (debugging Voronoi calculation info)
    Voids_result(int natoms, std::ostream &strm) : atom_surf(natoms), out_log(strm) {}
};

bool void_greater(const Voids_result::Void_info& v1, const Voids_result::Void_info& v2) { return v1.volume > v2.volume; }

// Cell key: two cells in periodic triangulation are equivalent if they can be superimposed by translation:
//   1) set(a0,a1,a2,a3) == set(a0',a1',a2',a3')                         have same set of atom ids
//   2) o(a1)-o(a1') == o(a2)-o(a2') == o(a3)-o(a3') == o(a4)-o(a4')     offsets of the corresponding atoms have identical translation
// Three offsets relative to the first vertex are used as a key after sorting them according to atom_ids and lexicographically by offset.
Cell_key cell_periodic_key(const Cell_handle c)
{
    boost::array<std::pair<int, Offset>, 4> a;   // pairs (a_i, o_i) for sorting them together

    for (int i = 0; i < 4; i++) {
        const Point_info &info = c->vertex(i)->info();
        a[i] = std::make_pair(info.atom_id, info.off);
    }
    std::sort(a.begin(), a.end());

    // periodic translation doesnt change lexicographical order of offsets
    // and their order matters when some atoms in cell have same ids
    const Offset o0 = a[0].second;  // offset of first atom
    Array_int_4 sorted_atom_ids = {{a[0].first, a[1].first, a[2].first, a[3].first}};
    Array_off_3 rel_offsets = {{a[1].second - o0, a[2].second - o0, a[3].second - o0}};

    return std::make_pair(sorted_atom_ids, rel_offsets);
}

inline std::ostream &operator<<(std::ostream &os, const Cell_key &k)
{
    os << "(" << k.first[0] << " " << k.first[1] << " " << k.first[2] << " " << k.first[3] << ") ";
    os << "[" << k.second[0] << "] [" << k.second[1] << "] [" << k.second[2] << "]";
    return os;
}

inline std::ostream &operator<<(std::ostream &os, const Point_info &inf)
{
    os << inf.atom_id << " [" << inf.off << "] ";
    return os;
}


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
// Return value: true if converged (periodic voids map can be built)
bool periodic_regular_triangulation_voids(const Wpi_container& points, Voids_result &res)
{
    typedef typename boost::adjacency_list<boost::vecS, boost::vecS, boost::undirectedS, Cell_handle> CellGraph;
    typedef typename boost::graph_traits<CellGraph>::vertex_descriptor GVertex;
    typedef typename boost::graph_traits<CellGraph>::vertex_iterator GVertex_iterator;
    typedef typename std::map<Cell_key, Cell_handle> Cell_periodic_map;
    enum { NOT_FROM_T3 = -1, COVERED_PRIMARY_CELL = -2, COPY_WITH_VERTEX_IN_DOMAIN = -3 };  // initially id()=NOT_FROM_T3

    Rt T;            // regular triangulation
    CellGraph G;     // Voronoi subgraph where cells_sqr_r > 0 and edges_sqr_r > 0 (cells and edges with voids)
    GVertex_iterator vi, vi_end;
    Cell_periodic_map map_cells;  // map periodic copies of cell to single cell representaion in primary domain

    res.out_log << "Number of input points for RT : " <<  points.size() << std::endl;
    T.insert(points.begin(), points.end());    // insert all points in a row (this is faster than one insert() at a time).
    T.infinite_vertex()->info().atom_id = -1;  // set special atom_id for dummy infinite triangulation vertex

    assert(T.dimension() == 3);
    /*res.out_log << "Is valid : " << T.is_valid() << std::endl;*/
    res.out_log << "Number of vertices in RT : " << T.number_of_vertices() << std::endl;
    res.out_log << "Number of cells : " << T.number_of_cells() << std::endl;
    res.out_log << "Number of finite cells : " << T.number_of_finite_cells() << std::endl;

    // Map periodic cell key to one of cell representations - primary domain cell
    Finite_cells_iterator cit, cit_end = T.finite_cells_end();
    for (cit = T.finite_cells_begin(); cit != cit_end; cit++) {
        bool has_vertex_in_domain = false;
        for (int i = 0; i < 4; i++) {
            const Offset &o = cit->vertex(i)->info().off;
            has_vertex_in_domain |= (o[0] == 0 && o[1] == 0 && o[2] == 0);
        }

        if (has_vertex_in_domain) {   // cell have at least one vertex in primary domain
            if (map_cells.insert(std::make_pair(cell_periodic_key(cit), cit)).second) {
                // cell discovered for the first time - it will be primary representation
                if (cit->weighted_circumcenter().weight() > 0)       // add as graph vertex if sqr_r > 0 (not fully covered cell)
                    cit->id() = add_vertex(Cell_handle(cit), G);     // store corresponding graph vertex descriptor in cell's id()
                else
                    cit->id() = COVERED_PRIMARY_CELL;                // mark as fully covered otherwise
            } else {
                cit->id() = COPY_WITH_VERTEX_IN_DOMAIN;
            }
        }
    }

    res.out_log << "Number of T3 periodic triangulation cells : " << map_cells.size() << std::endl;

    // Add edges connecting cells with void
    for (boost::tie(vi, vi_end) = vertices(G); vi != vi_end; ++vi) {
        const GVertex v1 = *vi;           // current graph vertex
        const Cell_handle c1 = G[v1];     // current cell, it's primary and with void

        for (int i = 0; i < 4; i++) {     // for each neighbor cell
            const Cell_handle c2 = c1->neighbor(i);
            Cell_handle c2_repr = c2;  // primary copy of neighbor cell: init as c2 for safety
            Cell_periodic_map::iterator pos;

            // Find copy of neighbour cell in primary domain
            // TODO: For very small/sparse/inhomogeneous systems requirement of existance of complete
            //       periodic copy of the cell is too strong even for 27-sheeted covering [Dolbilin97, Caroli09].
            //       There we could additionally check all facets of primary domain in order to find periodic neighbour.
            if (c2->id() >= 0 || c2->id() == COVERED_PRIMARY_CELL)  // if neig is in graph or marked covered => it's primary
                c2_repr = c2;
            else if ((pos = map_cells.find(cell_periodic_key(c2))) != map_cells.end())  // else try to find it in cells map
                c2_repr = pos->second;
            else
                return false;  // periodic image of neigbor cell not found, T3 triangulation did not converge

            if (c2_repr->id() >= 0) {  // if neighbor cell is also in void
                const GVertex v2 = c2_repr->id();   // neighbor vertex in graph to construct edge to
                const Weighted_point &p1 = c1->vertex((i+1)&3)->point();   // weighted points of common facet (c1,c2)
                const Weighted_point &p2 = c1->vertex((i+2)&3)->point();
                const Weighted_point &p3 = c1->vertex((i+3)&3)->point();
                const Point &wcc1 = c1->weighted_circumcenter().point();
                const Point &wcc2 = c2->weighted_circumcenter().point();
                bool on_same_side = (orientation(p1.point(), p2.point(), p3.point(), wcc1) ==  // do both ends lie on the same side
                                     orientation(p1.point(), p2.point(), p3.point(), wcc2));   // of corresponding cell facet?

                typename K::Compute_squared_radius_smallest_orthogonal_sphere_3 r_mouth;  // facet bottleneck squared radius
                if ((v2 > v1) && (on_same_side || r_mouth(p1, p2, p3) > 0))
                    add_edge(v1, v2, G);  // use ordering v2>v1 to rule out parallel edges
            }
        }
    }

    res.out_log << "Number of vertices in G : " << num_vertices(G) << std::endl;
    res.out_log << "Number of edges in G : " << num_edges(G) << std::endl;

    // Find connected components on this Voronoi subnetwork
    std::vector<int> component(num_vertices(G));
    int num_components = boost::connected_components(G, &component[0]);
    res.out_log << "Number of connected components (voids) : " << num_components << std::endl;

    // Reserve storage for results
    res.voids.resize(num_components);
    std::fill(res.atom_surf.begin(), res.atom_surf.end(), 0.0L);

    // Calculate total volume and surface of each component
    for (boost::tie(vi, vi_end) = vertices(G); vi != vi_end; ++vi) {
        const int comp_id = component[*vi];
        Cell_handle cell = G[*vi];               // current cell
        double Vc, Sc;                           // cell volume and surface
        Array_double_4 Sa;                       // per atom surface in cell
        Vc = cell_void_volume(cell, Sc, Sa);     // void volume of cell and its surface area

        res.voids[comp_id].volume += Vc;               // add to volume of component
        res.voids[comp_id].surface += Sc;
        for (int k = 0; k < 4; k++) {            // process four cell atoms
            int atom_id = cell->vertex(k)->info().atom_id;
            res.atom_surf[atom_id] += Sa[k];     // atom exposed surface
            res.voids[comp_id].atoms.insert(atom_id);  // add to set of cavity atoms
        }
    }

    res.domain_cells_volume = 0;  // sum volumes of primary domain cells. Should equal to box volume.
    for (Cell_periodic_map::iterator pos = map_cells.begin(); pos != map_cells.end(); ++pos) {
        res.domain_cells_volume += fabs(T.tetrahedron(pos->second).volume());
    }

    // Sort cavities by volume
    std::sort(res.voids.begin(), res.voids.end(), void_greater);

    return true;
}


// Builds periodically thickened configuraiton by adding periodic copies of atoms by the sides
int build_pbc_thickened_input(const CavConfig::Atoms_container &in_atoms, const CavConfig::Radii_container &in_radii,
                              const Array_double_3 &box, double shell_width, double r_scale, double r_probe, Wpi_container &out_points)
{
    int num_pbc_corrected = 0;
    Iso_cuboid domain(0, 0, 0, box[0], box[1], box[2]);
    PTraits ptraits = PTraits(domain);
    PTraits::Construct_point_3 periodic_point_construct = ptraits.construct_point_3_object();

    for (size_t atom_id = 0; atom_id < in_atoms.size(); atom_id++) {
        bool pbc_corrected = false;
        const Array_double_3 &atom = in_atoms[atom_id];
        Point_info info = { static_cast<int>(atom_id), Offset(0, 0, 0) };
        Weight weight = CGAL::square(in_radii[atom_id] * r_scale + r_probe);  // weight for regular triangulation is square radius

        // correct coordinates to assure that atom locates in primary domain
        K::FT x[3];
        for (int k = 0; k < 3; k++) {
            x[k] = atom[k];
            if (x[k] < 0 || x[k] >= box[k]) {
                x[k] = x[k] - floor(x[k] / box[k]) * box[k];
                pbc_corrected = true;
            }
        }
        if (pbc_corrected) num_pbc_corrected++;

        const Point corrected_p(x[0], x[1], x[2]);  // PBC corrected original point
        out_points.push_back(std::make_pair(Weighted_point(corrected_p, weight), info));

        // also add periodic images of atom if they fall within box*shell_width near primary domain
        for (int i = -1; i <= 1; i++)
            for (int j = -1; j <= 1; j++)
                for (int k = -1; k <= 1; k++) {
                    if (i == 0 && j == 0 && k == 0) continue;
                    const Point virt_point = periodic_point_construct(corrected_p, Offset(i, j, k));
                    info.off = Offset(i, j, k);
                    Weighted_point_with_info wpi = std::make_pair(Weighted_point(virt_point, weight), info);
                    if (fabs(virt_point[0]/box[0] - 0.5) <= shell_width + 0.5 &&
                        fabs(virt_point[1]/box[1] - 0.5) <= shell_width + 0.5 &&
                        fabs(virt_point[2]/box[2] - 0.5) <= shell_width + 0.5)
                    {
                         out_points.push_back(wpi);
                    }
                }
    }

    return num_pbc_corrected;
}


// Calculates periodically unwrapped coordinates for atoms forming cavity (for drawing that cavity)
void pbc_unwrap_cavity_atoms(const std::set<int> &cav_atoms, const CavConfig::Atoms_container &in_atoms,
                             const CavConfig::Radii_container &in_radii, const Array_double_3 &box,
                             double r_scale, double r_probe,
                             std::vector<Weighted_point> &out_points)
{
    const double ANGS_PER_NM = 10;  // since we output raw graphics for VMD (angst) from Gromacs file (nm)
    for (std::set<int>::const_iterator it = cav_atoms.begin(); it != cav_atoms.end(); ++it) {
        int atom_id = *it;

        Array_double_3 x = in_atoms[atom_id];
        const Array_double_3 &x0 = in_atoms[*cav_atoms.begin()];  // cluster around first atom
        Weight weight = CGAL::square((in_radii[atom_id] * r_scale + r_probe)*ANGS_PER_NM);

        for (int k = 0; k < 3; k++) {
            double dx = x[k] - x0[k];
            x[k] = x[k] - floor(dx / box[k] + 0.5) * box[k];
            /*if( dx > box[k]/2 ) x[k] = x[k] -  box[k];
              else if( dx < -box[k]/2 ) x[k] = x[k] + box[k];*/
        }

        const Point corrected_p(x[0]*ANGS_PER_NM, x[1]*ANGS_PER_NM, x[2]*ANGS_PER_NM);  // PBC corrected original point
        out_points.push_back(Weighted_point(corrected_p, weight));
    }
}


// Performs all calculations on current timestep data, writes results
bool process_conf(CavConfig &cfg)
{
    Wpi_container points;
    Voids_result res(cfg.atoms.size(), cfg.out_inf);  // reserve space for atoms surf,  pass main log stream
    int num_pbc_corrected;

    num_pbc_corrected = build_pbc_thickened_input(cfg.atoms, cfg.radii, cfg.box, cfg.shell_width, cfg.r_scale, cfg.r_probe, points);
    if (num_pbc_corrected)
        cfg.out_inf << "Number of atoms with PBC corrected coordinates : " << num_pbc_corrected << std::endl;

    // Try to build periodic connected empty cells graph.
    // Find void volumes and areas of connected components of that graph.
    if (!periodic_regular_triangulation_voids(points, res)) {
        cfg.out_inf << "Not converged with default thicken parameter : " << cfg.shell_width << std::endl;
        cfg.out_inf << "Trying full 27-sheeted covering..." << std::endl;
        points.clear();
        build_pbc_thickened_input(cfg.atoms, cfg.radii, cfg.box, 1.0, cfg.r_scale, cfg.r_probe, points);
        if (!periodic_regular_triangulation_voids(points, res)) {
            cfg.out_inf << "Couldn't build periodic void cell network even for 27-sheeted covering: " << std::endl;
            cfg.out_inf << "probably model is too small or contains large voids." << std::endl;
            cfg.out_inf << "algorithm fatal error!" << std::endl;
            return false;
        }
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
            pbc_unwrap_cavity_atoms(res.voids[i].atoms, cfg.atoms, cfg.radii, cfg.box, cfg.r_scale, cfg.r_probe, wp);
            Surf::write_cavity_surface(cfg.out_stl, wp.begin(), wp.end(), cfg.nSubdiv);
        }
    }

    Surf::write_footer(cfg.out_stl);

    cfg.out_inf << "Voids count : " << res.voids.size() << std::endl;
    const double box_volume = cfg.box[0] * cfg.box[1] * cfg.box[2];
    cfg.out_inf << "SumV/BoxV : " << res.domain_cells_volume/box_volume << std::endl;
    if (fabs(res.domain_cells_volume/box_volume - 1.0) > 1e-10) {
        cfg.out_inf << "Sum of cells volumes != box_volume : fatal error" << std::endl;
        return false;
    }

    long double total_void_volume = 0.0L;
    long double total_void_surface = 0.0L;
    long double total_void_surface_from_atoms = std::accumulate(res.atom_surf.begin(), res.atom_surf.end(), 0.0L);  // for check
    for (const auto &v: res.voids) {
        total_void_volume += v.volume;
        total_void_surface += v.surface;
    }

    cfg.out_evf << total_void_volume / box_volume << std::endl;
    // accumulate per atom surfaces
    std::transform(res.atom_surf.begin(), res.atom_surf.end(), cfg.atom_confs_surf.begin(), cfg.atom_confs_surf.begin(), std::plus<long double>());

    cfg.out_inf << std::setprecision(std::numeric_limits<double>::digits10);
    cfg.out_inf << "packing fraction : " << 1.0 - total_void_volume/box_volume << std::endl;
    cfg.out_inf << "box volume : " << box_volume << std::endl;
    cfg.out_inf << "voids total surface : " << total_void_surface << std::endl;
    cfg.out_inf << "voids surf atom sum : " << total_void_surface_from_atoms << std::endl;
    cfg.out_inf << "voids total volume  : " << total_void_volume << std::endl;
    cfg.out_inf << std::endl;

    return true;
}


int main(int argc, const char *argv[])
{
    CavConfig cfg;
    const char *run_control_file = "cavity_volumes_pbc.inp";

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
