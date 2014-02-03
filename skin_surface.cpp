// Copyright (C) 2014 Alexey Anikeenko
// Copyright (C) 1997 Utrecht University (The Netherlands), ETH Zurich (Switzerland),
// INRIA Sophia-Antipolis (France), Max-Planck-Institute Saarbruecken (Germany),
// and Tel-Aviv University (Israel).
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
#include <CGAL/Union_of_balls_3.h>
#include <CGAL/Skin_surface_3.h>
#include <CGAL/mesh_union_of_balls_3.h>
#include <CGAL/mesh_skin_surface_3.h>
#include <CGAL/subdivide_union_of_balls_mesh_3.h>
#include <CGAL/subdivide_skin_surface_mesh_3.h>
#include <CGAL/Polyhedron_3.h>
#include <vector>
#include <list>
#include <map>
#include <fstream>
#include "skin_surface.h"

namespace Surf {

enum { UNTAGGED_COMPONENT = -1 };

// Finds next vertex which is not yet tagged
template<class Polyhedron>
typename Polyhedron::Vertex_handle get_untagged_vertex(Polyhedron &p, std::map<typename Polyhedron::Vertex_handle, int> &comp)
{
    for (typename Polyhedron::Vertex_iterator it = p.vertices_begin(); it != p.vertices_end(); it++) {
        if (comp[it] == UNTAGGED_COMPONENT)
            return it;
    }
    return NULL;
}


// Tags connected component with id.
// return the size (number of vertices) of the component.
// Adapted from CGAL/HalfedgeDS_decorator.h
template<class Polyhedron>
unsigned int tag_component(typename Polyhedron::Vertex_handle pSeedVertex,
                           std::map<typename Polyhedron::Vertex_handle, int>& comp, int id)
{
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Halfedge_around_vertex_circulator Halfedge_around_vertex_circulator;
    unsigned int number_of_vertices = 0; // size (number of vertices) of the component

    std::list<Vertex_handle> vertices;
    vertices.push_front(pSeedVertex);
    while (!vertices.empty()) {
        Vertex_handle pVertex = vertices.front();
        vertices.pop_front();

        // Skip vertex if already done
        if (comp[pVertex] >= 0)
            continue;

        // Mark vertex done
        comp[pVertex] = id;
        number_of_vertices++;

        // Add vertex's "free" neighbors to the list
        Halfedge_around_vertex_circulator neighbor_cir, neighbor_end;
        neighbor_cir = pVertex->vertex_begin();
        neighbor_end = neighbor_cir;
        CGAL_For_all(neighbor_cir, neighbor_end) {
            Vertex_handle neighbor = neighbor_cir->opposite()->vertex();
            if (comp[neighbor] == UNTAGGED_COMPONENT)
                vertices.push_front(neighbor);
        }
    }

    return number_of_vertices;
}


// Computes all connected components of the polyhedron and erases the largest.
template<class Polyhedron>
void erase_largest_connected_components(Polyhedron &p)
{
    typedef typename Polyhedron::Vertex_handle Vertex_handle;
    typedef typename Polyhedron::Vertex_iterator Vertex_iterator;
    std::map<Vertex_handle, int> comp;
    std::vector<Vertex_handle> seeds;

    // Tag all mesh vertices as UNTAGGED_COMPONENT
    for (Vertex_iterator it = p.vertices_begin(); it != p.vertices_end(); it++) {
        comp[it] = UNTAGGED_COMPONENT;
    }

    Vertex_handle seed_vertex = NULL;
    int num_comps = 0;
    unsigned int small_comp_size = std::numeric_limits<unsigned int>::max();
    unsigned int large_comp_size = 0;
    Vertex_handle small_comp_seed = NULL;
    Vertex_handle large_comp_seed = NULL;

    while ((seed_vertex = get_untagged_vertex(p, comp)) != NULL) {
        // Tag component by id and compute its size (number of vertices)
        unsigned int number_of_vertices = tag_component<Polyhedron>(seed_vertex, comp, num_comps);
        // Update current records of the largest and smallest components
        if (number_of_vertices < small_comp_size) {
            small_comp_size = number_of_vertices;
            small_comp_seed = seed_vertex;
        }
        if (number_of_vertices > large_comp_size) {
            large_comp_size = number_of_vertices;
            large_comp_seed = seed_vertex;
        }
        seeds.push_back(seed_vertex);
        num_comps++;
    }

    if (num_comps != 2)
        std::cerr << " skin mesh components != 2 : " << num_comps << std::endl;

    // remove the largest component
    Vertex_handle vertex = large_comp_seed;
    if (vertex != small_comp_seed && vertex->halfedge() != NULL) // if not isolated vertex
        p.erase_connected_component(vertex->halfedge());

    // Alternative is to remove all components but smallest:
    //for(std::vector<Vertex_handle>::iterator si = seeds.begin(); si != seeds.end(); ++si)
    //  { erase if != smallest_seed }
}


// Writes STL file header
void write_header(std::ostream &out)
{
    out << "solid " << "cavities" << std::endl;
}


// Writes STL file footer
void write_footer(std::ostream &out)
{
    out << "endsolid " << "cavities" << std::endl;
}


// Writes facets of mesh polyhedron to STL file
template<class Polyhedron>
void write_facets_to_stl(const Polyhedron &P, std::ostream &out_stl)
{
    typedef typename Polyhedron::Facet_const_iterator Facet_const_iterator;
    //for all triangular facets
    for (Facet_const_iterator fi = P.facets_begin(); fi != P.facets_end(); ++fi) {
        typename Polyhedron::Halfedge_const_handle h = fi->halfedge();
        if (h->next()->next()->next() != h) {
            std::cerr << "Polyhedron facet is not triangle." << std::endl;
            return;
        }
        K::Point_3 p = h->vertex()->point();
        K::Point_3 q = h->next()->vertex()->point();
        K::Point_3 r = h->next()->next()->vertex()->point();

        K::Vector_3 n = CGAL::cross_product(q - p, r - p);
        n = n / CGAL::sqrt(n * n);

        out_stl << "    facet normal " << n << std::endl;
        out_stl << "      outer loop " << std::endl;
        out_stl << "        vertex " << p << std::endl;
        out_stl << "        vertex " << q << std::endl;
        out_stl << "        vertex " << r << std::endl;
        out_stl << "      endloop " << std::endl;
        out_stl << "    endfacet " << std::endl;
    }
}


// Constructs and writes 3D mesh of the cavity formed by the given weighted points.
// Based on CGAL's 3D Skin Surface Meshing.
void write_cavity_surface(std::ostream &out, WP_iterator begin, WP_iterator end, int nSubdiv)
{
    typedef CGAL::Union_of_balls_3<Traits> Union_of_balls_3;
    typedef CGAL::Polyhedron_3<K, CGAL::Skin_surface_polyhedral_items_3<Union_of_balls_3> > Polyhedron;
    Polyhedron p;

    Union_of_balls_3 union_of_balls(begin, end, Traits(), false);
    CGAL::mesh_union_of_balls_3(union_of_balls, p);
    
    // FIXME: assumes that the mesh of the union of balls with one cavity inside 
    // has two connected components and the largest is outer surface
    erase_largest_connected_components(p);

    CGAL::subdivide_union_of_balls_mesh_3(union_of_balls, p, nSubdiv);  // refine the mesh
    write_facets_to_stl(p, out);
}
} // namespace Surf
