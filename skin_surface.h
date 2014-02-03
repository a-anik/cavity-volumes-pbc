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

#ifndef CAVITYSKIN_SURFACE_H
#define CAVITYSKIN_SURFACE_H

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Skin_surface_3.h>
#include <iostream>
#include <vector>

// Namespace contains extra functions for 3D visualization of cavities.
namespace Surf {

typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
typedef CGAL::Skin_surface_traits_3<K>                      Traits;
typedef Traits::Weighted_point_3                            Weighted_point;
typedef Weighted_point::Weight                              Weight;
typedef std::vector<Weighted_point>                         WP_container;
typedef WP_container::iterator                              WP_iterator;

void write_header(std::ostream &out);
void write_footer(std::ostream &out);

// Constructs and writes 3D mesh of the cavity formed by the given weighted points
void write_cavity_surface(std::ostream &out, WP_iterator begin, WP_iterator end, int nSubdiv = 3);

} // namespace Surf
#endif /* CAVITYSKIN_SURFACE_H  */
