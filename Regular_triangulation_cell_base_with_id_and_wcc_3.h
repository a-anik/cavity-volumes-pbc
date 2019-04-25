// Copyright (c) 1999-2006  INRIA Sophia-Antipolis (France).
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

#ifndef CV_REGULAR_TRIANGULATION_CELL_BASE_WITH_ID_AND_WCC_3_H
#define CV_REGULAR_TRIANGULATION_CELL_BASE_WITH_ID_AND_WCC_3_H

#include <CGAL/license/Triangulation_3.h>

#include <CGAL/basic.h>
#include <CGAL/triangulation_assertions.h>
#include <CGAL/Regular_triangulation_cell_base_3.h>

#include <CGAL/constructions/kernel_ftC3.h>

namespace CGAL {

// Cell of a regular triangulation of any dimension <=3,
// with int_id info, and lazy computation of its weighted circumcenter.
// Based on CGAL/Triangulation_cell_base_with_info_3.h,
//          CGAL/Triangulation_cell_base_with_circumcenter_3.h,
//          CGAL/Regular_triangulation_cell_base_3.h

template < typename GT, typename Cb = Regular_triangulation_cell_base_3<GT> >
class Regular_triangulation_cell_base_with_id_and_wcc_3
  : public Cb
{
  typedef GT                                           Geom_traits;
  typedef typename Geom_traits::FT                     FT;
  typedef typename Geom_traits::Weighted_point_3       Weighted_point;
  typedef typename Weighted_point::Point               Bare_point;
  Weighted_point _wcc;             // weighted cirumcenter
  bool _wcc_constructed;
  int  _id;                        // cell info

  void invalidate_circumcenter() { _wcc_constructed = false; }

public:
  typedef typename Cb::Vertex_handle                   Vertex_handle;
  typedef typename Cb::Cell_handle                     Cell_handle;

  template < typename TDS2 >
  struct Rebind_TDS {
    typedef typename Cb::template Rebind_TDS<TDS2>::Other                Cb2;
    typedef Regular_triangulation_cell_base_with_id_and_wcc_3<GT, Cb2>   Other;
  };

  Regular_triangulation_cell_base_with_id_and_wcc_3()
    : Cb(), _wcc_constructed(false), _id(-1) {}

  Regular_triangulation_cell_base_with_id_and_wcc_3(Vertex_handle v0,
                                                    Vertex_handle v1,
                                                    Vertex_handle v2,
                                                    Vertex_handle v3)
    : Cb(v0, v1, v2, v3), _wcc_constructed(false), _id(-1) {}

  Regular_triangulation_cell_base_with_id_and_wcc_3(Vertex_handle v0,
                                                    Vertex_handle v1,
                                                    Vertex_handle v2,
                                                    Vertex_handle v3,
                                                    Cell_handle   n0,
                                                    Cell_handle   n1,
                                                    Cell_handle   n2,
                                                    Cell_handle   n3)
    : Cb(v0, v1, v2, v3, n0, n1, n2, n3), _wcc_constructed(false), _id(-1) {}


  // We must override the functions that modify the vertices.
  // And if the point inside a vertex is modified, we fail,
  // but there's not much we can do for this now.
  void set_vertex(int i, Vertex_handle v)
  {
      invalidate_circumcenter();
      Cb::set_vertex(i, v);
  }

  void set_vertices()
  {
      invalidate_circumcenter();
      Cb::set_vertices();
  }

  void set_vertices(Vertex_handle v0, Vertex_handle v1,
                    Vertex_handle v2, Vertex_handle v3)
  {
      invalidate_circumcenter();
      Cb::set_vertices(v0, v1, v2, v3);
  }

  const Weighted_point& weighted_circumcenter()
  {
      if (!_wcc_constructed) {
          const Weighted_point & p = this->vertex(0)->point();
          const Weighted_point & q = this->vertex(1)->point();
          const Weighted_point & r = this->vertex(2)->point();
          const Weighted_point & s = this->vertex(3)->point();
          FT x, y, z, w;
          weighted_circumcenterC3(p.x(), p.y(), p.z(), p.weight(),
                                  q.x(), q.y(), q.z(), q.weight(),
                                  r.x(), r.y(), r.z(), r.weight(),
                                  s.x(), s.y(), s.z(), s.weight(),
                                  x, y, z, w);
          _wcc = Weighted_point(Bare_point(x, y, z), w);
      }
      return _wcc;
  }

  int id() const { return _id; }

  int & id() { return _id; }
};

}  // namespace CGAL

#endif /* CV_REGULAR_TRIANGULATION_CELL_BASE_WITH_ID_AND_WCC_3_H */
