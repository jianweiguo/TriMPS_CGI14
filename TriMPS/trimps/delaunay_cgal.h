/*
 *  _____ _____ ________  _
 * /  __//  __//  __/\  \//
 * | |  _|  \  |  \   \  / 
 * | |_//|  /_ |  /_  /  \ 
 * \____\\____\\____\/__/\\
 *
 * Graphics Environment for EXperimentations.
 *  Copyright (C) 2006 INRIA - Project ALICE
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
 *
 *  If you modify this software, you should include a notice giving the
 *  name of the person performing the modification, the date of modification,
 *  and the reason for such modification.
 *
 *  Contact: 
 *
 *     ALICE Project - INRIA
 *     INRIA Lorraine, 
 *     Campus Scientifique, BP 239
 *     54506 VANDOEUVRE LES NANCY CEDEX 
 *     FRANCE
 *
 *  Note that the GNU General Public License does not permit incorporating
 *  the Software into proprietary programs. 
 */

#ifndef __DELAUNAY_CGAL__
#define __DELAUNAY_CGAL__

#include "geometry.h"

//______________________________________________________________________________________
// CGAL stuff

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_euclidean_traits_2.h>
#include <CGAL/Regular_triangulation_filtered_traits_2.h>
#include <CGAL/Regular_triangulation_2.h>
#include <CGAL/Triangulation_vertex_base_with_info_2.h>
#include <CGAL/Lipschitz_sizing_field_2.h>

namespace Geex {


	// ------------- My stuff to extend CGAL --------------------------

	class VertexDecoration {
	public:
		VertexDecoration()
			:index(-1), dual_intersects_boundary(false), locked(false), radius(0), type(0)
		{}
	public:
		int    index ;
		bool   dual_intersects_boundary ;
		bool   locked ;
		double radius ;
		int    type ; // used for plant simulation
	} ;

	class CellDecoration {
	public:
		CellDecoration() :infinite(false), dual_outside(false), dual(0, 0), has_gap(false), visited(false), gap_index(-1)
		{}
	public:
		bool infinite ;
		bool dual_outside ;
		vec2 dual ;
		bool has_gap ;
		int  gap_index ;
		bool visited ;
	} ;


	// ------------- CGAL stuff ---------------------------------

	// ----------------------- A CGAL::Vertex with decoration ------------------
	template < class Gt, class Vb = CGAL::Regular_triangulation_vertex_base_2<Gt> >
	class Vertex : public  Vb, public VertexDecoration  {
		typedef Vb superclass;
	public:
		typedef typename Vb::Vertex_handle      Vertex_handle;
		typedef typename Vb::Face_handle        Face_handle;
		typedef typename Vb::Point              Point;
        
		template < typename TDS2 >
		struct Rebind_TDS {
			typedef typename Vb::template Rebind_TDS<TDS2>::Other Vb2;
			typedef Vertex<Gt,Vb2> Other;
		} ;
        
	public:
		Vertex() : superclass() {}
		Vertex(const Point & p) : superclass(p) {}
		Vertex(const Point & p, Face_handle f) : superclass(f,p) {}
		Vertex(Face_handle f) : superclass(f) {}
	} ;


	// ----------------------- A CGAL::Cell with decoration ------------------
	template < class Gt, class Cb = CGAL::Regular_triangulation_face_base_2<Gt> >
	class Cell : public Cb, public CellDecoration {
		typedef Cb superclass;
	public:
		typedef typename Cb::Vertex_handle      Vertex_handle;
		typedef typename Cb::Face_handle        Face_handle;
		template < typename TDS2 >
		struct Rebind_TDS {
			typedef typename Cb::template Rebind_TDS<TDS2>::Other Cb2;
			typedef Cell<Gt,Cb2> Other;
		};


		Cell() : superclass() {  }
		Cell(
			Vertex_handle v0, Vertex_handle v1, Vertex_handle v2
		) : superclass(v0,v1,v2) { }
            
		Cell(
			Vertex_handle v0, Vertex_handle v1, Vertex_handle v2, 
			Face_handle n0, Face_handle n1, Face_handle n2 
		) : superclass(v0,v1,v2,n0,n1,n2) { }

	} ;

	typedef CGAL::Exact_predicates_inexact_constructions_kernel K;
	typedef K::Point_2 Point_2;
	typedef K::Point_2 Circle_2;
	typedef K::Iso_rectangle_2 Iso_rectangle_2;

	typedef double Weight;
	typedef CGAL::Regular_triangulation_euclidean_traits_2<K,Weight>  Gt;
	typedef Vertex<Gt> Vb ;
	typedef Cell<Gt>   Cb ;
	typedef CGAL::Triangulation_data_structure_2<Vb,Cb> TDS;
	typedef CGAL::Regular_triangulation_2<Gt, TDS> Regular;


    typedef K::Point_2    Point;
    typedef K::Vector_2   CGAL_Vector;
	typedef Gt::Weighted_point_2 WtPoint ;
	typedef K::Ray_2      Ray_2;
	typedef K::Vector_2   Vector_2;
	typedef K::Segment_2  Segment_2 ;
	typedef K::Line_2     Line_2 ;

    class DelaunayBase : public Regular {
    public:
    } ;

	// a cluster of faces belong to a connected gap
	class GapFaceCluster : public std::vector<DelaunayBase::Face_handle> {
	public:
		GapFaceCluster(bool is_boundary=false) : is_boundary_(is_boundary) {
		}
	public:
		bool is_boundary_ ;
	} ;

    inline Point to_cgal(const vec2& p) { return Point(p.x, p.y) ; }
	inline vec2  to_geex(const Point& p) { return vec2(p.x(), p.y()) ; }
	inline vec2  center(DelaunayBase::Face_handle f) {
		return (to_geex(f->vertex(0)->point()) + to_geex(f->vertex(1)->point()) + to_geex(f->vertex(2)->point()))/3.0 ;
	}

	/*
	 * sizing field
	**/
	typedef CGAL::Lipschitz_sizing_field_2<K> Lipschitz_sizing_field;
}


#endif
