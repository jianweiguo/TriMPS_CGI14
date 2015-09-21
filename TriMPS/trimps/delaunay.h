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

#ifndef __DELAUNAY__
#define __DELAUNAY__

//#define NOMINMAX 

#include "delaunay_cgal.h"
#include "polygons.h"
#include "gap_primitive.h"
#include "grid.h"

#include <Geex/combinatorics/map.h>
#include <Geex/combinatorics/map_builder.h>
#include <Geex/combinatorics/map_editor.h>
#include <Geex/combinatorics/map_properties.h>

namespace Geex {
	class GapCompute ;
	class GridCell ;
	enum  WeightScheme {WT_UNIFORM, WT_RANDOM, WT_FUNCTION, WT_IMAGE, WT_AIS, WT_LFS} ;
//______________________________________________________________________________________

	class cellPair{
	public:
		cellPair(vec2 sample, vec2 center_sample, double cen_radius, int_two id, int_two center_id) 
			: sample_(sample), center_sample_(center_sample), cen_radius_(cen_radius), id_(id), center_id_(center_id) {
				dist_ = distance2(sample_, center_sample_)- cen_radius_*cen_radius_;
		}

		bool operator > (const cellPair& rhs) const {
			return dist_ > rhs.dist_ ;
		}
		bool operator < (const cellPair& rhs) const {
			return dist_ < rhs.dist_ ;
		}
	public:
		vec2 sample_;
		vec2 center_sample_;
		double cen_radius_;
		int_two id_;
		int_two center_id_;
		double dist_;
	};
//______________________________________________________________________________________________
    class Delaunay : public DelaunayBase {
        typedef DelaunayBase baseclass ;

    public:
        Delaunay() ;
		~Delaunay() ;

        void save(const std::string& filename) ;
        void load(const std::string& filename) ;
		void save_ps(const std::string& filename);
        void load_boundary(const std::string& filename) ;
        void set_non_convex_mode(bool x) { non_convex_mode_ = x ; }
		void normalize_boundary() ;
        void get_bbox(
            real& x_min, real& y_min, real& z_min,
            real& x_max, real& y_max, real& z_max
        ) ;

        bool in_boundary(const vec2& p) { 
            return non_convex_mode_ ? boundary_.contains(p) : boundary_convex_.contains(p) ; 
        }
		bool on_boundary(const vec2& p) { 
			return boundary_.onBounary(p) ; 
        }
        bool& insert_boundary() { return insert_boundary_ ; }

        // ------------------------------------ Delaunay------------------------------ 

        int nb_vertices() const { return all_vertices_.size() ; }
        void clear() ;
        void begin_insert() ;
        Vertex_handle insert(const vec2& p, real w=0) ;
        Vertex_handle nearest(const vec2& p) ;
        double remove(const vec2& p) ;
        void end_insert(bool redraw = true) ;

        std::vector<Vertex_handle>& vertices() { return all_vertices_ ; }

		vec2 vertex(unsigned int i) const {
			return vec2(
			all_vertices_[i]->point().x(),
			all_vertices_[i]->point().y()
			) ;
		}

        // ------------------------------------ Delaunay combinatorics ----------------

        vec2 dual(Face_handle c) { return c->dual ; }
        bool dual_cell_intersects_boundary(Vertex_handle v) {
            return v->dual_intersects_boundary ;
        }
 
        Line<real> get_dual_line(Edge e) {
			vec2 p1 = to_geex(e.first->vertex(ccw(e.second))->point()) ;
			vec2 p2 = to_geex(e.first->vertex( cw(e.second))->point()) ;
			double w1 = e.first->vertex(ccw(e.second))->point().weight() ;
			double w2 = e.first->vertex( cw(e.second))->point().weight() ;
			double t = 0.5*(w2-w1)/distance(p1, p2) ;
			return Line<real>(0.5*(p1+p2)+t*normalize(p1-p2), p2-p1) ;
        }

		vec2 get_mid_point(Point& P1, Point& P2) {
			vec2 p1 = to_geex(P1) ;
			vec2 p2 = to_geex(P2) ;
			double w1 = P1.weight() ;
			double w2 = P2.weight() ;
			double t = 0.5*(w2-w1)/distance(p1, p2) ;
			return 0.5*(p1+p2)+t*normalize(p1-p2) ;
		}

        Polygon2* dual_convex_clip(Vertex_handle v, bool close = true) ;
        int dual_facet_degree(Vertex_handle v) ;

        // -------------- Sampling---------------
		void generate_poisson_disk() ;
		void fill_gaps() ;
		void fill_vertex_gaps() ;
		void fill_face_gaps() ;
		void delete_grid() ;
		void update_grid(double radius) ;
		int  sample_domain_vertices() ;
		int  sample_domain_boundary() ;
		void update_edge_gaps(std::vector<PolygonEdge>& edgegaps) ;
		void sample_edge_gaps(std::vector<PolygonEdge>& edgegaps) ;
		void optimize_boundary_samples() ;
		void get_boundary_primal_edges(std::set<std::pair<int, int>>& primal_edges) ;
		void fill_edge_gaps() ;

		void sample_init() ;
		int  sample_grid_init(int nmaxfail=400) ;
		int  sample_grid_init_random(int nmaxfail=400) ;
		void insert_init_samples() ;
		bool is_hit(int u, int v, vec2 p, double r2, double gridlen) ;
		bool is_hit_fill(int u, int v, vec2 p, double r2, double gridlen) ;
		real weight(const vec2& p) ;
		real weight_fill(const vec2& p, int u, int v, double gridlen) ;
		void mark_hit_cells(int u, int v, double gridlen, vec2 p, double r2, int idx) ;
		bool is_covered(int u, int v, double gridlen, vec2 p, double r2) ;

		// gap computation
		bool has_gap(Face_handle f) ;
		void update_gaps() ;
		int  sample_vertex_gaps(int nmaxfail=400) ;
		int  sample_face_gaps(int nmaxfail=400) ;

		int  cluster_gap_cells() ;
		void update_cell_gaps() ;
		void update_vertex_gaps() ;

		void remark_grid_cells() ;
		//
		double disk_density(const vec2& p, double cur_r) ;
		double update_disk_radius(const vec2& p, double cur_r, double K) ;
		//
		std::vector<GapFaceCluster>& gap_faces() { return gap_faces_ ; }
		std::vector<FaceGap>& facegaps() { return face_gaps_ ; }
		std::vector<VertexGap>& vgaps()  { return vgaps_ ; } 
		GridCell ** grid()   { return grid_ ; }
		int& x_res()         { return x_res_ ; }
		int& y_res()         { return y_res_ ; }
		double& radius_min() { return radius_min_ ; }
		double& radius_max() { return radius_max_ ; }
		double& weight_lp()  { return weight_lp_ ;  }
		double& lfs_weight() { return lfs_weight_ ; }
		const double  gridsize()   { return gridsize_ ; }
		bool&   sample_edges()    { return sample_boundary_ ; }
		bool&   sample_vertices() { return sample_vertices_ ; }
		double& ratio_max_radius() { return ratio_max_radius_ ; }
		std::vector<vec2>& samples() { return samples_ ; }

		void init_sizing_field() ;
		Lipschitz_sizing_field *field() { return field_ ; } 
		void classify_gridcells() ;
		WeightScheme& sample_mode() { return weight_scheme_ ; } 

		void resetDelaunay();
		vec2 grid_res(){ return vec2(x_res_ , y_res_) ;}
		
		//---------------------Triangulating-----------------------------
		void boundary_sample_cluster();
		void maximal_sampling();
		bool is_hit_free(int u, int v, vec2 p, double r2, double gridlen);
		bool diskFreeCheck(vec2 p, int u, int v);
		bool subcell_is_covered(int x_index, int y_index, int iterationLevel, int base_x_index, int base_y_index);
		void proxy_clustering();
		void proxy_clustering_vary();
		void generate_triangle();
		void generate_triangle_vary();
		void show_original_triangle();
		void show_non_delanuary();
		void edge_flip();
		double determinant(double a[], int n);
		void show_non_rt();
		void save_ps_new(const std::string& filename);
	    vec3 ** cluster_colors()   { return cluster_colors_ ; }
		Geex::Map*          map()      { return triangleMesh_ ; }

		//double time_cluster_, time_triangle_, time_check_, time_flip_;

    protected:
        std::vector<Vertex_handle> all_vertices_ ;
        bool non_convex_mode_ ;

        Polygon2 boundary_ ;
        Convex boundary_convex_ ;
		Polygon2 ping_ ;
        Polygon2 pong_ ;
        bool cached_bbox_ ;
        double x_min_, y_min_, x_max_, y_max_ ;
        bool opened_ ;
        bool insert_boundary_ ;

		GapCompute *gap_compute_ ;
		GridCell **grid_ ;
		double   radius_min_ ;
		double   radius_max_ ;
		double   lfs_min_, lfs_max_ ;
		double   weight_lp_ ;
		int      x_res_ ;
		int      y_res_ ;
		double   gridsize_ ;
		int sample_number_;
		std::vector<double> weights_ ;
		std::vector<vec2>   samples_ ;
		std::vector<int>    colors_ ;
		std::vector<vec3>   cell_colors_ ;
		
		std::vector<GapFaceCluster> gap_faces_ ;
		std::vector<VertexGap>      vgaps_ ;
		std::vector<FaceGap>        face_gaps_ ;
		bool     sample_boundary_ ;
		bool     sample_vertices_ ;
		double                  total_density_ ;
		Lipschitz_sizing_field *field_ ;
		double                  Lipschitz_K_ ;
		WeightScheme            weight_scheme_ ;
		double                  lfs_weight_ ;
		double                  ratio_max_radius_ ;

		vec3 **cluster_colors_;
		Geex::Map* triangleMesh_;
    } ;

//______________________________________________________________________________________

#define FOR_EACH_VERTEX_M(TCLASS, INSTANCE, IT)                       \
    for(                                                            \
        TCLASS::Vertex_iterator IT = (INSTANCE)->finite_vertices_begin() ; \
        IT != (INSTANCE)->finite_vertices_end(); IT++                      \
    )

#define FOR_EACH_FACE(TCLASS, INSTANCE, IT)                      \
    for(                                                         \
        TCLASS::Face_iterator IT = (INSTANCE)->faces_begin() ;   \
        IT != (INSTANCE)->faces_end(); IT++                      \
    )
#define FOR_EACH_FINITE_FACE(TCLASS, INSTANCE, IT)                      \
    for(                                                         \
        TCLASS::Finite_faces_iterator IT = (INSTANCE)->finite_faces_begin() ;   \
        IT != (INSTANCE)->finite_faces_end(); IT++                      \
    )
    
#define FOR_EACH_FINITE_EDGE(TCLASS, INSTANCE, IT)                      \
	for(                                                         \
	TCLASS::Finite_edges_iterator IT = (INSTANCE)->finite_edges_begin() ;   \
	IT != (INSTANCE)->finite_edges_end(); IT++                      \
	)
}

#endif
