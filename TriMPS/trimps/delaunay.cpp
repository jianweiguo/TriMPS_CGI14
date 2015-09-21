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
 *  You should have received a copy of the GNU General Public License<
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

#include "delaunay.h"
#include "gap_compute.h"
#include <fstream>
#include <queue>

#ifndef SQR(x)
#define SQR(x) ((x)*(x))
#endif 

namespace Geex {

	Delaunay::Delaunay() {
        non_convex_mode_ = false ;
        cached_bbox_ = false ;
        opened_ = false ;
        insert_boundary_ = false ;
		grid_ = nil ;
		radius_min_ = 0.01 ;
		radius_max_ = 2*radius_min_ ;
		weight_lp_  = 2.0;
		total_density_ = 0 ;
		field_ = nil ;
		Lipschitz_K_ = 0.1 ;
		weight_scheme_ = WT_UNIFORM ;
		lfs_weight_ = 1.0 ;
		gridsize_ = 0 ;
		x_res_ = 0 ;
		y_res_ = 0 ;
		sample_boundary_ = true ;
		sample_vertices_ = true ;
		ratio_max_radius_ = 16 ;
		gap_compute_ = new GapCompute(this) ;

		//
		sample_number_=0;
		cluster_colors_=nil;
		triangleMesh_=new Geex::Map;
    }

	Delaunay::~Delaunay() {
		if(field_!=nil) 
			delete field_ ;
		if(grid_!=nil)  
			delete_grid() ;
		delete gap_compute_ ;
	}
	void Delaunay::resetDelaunay(){
		sample_number_=0;
		cluster_colors_=nil;
		triangleMesh_=new Geex::Map;
		grid_ = nil ;
	}

	static float color1[4] = {
		0.0f, 0.5f, 0.0f, 0.3f
	} ;
	static float color2[4] = {
		0.9f, 0.0f, 0.0f, 0.7f
	} ;
	static float blue[4] = {
		0.1f, 0.7f, 0.8f, 1.0f
	} ;

    void Delaunay::load_boundary(const std::string& filename) {
        cached_bbox_ = false ;
        boundary_.clear() ;
        boundary_convex_.clear() ;
        boundary_.load(filename) ;
        //boundary_.normalize() ;
		normalize_boundary() ;

        if(!non_convex_mode_) {
            for(unsigned int i=0; i<boundary_.size(); i++) {
                boundary_convex_.push_back(boundary_[i].line()) ;
            }
        }

        double x_min, y_min, z_min ;
        double x_max, y_max, z_max ;
        get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
        double dx = x_max - x_min ;
        double dy = y_max - y_min ;
        double diag = sqrt(dx*dx + dy*dy) ;

		//
		init_sizing_field() ;
		update_grid(radius_min_) ;
		std::cout << "domain area " << boundary_.area() << std::endl ;
    }

    void Delaunay::get_bbox(
        real& x_min, real& y_min, real& z_min,
        real& x_max, real& y_max, real& z_max
    ) {
        z_min = 0.0 ; z_max = 1.0 ;
        if(!cached_bbox_) {
            x_min_ = y_min_ =  1e30 ;
            x_max_ = y_max_ = -1e30 ;
            if(boundary_.size() == 0) {
                x_min_ = 0.0 ; y_min_ = 0.0 ;
                x_max_ = 1.0 ; y_max_ = 1.0 ;
            }
            for(unsigned int i=0; i<boundary_.size(); i++) {
                for(unsigned int j=0; j<2; j++) {
                    x_min_ = gx_min(x_min_, boundary_[i].vertex[j].x) ;
                    y_min_ = gx_min(y_min_, boundary_[i].vertex[j].y) ;
                    x_max_ = gx_max(x_max_, boundary_[i].vertex[j].x) ;
                    y_max_ = gx_max(y_max_, boundary_[i].vertex[j].y) ;
                }
            }
            cached_bbox_ = true ;
        }
        x_min = x_min_ ; y_min = y_min_ ;
        x_max = x_max_ ; y_max = y_max_ ;
    }

	void Delaunay::normalize_boundary() {
		double x_min, y_min, z_min ;
        double x_max, y_max, z_max ;
        get_bbox(x_min, y_min, z_min, x_max, y_max, z_max) ;
        double dx = x_max - x_min ;
        double dy = y_max - y_min ;
		double scale = 1.0/max(dx, dy) ;

		for(unsigned int i=0; i<boundary_.size(); ++i) {
			boundary_[i].vertex[0] = scale*(boundary_[i].vertex[0]-vec2(x_min,y_min)) ;
			boundary_[i].vertex[1] = scale*(boundary_[i].vertex[1]-vec2(x_min,y_min)) ;			
		}
		cached_bbox_ = false ;
	}

    inline vec2 random_v() {
        return vec2( 
            Numeric::random_float64(),
            Numeric::random_float64()
        ) ;
    }

	// ------------------------------------ Grid------------------------------------------------
	void Delaunay::delete_grid() {
		if(grid_!=nil) {
			for(int i=0; i<x_res_; ++i) 
				delete [] grid_[i] ;
			delete [] grid_ ;	
		}
		grid_ = nil ;
	}

	void Delaunay::update_grid(double radius) {

		radius_min_ = radius ;
		gridsize_ = radius/sqrt(2.0) ;
		double x = (x_max_-x_min_)/gridsize_ ;
		double y = (y_max_-y_min_)/gridsize_ ;
		int x_res = x-floor(x) > 1e-3 ? ceil(x) : floor(x);
		int y_res = y-floor(y) > 1e-3 ? ceil(y) : floor(y);

		if(x_res_!=x_res || y_res_!=y_res) {
			delete_grid() ;
			x_res_ = x_res ;
			y_res_ = y_res ;
		}

		if(grid_==nil) {
			grid_ = new GridCell*[x_res_] ;
			for(int i=0; i<x_res_; ++i) {
				grid_[i] = new GridCell[y_res_] ;
			}
		}

		for(int i=0; i<x_res_; ++i) {
			for(int j=0; j<y_res_; ++j) {
				grid_[i][j].valid = true ;
				grid_[i][j].type  = CELL_OUTSIDE ;
				grid_[i][j].indices.clear() ;
				//New added
				grid_[i][j].sampled=false;
				grid_[i][j].point.x=-10;
				grid_[i][j].point.y=-10;
				grid_[i][j].clusterIndex.x=-1;
				grid_[i][j].clusterIndex.y=-1;
				grid_[i][j].locked=false;
				grid_[i][j].radius=0.0;

				grid_[i][j].visited=false;
				grid_[i][j].edge_number=0;
				//End added
			}
		}
		// in/out classify
		classify_gridcells() ;
	}

	void Delaunay::classify_gridcells() {
		for(int i=0; i<x_res_; ++i) {
			for(int j=0; j<y_res_; ++j) {
				grid_[i][j].type = CELL_OUTSIDE ;
				grid_[i][j].valid = false ;
			}
		}

		for(unsigned int i=0; i<boundary_.size(); ++i) {
			int u0 = boundary_[i].vertex[0].x/gridsize_ ;
			int v0 = boundary_[i].vertex[0].y/gridsize_ ;
			int u1 = boundary_[i].vertex[1].x/gridsize_ ;
			int v1 = boundary_[i].vertex[1].y/gridsize_ ;
			gx_clamp(u0, 0, x_res_-1) ;
			gx_clamp(v0, 0, y_res_-1) ;
			gx_clamp(u1, 0, x_res_-1) ;
			gx_clamp(v1, 0, y_res_-1) ;
			if(u0>u1) gx_swap(u0, u1) ;
			if(v0>v1) gx_swap(v0, v1) ;
			vec2 n = normalize(boundary_[i].vertex[1] - boundary_[i].vertex[0]) ;
			Line<real> L(boundary_[i].vertex[0], vec2(-n.y, n.x)) ;
			//grid_[u0][v0].type = CELL_BOUNDARY ;
			//grid_[u1][v1].type = CELL_BOUNDARY ;
			for(int ii=u0; ii<=u1; ++ii) {
				for(int jj=v0; jj<=v1; ++jj) {
					bool s1 = L.side(gridsize_*vec2(ii,   jj)) > 0;
					bool s2 = L.side(gridsize_*vec2(ii+1, jj)) > 0;
					bool s3 = L.side(gridsize_*vec2(ii+1, jj+1)) > 0 ;
					bool s4 = L.side(gridsize_*vec2(ii,   jj+1)) > 0 ;
					if(!(s1&&s2&&s3&&s4 || !s1&&!s2&&!s3&&!s4)) {
						grid_[ii][jj].type = CELL_BOUNDARY ;
						grid_[ii][jj].valid = true ;
					}
				}
			}
		}

		// find an inside point
		bool find=false ;
		int  u=0, v=0 ;
		while(!find) {
			vec2 pin = random_v() ;
			u = pin.x/gridsize_ ;
			v = pin.y/gridsize_ ;
			gx_clamp(u, 0, x_res_-1) ;
			gx_clamp(v, 0, y_res_-1) ;
			if(grid_[u][v].type == CELL_OUTSIDE && in_boundary(pin)) {
				find = true ;
				break ;
			}
		}
		
		// floodfill inner cells
		typedef std::pair<int,int> CellID ;
		std::queue<CellID> Q ;
		int offset[][2] = {{-1,0}, {1,0}, {0,-1}, {0,1}} ;
		Q.push(CellID(u, v)) ;
		while(!Q.empty()) {
			CellID top = Q.front() ;
			Q.pop() ;
			u = top.first ; 
			v = top.second ; 
			if(grid_[u][v].type == CELL_OUTSIDE) {
				grid_[u][v].type = CELL_INSIDE ;
				grid_[u][v].valid = true ;
				for(int i=0; i<4; ++i) {
					int uu = u+offset[i][0] ;
					int vv = v+offset[i][1] ;
					if(uu>=0 && uu<x_res_ && vv>=0 && vv<y_res_ && grid_[uu][vv].type == CELL_OUTSIDE)
						Q.push(CellID(uu,vv)) ;
				}
			}
		} ; // end while

		// output max/min lfs value
		lfs_max_=0 ;
		lfs_min_=1e10 ;
		for(int i=0; i<x_res_; ++i) {
			for(int j=0; j<y_res_; ++j) {
				vec2 p = gridsize_ * vec2(i, j) ;
				double lfs = field_->query(to_cgal(p)) ;
				lfs_max_ = gx_max(lfs_max_, lfs) ;
				lfs_min_ = gx_min(lfs_min_, lfs) ;
			}
		}
		//std::cout << "lfs range [" << lfs_min_ <<", " << lfs_max_ << "]" <<std::endl ;
	} 

	void Delaunay::init_sizing_field() {
		if(field_!=nil) delete field_ ;
		std::list<Point> points ;
		for(unsigned int i=0; i<boundary_.size(); ++i) {
			points.push_back(to_cgal(boundary_[i].vertex[0])) ;
			points.push_back(to_cgal(boundary_[i].mid_point())) ;
		}
		field_ = new Lipschitz_sizing_field(points.begin(), points.end(), Lipschitz_K_) ;
	}

	// ----------------------------------- Sampling ---------------------------------------------------
	void Delaunay::generate_poisson_disk() {
		clear() ;
		update_grid(radius_min_) ;

		if(sample_boundary_) {
			sample_domain_boundary() ;
			optimize_boundary_samples();
			update_gaps() ;
		}		

		// init sampling without gap filling
		sample_init() ;

		double density = 0 ;
		for(unsigned int i=0; i<samples_.size(); ++i) {
			density += M_PI*weights_[i]*0.25 ;
		}

		unsigned int order=0;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE ){
					grid_[i][j].vertexOrder=order;
					order++;
					sample_number_++;
				}
			}
		}
	}

	void Delaunay::sample_init() {
		update_gaps() ;

		if(weight_scheme_== WT_RANDOM) {
			sample_grid_init_random() ;
		}else {
			sample_grid_init(400) ;
		}
		insert_init_samples() ;
	}

	int Delaunay::sample_grid_init(int nmaxfail) {
		int    npoint = 0 ;
		int    nmaxthrow = 5*x_res_*y_res_ ;
		int    nminthrow = x_res_*y_res_/16 ;
		int    nthrow = 0 ;
		int    nmiss = 0 ;
		
		std::vector<GridCell> cells ;
		std::vector<double>   weights ;
		std::vector<double>   pdf ;
		double total_w = 0 ;

		for(unsigned int i=0; i<x_res_; ++i) {
			for(unsigned int j=0; j<y_res_; ++j) {
				if(grid_[i][j].type != CELL_OUTSIDE) {
					grid_[i][j].u = i ;
					grid_[i][j].v = j ;
					cells.push_back(grid_[i][j]) ;
					double wij = 1./weight(gridsize_*vec2(i+0.5, j+0.5)) ;
					weights.push_back(wij) ;
					total_w += wij ;
					pdf.push_back(total_w) ;
				}
			}
		}

		for(unsigned int i=0; i<pdf.size(); ++i) {
			pdf[i]/=total_w ; 
		}

		while(nmiss<nmaxfail && nthrow<nmaxthrow || nthrow<nminthrow) {
			int  cidx = random_polygon_by_area(pdf, Numeric::random_float64()) ; 
			int  u = cells[cidx].u ;
			int  v = cells[cidx].v ;
			vec2 p = gridsize_*(vec2(u, v)+random_v()) ; ; 
			double w = weight(p) ;//* SQR(1+0.5*Numeric::random_float64()) ;

			if(w == 0 || grid_[u][v].type == CELL_OUTSIDE 
				|| (grid_[u][v].type == CELL_BOUNDARY && !in_boundary(p))) 
			{ continue ; }

			if(is_hit(u, v, p, w, gridsize_)) {
				mark_hit_cells(u, v, gridsize_, p, w, samples_.size()) ;
				grid_[u][v].sampled=true;
				grid_[u][v].point=p;
				grid_[u][v].radius=sqrt(w);
				samples_.push_back(p) ;
				weights_.push_back(w) ;
				colors_.push_back(samples_.size()-1) ;
				nmiss = 0 ;
			} else {
				nmiss ++ ;
			}

			nthrow ++ ;
		} ;

		return samples_.size() ;
	}

	int Delaunay::sample_grid_init_random(int nmaxfail) {
		int    npoint = 0 ;
		int    nmaxthrow = 5*x_res_*x_res_ ;
		int    nminthrow = x_res_*y_res_/16 ;
		int    nthrow = 0 ;
		int    nmiss = 0 ;

		std::vector<GridCell> cells ;
		std::vector<double>   weights ;
		std::vector<double>   pdf ;
		double total_w = 0 ;

		while(nmiss<nmaxfail && nthrow<nmaxthrow || nthrow<nminthrow) {
			int u = rand()%x_res_ ;
			int v = rand()%x_res_ ;
			vec2 p = gridsize_*(vec2(u, v)+random_v()) ; ; 
			double w = weight(p) ;//* SQR(1+0.5*Numeric::random_float64()) ;

			if(w == 0 || grid_[u][v].type == CELL_OUTSIDE 
				|| (grid_[u][v].type == CELL_BOUNDARY && !in_boundary(p))) 
			{ continue ; }

			if(is_hit(u, v, p, w, gridsize_)) {
				mark_hit_cells(u, v, gridsize_, p, w, samples_.size()) ;
				grid_[u][v].sampled=true;
				grid_[u][v].point=p;
				grid_[u][v].radius=sqrt(w);
				samples_.push_back(p) ;
				weights_.push_back(w) ;
				colors_.push_back(samples_.size()-1) ;
				nmiss = 0 ;
			} else {
				nmiss ++ ;
			}

			nthrow ++ ;
		} ;
		std::cout << "Initial sampling of " << samples_.size() << " points " << std::endl ;

		return samples_.size() ;
	}

	/*
	* sample the sharp features
	*/
	int Delaunay::sample_domain_vertices() {
		begin_insert() ;
		for(unsigned int i=0; i<boundary_.size(); ++i) {
			vec2 p = boundary_[i].vertex[0] ;
			int u = p.x/gridsize_ ;
			int v = p.y/gridsize_ ;
			double w = weight(p) ;
			if(is_hit(u, v, p, w, gridsize_)) {
				Vertex_handle vh = insert(p, w) ; //insert_vertex_periodic(p0, exclude_radius) ;
				vh->index = all_vertices_.size() ;
				vh->locked = true ;
				all_vertices_.push_back(vh) ;
				grid_[u][v].sampled=true;
				grid_[u][v].point=p;
				grid_[u][v].radius=sqrt(w);
				mark_hit_cells(u, v, gridsize_, p, w, samples_.size()) ;
				samples_.push_back(p) ;
				weights_.push_back(w) ;				
			}
		}
		end_insert() ;
		return samples_.size() ;
	}

	bool clip_seg_by_circle(PolygonEdge& e, vec2& c, double radius, std::vector<PolygonEdge>& eout) {
		bool intersect0 = distance(e.vertex[0], c) > radius ;
		bool intersect1 = distance(e.vertex[1], c) > radius ;
		double r2 = SQR(radius) ;

		if(!intersect0 && !intersect1) {
			return false ;
		}
		else if(intersect0 && !intersect1) {
			vec2 pv0 = c+radius*normalize(e.vertex[0]-c) ;
			vec2 v01 = normalize(e.vertex[1] - e.vertex[0]) ;
			Line<double> L(e.vertex[0], vec2(-v01.y, v01.x)) ;
			double h = L.side(c) ;
			double t = sqrt(distance2(e.vertex[0], c)-SQR(h)) - sqrt(r2-SQR(h)) ;
			//vec2 p01 = e.vertex[0] + t*v01 ;
			gx_assert(!Numeric::is_nan(t)) ;
			PolygonEdge newe ; 
			newe.vertex[0] = e.vertex[0] ; // c+radius*normalize(e.vertex[0]-c) ;
			newe.vertex[1] = e.vertex[0] + t*v01 ;
			eout.push_back(newe) ;
		}
		else if(!intersect0 && intersect1) {
			vec2 pv1 = c+radius*normalize(e.vertex[1]-c) ;
			vec2 v01 = normalize(e.vertex[1] - e.vertex[0]) ;
			Line<double> L(e.vertex[0], vec2(-v01.y, v01.x)) ;
			double h = L.side(c) ;
			double t = sqrt(distance2(e.vertex[1], c)-SQR(h)) - sqrt(r2-SQR(h)) ;
			//vec2 p10 = e.vertex[1] - t*v01 ;
			gx_assert(!Numeric::is_nan(t)) ;
			PolygonEdge newe ; 
			newe.vertex[0] = e.vertex[1] - t*v01 ; //c+radius*normalize(e.vertex[1]-c) ;
			newe.vertex[1] = e.vertex[1] ;
			eout.push_back(newe) ;
		}
		else {
			vec2 v01 = normalize(e.vertex[1] - e.vertex[0]) ;
			Line<double> L(e.vertex[0], vec2(-v01.y, v01.x)) ;
			double h = L.side(c) ;
			bool t0 = dot(c-e.vertex[0], v01)>0 ;
			bool t1 = dot(c-e.vertex[1], v01)>0 ;

			if(t0 && t1 || !t0 && !t1) {
				eout.push_back(e) ;
			}
			else{ // two gaps
				if(fabs(h)>radius) {
					eout.push_back(e) ;
				} else {
					vec2   mp = c+h*vec2(v01.y, -v01.x);
					double t = sqrt(r2-SQR(h)) ;
					gx_assert(!Numeric::is_nan(t)) ;

					PolygonEdge e0, e1 ;
					e0.vertex[0] = e.vertex[0] ;
					e0.vertex[1] = mp-t*v01 ;//c+radius*normalize(e.vertex[0]-c) ;
					eout.push_back(e0) ;

					e1.vertex[1] = e.vertex[1] ;
					e1.vertex[0] = mp+t*v01 ; //c+radius*normalize(e.vertex[1]-c) ;
					eout.push_back(e1) ;
				}
			}
		}
		return true ;
	}

	int  Delaunay::sample_domain_boundary()  {
		std::vector<PolygonEdge> edgegaps ;
		std::vector<double> epdf(boundary_.size(), 0) ;
		double total = 0 ;

		if(sample_vertices_) {
			sample_domain_vertices() ;
			//return 1;
		}

		for(int i=0; i<boundary_.size(); ++i) {
			double ew = weight(boundary_[i].mid_point()) ; //field_->query(to_cgal(boundary_[i].mid_point())) ;
			total += boundary_[i].length()/ew ; //SQR(ew) ;
			epdf[i] = total ;
			edgegaps.push_back(boundary_[i]) ;
		}

		for(int i=0; i<epdf.size(); ++i) {
			epdf[i] /= total ;
		}

		begin_insert() ;

		int eidx, u, v; 
		double rnd, w;
		vec2 p0 ;
		Vertex_handle v0 ;

		// step1: random sample edges
		int nmiss = 0 ;
		int npoint = 0 ;
		int nthrow = 0 ;
		int size = epdf.size() ;
		int nmaxthrow = gx_max(5*size, 400) ;
		int nminthrow = gx_max(size/16, 100) ;

		while(nmiss<50 && nthrow<nmaxthrow || nthrow<nminthrow) {
			eidx = random_polygon_by_area(epdf, Numeric::random_float64()) ; 
			rnd = Numeric::random_float64() ;
			p0 = (1-rnd)*edgegaps[eidx].vertex[0] + rnd*edgegaps[eidx].vertex[1] ;
			w = weight(p0) ;
			u = p0.x/gridsize_ ;
			v = p0.y/gridsize_ ;
			gx_clamp(u, 0, x_res_-1) ;
			gx_clamp(v, 0, y_res_-1) ;

			if(is_hit(u, v, p0, w, gridsize_)) {
				Vertex_handle vh = insert(p0, w) ; //insert_vertex_periodic(p0, exclude_radius) ;
				vh->index = all_vertices_.size() ;
				all_vertices_.push_back(vh) ;
			    grid_[u][v].sampled = true ;
				grid_[u][v].point  = p0 ;
				grid_[u][v].radius=sqrt(w);
				this->mark_hit_cells(u, v, gridsize_, p0, w, vh->index) ;
				//samples_.push_back(vec2w(p0, w)) ;
				samples_.push_back(p0) ;
				weights_.push_back(w) ;
				nmiss = 0 ;
				npoint++ ;
			}
			else {
				nmiss ++ ;
				nthrow++ ;
			}
		} ;

		//step2: fill edge gaps
		do {
			update_edge_gaps(edgegaps) ;
			sample_edge_gaps(edgegaps) ;
		} while(edgegaps.size()>0) ;
		
		end_insert() ;
		
		return all_vertices_.size() ;
	}

	void Delaunay::update_edge_gaps(std::vector<PolygonEdge>& edgegaps) {
		// update polygon epdf
		std::vector<PolygonEdge> newedges ;
		double total = 0 ;
		for(int j=0; j<samples_.size(); ++j) {
			newedges.clear() ;				
			for(int i=0; i<edgegaps.size(); ++i) {
				clip_seg_by_circle(edgegaps[i], samples_[j], sqrt(weights_[j]), newedges) ;
			}
			edgegaps.swap(newedges) ;
		}
	}

	void Delaunay::sample_edge_gaps(std::vector<PolygonEdge>& edgegaps) {
		std::vector<double> epdf(edgegaps.size(), 0) ;
		double total = 0 ;

		if(edgegaps.size()==0) return ;

		for(int i=0; i<edgegaps.size(); ++i) {
			double w = weight(edgegaps[i].mid_point()) ;
			total += edgegaps[i].length()/w ;
			epdf[i] = total ;
		}
		for(int i=0; i<epdf.size(); ++i) {
			epdf[i]/=total ;
		}

		do {
			int nmiss = 0 ;
			int npoint = 0 ;
			int nthrow = 0 ;
			int size = epdf.size() ;
			int nmaxthrow = gx_max(5*size, 400) ;
			int nminthrow = gx_max(size/16, 100) ;

			while(nmiss<20 && nthrow<nmaxthrow || nthrow<nminthrow) {
				int    eidx = random_polygon_by_area(epdf, Numeric::random_float64()) ; 
				double rnd = Numeric::random_float64() ;
				vec2   p0 = (1-rnd)*edgegaps[eidx].vertex[0] + rnd*edgegaps[eidx].vertex[1] ;
				int    u = p0.x/gridsize_ ;
				int    v = p0.y/gridsize_ ;
				gx_clamp(u, 0, x_res_-1) ;
				gx_clamp(v, 0, y_res_-1) ;
				double w = weight_fill(p0, u, v, gridsize_) ;

				if(is_hit_fill(u, v, p0, w, gridsize_)) {						
					Vertex_handle vh = insert(p0, w) ; //insert_vertex_periodic(p0, exclude_radius) ;
					vh->index = all_vertices_.size() ;
					all_vertices_.push_back(vh) ;
					grid_[u][v].sampled = true ;
					grid_[u][v].point  = p0 ;
					grid_[u][v].radius=sqrt(w);
					this->mark_hit_cells(u, v, gridsize_, p0, w, vh->index) ;
//						samples_.push_back(vec2w(p0, w)) ;
					samples_.push_back(p0) ;
					weights_.push_back(w) ;
					nmiss = 0 ;
					npoint++ ;
				}
				else {
					nmiss ++ ;
					nthrow++ ;
				}
			} ;

			if(npoint==0) 
				break ;
		} while(edgegaps.size()>0) ;
	}

	void Delaunay::mark_hit_cells(int u, int v, double gridlen, vec2 p, double r2, int idx) {
		int nei = int(sqrt(r2)/gridlen+2);
		for(int i=u-nei; i<=u+nei; ++i) {
			if(i<0 || i>=x_res_) continue ;
			for(int j=v-nei; j<=v+nei; ++j) {
				if(j<0 || j>=y_res_) continue ;
				if(is_covered(i, j, gridlen, p, r2)) {
					grid_[i][j].valid = false ;
					grid_[i][j].indices.push_back(idx) ;
				}
			}
		}
	}

	// test whether a cell is completely covered by a circle
	bool Delaunay::is_covered(int u, int v, double gridlen, vec2 p, double r2) {
		vec2 c(gridlen*(u+0.5), gridlen*(v+0.5)) ;
		double d2 = SQR(fabs(c.x-p.x)+gridlen*0.5) + SQR(fabs(c.y-p.y)+gridlen*0.5) ;
		return d2 < r2 ;
	}

	// optimize edge samples
	void Delaunay::optimize_boundary_samples() {
		int iter=0, max_iter=20 ;

		do {
			// collect primal boundary edges
			std::set<std::pair<int, int>> primal_edges ;
			std::set<Vertex_handle> to_remove ;

			get_boundary_primal_edges(primal_edges) ;

			// mark long edges
			for(std::set<std::pair<int, int>>::iterator it=primal_edges.begin(); it!=primal_edges.end(); ++it) {
				vec2 p0 = to_geex(all_vertices_[it->first]->point()) ;
				vec2 p1 = to_geex(all_vertices_[it->second]->point()) ;
				double r0 = all_vertices_[it->first]->radius ;
				double r1 = all_vertices_[it->second]->radius ;

				//if(distance2(p0, p1)>0.75*SQR(r0+r1)) {
				if(distance2(p0, p1)>0.65*SQR(r0+r1)) {
					if(!all_vertices_[it->first]->locked) {
						to_remove.insert(all_vertices_[it->first]) ;
					}
					if(!all_vertices_[it->second]->locked) {
						to_remove.insert(all_vertices_[it->second]) ;
					}
				}
			}		

			if(to_remove.size()==0) // we cannot return here since we have to fix all the vertices in the end
				break ;

			// remove long edges
			begin_insert() ;
			for(std::set<Vertex_handle>::iterator it=to_remove.begin(); it!=to_remove.end(); ++it) {
				baseclass::remove(*it) ;
			}
			end_insert() ;

			remark_grid_cells() ;

			// resample boundary
			fill_edge_gaps() ;

			iter ++ ;
		} while(iter < max_iter) ;

		// lock optimized boundary vertices
		if(iter < max_iter) {
			for(unsigned int i=0; i<all_vertices_.size(); ++i) {
				all_vertices_[i]->locked = true ;
			}
		}
		std::cout << "Edge sampling optimization Done!" << std::endl ;
	}

	void Delaunay::get_boundary_primal_edges(std::set<std::pair<int, int>>& primal_edges) {
		primal_edges.clear() ;
		for(unsigned int i=0; i<all_vertices_.size(); ++i) {
			if(!all_vertices_[i]->dual_intersects_boundary) 
				continue ;
			Polygon2* P = dual_convex_clip(all_vertices_[i], false) ;
			for(unsigned int j=0; j<P->size(); ++j) {
				std::set<int>& bis = (*P)[j].vertex[0].bisectors ;
				if(bis.size()>0) {
					for(std::set<int>::iterator it=bis.begin(); it!=bis.end(); ++it) {
						gx_assert(*it>=0) ;
						if(*it < i) {
							primal_edges.insert(std::pair<int, int>(*it, i)) ;
						} else {
							primal_edges.insert(std::pair<int, int>(i, *it)) ;
						}
					}
				}
			}
		}
	}

	void Delaunay::fill_edge_gaps() {
		// fill edge gaps
		std::vector<PolygonEdge> edgegaps ;
		for(int i=0; i<boundary_.size(); ++i) {
			edgegaps.push_back(boundary_[i]) ;
		}

		begin_insert() ;
		do {
			update_edge_gaps(edgegaps) ;
			sample_edge_gaps(edgegaps) ;
		} while(edgegaps.size()>0) ;
		end_insert() ;
	}

	void Delaunay::remark_grid_cells() {
		// update samples, update grid
		update_grid(radius_min_) ;

		for(unsigned int i=0; i<samples_.size(); ++i) {
			int u = samples_[i].x/gridsize_ ;
			int v = samples_[i].y/gridsize_ ;
			grid_[u][v].sampled = true ;
			grid_[u][v].point.x = samples_[i].x ;
			grid_[u][v].point.y = samples_[i].y ;
			grid_[u][v].radius=sqrt(weights_[i]);
			mark_hit_cells(u, v, gridsize_, samples_[i], weights_[i], i) ;
		}
	}

	/*
	* fill gaps
	*/
	void Delaunay::fill_gaps() {
		std::cerr << "Filling gaps..." <<std::endl ;
		int nb_gaps = 0 ;
		int nb_sample = 0 ;

		do {
			update_gaps() ;
			if(sample_boundary_) {
				nb_gaps = face_gaps_.size() ;
				nb_sample = sample_face_gaps() ;
			} else {
				nb_gaps = vgaps_.size() ;
				nb_sample = sample_vertex_gaps() ;
			}
		} while (nb_sample > 0) ;	
		//new added
		sample_number_=0;
		unsigned int order=0;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE ){
					grid_[i][j].vertexOrder=order;
					order++;
					sample_number_++;
				}
			}
		}
		std::cerr<<std::endl ;
		std::cerr << "Sample radius: " << radius_min_ <<std::endl ;
		std::cerr << "Number of sampled points: " << sample_number_ <<std::endl ;
	}

	void Delaunay::fill_vertex_gaps() {
		int nb_gaps = 0 ;
		int nb_sample = 0 ;
		do {
			update_vertex_gaps() ;
			nb_gaps = vgaps_.size() ;
			nb_sample = sample_vertex_gaps() ;
		} while (nb_sample > 0) ;		
	}

	void Delaunay::fill_face_gaps() {
		int nb_gaps = 0 ;
		int nb_sample = 0 ;
		do {
			update_gaps() ;
			nb_gaps = vgaps_.size() ;
			nb_sample = sample_face_gaps() ;
		} while (nb_sample > 0) ;		
	}

	bool Delaunay::is_hit(int u, int v, vec2 p, double r2, double gridlen) {
		int nei = int(sqrt(r2)/gridlen+2);
		
		if(!grid_[u][v].valid || grid_[u][v].sampled) { return false ; }

		for(int i=u-nei; i<=u+nei; ++i) {
			if(i<0 || i>=x_res_) continue ;
			for(int j=v-nei; j<=v+nei; ++j) {
				if(j<0 || j>=y_res_) continue ;
				if(grid_[i][j].sampled && !grid_[i][j].valid && grid_[i][j].type!=CELL_OUTSIDE) {
//					gx_assert(grid_[i][j].indices.size() > 0) ;
					std::vector<int>& indices = grid_[i][j].indices ;
					for(int k=0; k<indices.size(); ++k) {
						if(distance2(p, samples_[indices[k]]) < max(r2, weights_[indices[k]]))
							return false ;
					}
				}
			}
		}
		return true ;
	}

	bool Delaunay::is_hit_fill(int u, int v, vec2 p, double r2, double gridlen) {
		if(!grid_[u][v].valid || grid_[u][v].sampled) { return false ; }
		
		int nei = int(sqrt(r2)/gridlen+2);
		for(int i=u-nei; i<=u+nei; ++i) {
			if(i<0 || i>=x_res_) continue ;
			for(int j=v-nei; j<=v+nei; ++j) {
				if(j<0 || j>=y_res_) continue ;
				if(!grid_[i][j].valid && grid_[i][j].type!=CELL_OUTSIDE) {
//					gx_assert(grid_[i][j].indices.size() > 0) ;
					std::vector<int>& indices = grid_[i][j].indices ;
					for(int k=0; k<indices.size(); ++k) {
						// ------ here is different from initial sampling
						if(distance2(p, samples_[indices[k]]) < weights_[indices[k]])
							return false ;
					}
				}
			}
		}
		return true ;
	}

	double Delaunay::disk_density(const vec2& p, double r) {
		int u = p.x/gridsize_ ;
		int v = p.y/gridsize_ ;
		int nei = r/gridsize_ ; // don't plus 1 here
		double density = 0 ;
		int nb=0;
		for(int i=max(0, u-nei); i<=min(x_res_-1, u+nei); ++i) {
			for(int j=max(0, v-nei); j<=min(y_res_-1, v+nei); ++j) {
				//density += 1.0-brightness(gridsize_*vec2(i+0.5, j+0.5), image_) ;
				nb++ ;
			}
		}
		return density/nb*0.25*M_PI*r*r*255 ;
	}

	double Delaunay::update_disk_radius(const vec2& p, double cur_r, double K) {
		double density = disk_density(p, cur_r) ;
		double ratio = sqrt(K/density) ;
		double newr ;
		int    step = 0 ;
		double eps = 0.01 ;
		while(fabs(ratio-1)>eps && step<6) {
			newr = ratio * cur_r ;
			ratio = K/disk_density(p, newr) ;
			step ++ ;
		} ;
		return max(newr, radius_min_) ;
	}
	
	real Delaunay::weight(const vec2& p) {
		radius_max_ = ratio_max_radius_ * radius_min_ ;

		switch(weight_scheme_) {
		case WT_UNIFORM: 
			return SQR(radius_min_) ;
		case WT_RANDOM:
			return SQR(radius_min_*(1+(ratio_max_radius_-1)*Numeric::random_float64())) ;
		case WT_IMAGE: {
		}
		case WT_AIS: {
			double scale = radius_max_-radius_min_ ;
			double K = total_density_ / 10000;
		}
		case WT_FUNCTION:
		{
			// The density function in CCVT paper
			real x = (p.x-0.5)*2 ;
			real y = (p.y-0.5)*2 ;
			real d = exp(-20*(SQR(x)+SQR(y)))+0.2*SQR(sin(M_PI*x)*sin(M_PI*y)) ;
			double r = max(1.0/sqrt(d)*radius_min_, radius_min_) ;
			return SQR(r) ;

			//// linear map
			//double r = radius_min_+pow((1-p.x/(x_max_-x_min_)), weight_lp_)*(radius_max_-radius_min_) ;
			//return r*r ;
		}
		case WT_LFS: {
			//double lfs = lfs_weight_*field_->query(to_cgal(p)) ;
			//double r = max(lfs, radius_min_) ;
			double lfs = field_->query(to_cgal(p)) ;
			double r = gx_max(radius_min_, radius_min_*(ratio_max_radius_-1)*(lfs-lfs_min_)/(lfs_max_-lfs_min_)) ;
			//double r = gx_max(radius_min_, lfs) ;
			return r*r ;
		}
		default:
			return 0 ;
		}
	}

	real Delaunay::weight_fill(const vec2& p, int u, int v, double gridlen) {
		real r2 = weight(p) ;
		real w = r2 ; 
		std::set<int> neighbors ;

		if(w==0.0) return w ;

		// TODO: check the problem
		//Vertex_handle vnst = nearest_power_vertex(to_cgal(p)) ;
		//real t = min(r2, distance2(p, to_geex(vnst->point()))) ;

		int nei = int(sqrt(r2)/gridlen+2);
		for(int i=u-nei; i<=u+nei; ++i) {
			if(i<0 || i>=x_res_) continue ;
			for(int j=v-nei; j<=v+nei; ++j) {
				if(j<0 || j>=y_res_) continue ;
				if(!grid_[i][j].valid && grid_[i][j].type!=CELL_OUTSIDE) {
					gx_assert(grid_[i][j].indices.size() > 0) ;
					std::vector<int>& indices = grid_[i][j].indices ;
					for(int k=0; k<indices.size(); ++k) {
						neighbors.insert(indices[k]) ;
					}
				}
			}
		}

		for(std::set<int>::iterator it=neighbors.begin(); it!=neighbors.end(); ++it) {
			w = min(distance2(p, samples_[*it]), w) ;
		}

		return w ;
	}

	/*
	* when we insert init sample, the number of vertices might not be 0 because the boundray sampling
	*/
	void Delaunay::insert_init_samples() {
		begin_insert() ;
		int n0 = all_vertices_.size() ; 
		for(unsigned int i=n0; i<samples_.size(); ++i) {
			insert(samples_[i], weights_[i]) ;
		}
		end_insert() ;
	}

	void Delaunay::update_gaps() {
		cluster_gap_cells() ;  // cluster gap cells
	
		if(sample_boundary_) { // for viusal debugging
			update_cell_gaps() ;
		} 
		else { // only vertex gaps is used for now.		
			update_vertex_gaps() ;
		}

	}

	bool Delaunay::has_gap(Delaunay::Face_handle f) {
		return distance2(f->dual, to_geex(f->vertex(0)->point())) > f->vertex(0)->point().weight()
			|| distance2(f->dual, to_geex(f->vertex(1)->point())) > f->vertex(1)->point().weight()
			|| distance2(f->dual, to_geex(f->vertex(2)->point())) > f->vertex(2)->point().weight() ;
	}

	int  Delaunay::cluster_gap_cells() {
		int cur_index = 0 ;

		gap_faces_.clear() ;
		FOR_EACH_FACE(Delaunay, this, f) {
			f->visited = false ;
			//f->has_gap = is_infinite(f) ? false : has_gap(f) ; // infinite facet has no gap
		}		

		FOR_EACH_FACE(Delaunay, this, f) {
			if(f->has_gap/*has_gap(f)*/ && !f->visited && !is_infinite(f) && in_boundary(center(f))) {
				GapFaceCluster          gap ;
				std::queue<Face_handle> Q ;

				Q.push(f) ;
				while(!Q.empty()) {
					Face_handle top = Q.front() ;
					Q.pop() ;
					if(!top->visited) {
						top->visited = true ;
						top->gap_index = cur_index ;
						gap.push_back(top) ;
						gap.is_boundary_ |= top->dual_outside ; // indicate whether a gap connect to boundary or not

						for(int i=0; i<3; ++i) {
							Face_handle fn = top->neighbor(i) ;
							if(fn->has_gap && !fn->visited && !is_infinite(fn) && in_boundary(center(fn))) {
								Vertex_handle v0 = top->vertex(top->cw(i)) ;
								Vertex_handle v1 = top->vertex(top->ccw(i)) ;
								vec2 p0 = to_geex(v0->point()) ;
								vec2 p1 = to_geex(v1->point()) ;
								if(distance(p0, p1)>sqrt(v0->point().weight())+sqrt(v1->point().weight())) {
									Q.push(fn) ;
								}
								else if(same_side(p0, p1, top->dual, fn->dual)) {
									Q.push(fn) ;
								}
							}
							//if(has_gap(fn) && !top->neighbor(i)->visited  && !is_infinite(fn)) {
							//	Q.push(fn) ;
							//}
						}
					}
				} ;

				gap_faces_.push_back(gap) ;
				cur_index ++ ;
			}
		}
		//std::cout << "number of gaps: " << gap_faces_.size() << std::endl ; 
		return gap_faces_.size() ;
	}

	void Delaunay::update_cell_gaps() {
		face_gaps_.clear() ;
		for(unsigned int i=0; i<gap_faces_.size(); ++i) {
			GapFaceCluster& cluster = gap_faces_[i] ;
			// compute gaps connected to boundary
			if(cluster.is_boundary_) {
				for(unsigned int j=0; j<cluster.size(); ++j) {
					Face_handle f = cluster[j] ;
					for(int k=0; k<3; ++k) {
						Vertex_handle v = f->vertex(k) ;
						if(v->dual_intersects_boundary) {
						} 
						else {
						}
					}
				}
			}
			// compute inner gaps
			else {
				for(unsigned int j=0; j<cluster.size(); ++j) {
					gap_compute_->compute_face_gap(cluster[j], face_gaps_) ;
				}
			}
		}
	}

	void Delaunay::update_vertex_gaps() {
		vgaps_.clear() ;
		FOR_EACH_VERTEX_M(Delaunay, this, it) {
			gap_compute_->compute_vertex_gaps(it, vgaps_) ;
		}
		std::cout << "number of vertex gaps " << vgaps_.size() << std::endl ;
	}

	int  Delaunay::sample_face_gaps(int nmaxfail) {
		unsigned int ngaps = face_gaps_.size() ;
		if(ngaps==0) return 0 ;
		
		std::vector<double> area(ngaps, 0) ;
		std::vector<double> cpdf(ngaps, 0) ;
		double total_area = 0 ;

		for(unsigned int i=0; i<ngaps; ++i) {
			area[i] = face_gaps_[i].area()/weight(face_gaps_[i].center()) ;
			total_area += area[i] ;
			cpdf[i] = total_area ;
		}
		
		for(unsigned int i=0; i<ngaps; ++i) {
			cpdf[i] /= total_area ;
		}

		int    nmiss = 0 ;
		int    npoint = 0 ;
		int    nthrow = 0 ;
		int size = face_gaps_.size() ;
		int nmaxthrow = max(3*size, 400) ;
		int nminthrow = max(size, 100) ;
		int nv0 = samples_.size() ;

		while(nmiss<nmaxfail && nthrow<nmaxthrow || nthrow<nminthrow) {
			int idx = random_polygon_by_area(cpdf, Numeric::random_float64()) ; //Random::UniformRandom()) ;
			// vec2 p = face_gaps_[idx].random_point() ;
			vec2 p ;
			if(!face_gaps_[idx].random_point(p))
				continue ;
			int  u = p.x/gridsize_ ;
			int  v = p.y/gridsize_ ;
			gx_clamp(u, 0, x_res_-1) ;
			gx_clamp(v, 0, y_res_-1) ;

			if(!in_boundary(p)) 
				continue ;

			/*
			** in the gap sampling stage, we have change the heuristic of hit test.
			** The hit is success if the new point is not covered by the old ones.
			** the radius of the new point is descided by the distance to the existing
			** points. 
			*/ 
			double w = weight_fill(p, u, v, gridsize_) ;
			if(w==0.0) continue ; // 
			if(is_hit_fill(u, v, p, w, gridsize_)) {
				mark_hit_cells(u, v, gridsize_, p, w, samples_.size()) ;
				grid_[u][v].sampled = true ;
				grid_[u][v].point  = p ;
				grid_[u][v].radius=sqrt(w);
				samples_.push_back(p) ;
				weights_.push_back(w) ;
				colors_.push_back(samples_.size()-1) ;
				//new_points.push_back(vec2w(p,w)) ;
				nmiss = 0 ;
			} else {
				nmiss ++ ;
			}

			nthrow ++ ;
		} ;

		// to do: accelerate insertion by indicate the facet
		// this can be done in constant time since we know the location of new sample	
		begin_insert() ;
		for(unsigned int i=nv0; i<samples_.size(); ++i) {
			Vertex_handle vh ;
			vh = insert(samples_[i], weights_[i]) ;
		}
		end_insert() ;

		return samples_.size()-nv0 ;
	}

	int  Delaunay::sample_vertex_gaps(int nmaxfail) {
		unsigned int ngaps = vgaps_.size() ;
		std::vector<double> area ;
		std::vector<double> cpdf ;
		double total_area = 0 ;
		area.assign(ngaps, 0.0) ;
		cpdf.assign(ngaps, 0.0) ;

		for(unsigned int i=0; i<ngaps; ++i) {
//			gx_assert(in_boundary(vgaps_[i].center())) ;
			area[i] = vgaps_[i].area()/weight(vgaps_[i].center()) ;

			if(area[i]<0) { // ignore this gap
				std::cout << "gap " << i << ": negative area: " << area[i] << std::endl ;
				area[i] = 0 ;
			}
			if(!in_boundary(vgaps_[i].center())) {// ignore this gap
				std::cerr << "gap " << i << ": center outside: " << vgaps_[i].center() << std::endl ;
				area[i] = 0 ;
			}
			total_area += area[i] ;
			cpdf[i] = total_area ;

		}
		
		// TODO: fix the bug
		if(total_area <= 0) return 0 ;

		for(unsigned int i=0; i<ngaps; ++i) {
			cpdf[i] /= total_area ;
		}

//		double gridlen = 1.0/gridres_ ;
		int    nmiss = 0 ;
		int    npoint = 0 ;
		int    nthrow = 0 ;
		int size = vgaps_.size() ;
		int nmaxthrow = 3*size ;
		int nminthrow = max(size, size/16) ;
//		std::vector<vec2w> new_points ;
		int nv0 = samples_.size() ;

		while(nmiss<nmaxfail && nthrow<nmaxthrow || nthrow<nminthrow) {
			int idx = random_polygon_by_area(cpdf, Numeric::random_float64()) ; //Random::UniformRandom()) ;
			vec2 p = vgaps_[idx].random_point() ;
			int  u = p.x/gridsize_ ;
			int  v = p.y/gridsize_ ;
			gx_clamp(u, 0, x_res_-1) ;
			gx_clamp(v, 0, y_res_-1) ;

			if(!in_boundary(p)) 
				continue ;

			/*
			** in the gap sampling stage, we have change the heuristic of hit test.
			** The hit is success if the new point is not covered by the old ones.
			** the radius of the new point is descided by the distance to the existing
			** points. 
			*/ 
			double w = weight_fill(p, u, v, gridsize_) ;
			if(w==0.0) continue ; // 
			if(is_hit_fill(u, v, p, w, gridsize_)) {
				mark_hit_cells(u, v, gridsize_, p, w, samples_.size()) ;
				grid_[u][v].sampled = true ;
				grid_[u][v].point  = p ;
				grid_[u][v].radius=sqrt(w);
				samples_.push_back(p) ;
				weights_.push_back(w) ;
				colors_.push_back(samples_.size()-1) ;
				//new_points.push_back(vec2w(p,w)) ;
				nmiss = 0 ;
			} else {
				nmiss ++ ;
			}

			nthrow ++ ;
		} ;

		// to do: accelerate insertion by indicate the facet
		// this can be done in constant time since we know the location of new sample	
		begin_insert() ;
		for(unsigned int i=nv0; i<samples_.size(); ++i) {
			Vertex_handle vh ;
			vh = insert(samples_[i], weights_[i]) ;
		}
		end_insert() ;

		return samples_.size()-nv0 ;
	} 

    // ------------------------------------ Delaunay-------------------------- 
    void Delaunay::clear() {
        baseclass::clear() ;
        all_vertices_.clear() ;
		samples_.clear() ;
		weights_.clear() ;
		gap_faces_.clear() ;
    }

    Delaunay::Vertex_handle Delaunay::insert(const vec2& p, real w) {
        gx_assert(opened_) ;
        if(Numeric::is_nan(p.x) || Numeric::is_nan(p.y)) {
            std::cerr << "Nan !" << std::endl ;
            return 0 ;
        }
		Point wp(to_cgal(p), w) ;
        Vertex_handle result = baseclass::insert(wp) ;
		result->radius = sqrt(w) ;
        result->locked = false ;
        return result ;
    }

    Delaunay::Vertex_handle Delaunay::nearest(const vec2& p) {
        Delaunay::Face_handle f = locate(to_cgal(p)) ;
        Delaunay::Vertex_handle result = f->vertex(0) ;
        double dist = (to_geex(result->point()) - p).length2() ;
        for(unsigned int i=1; i<3; i++) {
            double cur_dist = (to_geex(f->vertex(i)->point()) - p).length2() ;
            if(cur_dist < dist) {
                dist = cur_dist ;
                result = f->vertex(i) ;
            }
        }
        return result ;
    }

    double Delaunay::remove(const vec2& p) {
        if(number_of_vertices() <= 3) { return -1 ; }
        gx_assert(opened_) ;
        
		Face_handle f = locate(to_cgal(p)) ;
		double w = 0 ;
        double min_d = 1e30 ;
        Vertex_handle v = 0 ;

        for(unsigned int i=0; i<3; i++) {
            if(!is_infinite(f->vertex(i))) {
                double cur_d = (to_geex(f->vertex(i)->point()) - p).length2() ;
                if(cur_d < min_d) {
                    min_d = cur_d ;
                    v = f->vertex(i) ;
                }
            }
        }
		w = v->point().weight() ;
        baseclass::remove(v) ;
		return w ;
    }

    void Delaunay::begin_insert() { opened_ = true ; }

    void Delaunay::end_insert(bool redraw) {
        all_vertices_.clear() ;
        for(All_vertices_iterator it = all_vertices_begin() ; it != all_vertices_end() ; it++) {
            it->dual_intersects_boundary = false ;
            it->index = -1 ;
			it->radius = sqrt(it->point().weight()) ;
        }

        for(All_faces_iterator it = all_faces_begin() ; it != all_faces_end() ; it++) {
            if(is_infinite(it)) {
                it->infinite = true ;
                it->dual = vec2(0, 0) ;
                it->dual_outside = true ;
				it->has_gap = false ;
            } else {
                it->infinite = false ;
                it->dual = to_geex(baseclass::dual(it)) ;
                it->dual_outside = !in_boundary(it->dual) ;
				it->has_gap = has_gap(it) ;
            }

            if(it->dual_outside) {
                it->vertex(0)->dual_intersects_boundary = true ;
                it->vertex(1)->dual_intersects_boundary = true ;
                it->vertex(2)->dual_intersects_boundary = true ;
            }
        }        
        
        all_vertices_.clear() ;
        int cur_index = 0 ;
        for(Vertex_iterator it = finite_vertices_begin(); it != finite_vertices_end() ; it++) {
            it->index = cur_index ;
            all_vertices_.push_back(it) ;
            cur_index++ ;
        }
        opened_ = false ;

		// FIX: the problem caused by CGAL insert/remove
		samples_.resize(all_vertices_.size()) ;
		weights_.resize(all_vertices_.size()) ;
		for(unsigned int i=0; i<all_vertices_.size(); ++i) {
			samples_[i] = to_geex(all_vertices_[i]->point()) ;
			weights_[i] = all_vertices_[i]->point().weight() ;
		}
    }

    Polygon2* Delaunay::dual_convex_clip(Vertex_handle v, bool close) {
        Polygon2* from = &boundary_ ;
        Polygon2* to = &pong_ ;

        std::vector<Edge> edges ;
        Edge_circulator it = incident_edges(v) ;
        do {
            edges.push_back(*it) ;
            it++ ;
        } while(it != incident_edges(v)) ;

        // No need to have special cases for infinite vertices: they
        // do not yield any clipping plane !
        for(unsigned int i=0; i<edges.size(); i++) {
            Edge e = edges[i] ;
            if(is_infinite(e)) { continue ; }
            Geex::Line<real> L = get_dual_line(e) ;
            int E = e.first->vertex(ccw(e.second))->index ;
            gx_assert(E != v->index) ;
            from->convex_clip(*to, L, E, close) ;
            if(from == &boundary_) {
                from = &pong_ ;
                to = &ping_ ;
            } else {
                gx_swap(from, to) ;
            }
        }
        return from ;
    }

    int Delaunay::dual_facet_degree(Vertex_handle v) {
        int result = 0 ;
        if(dual_cell_intersects_boundary(v)) { 
            Polygon2* P = dual_convex_clip(v) ;
            result = P->size() ;
        } else {
            Face_circulator it = incident_faces(v) ;
            do {
                result++ ;
                it++ ;
            } while(it != incident_faces(v)) ;
        } 
        return result ;
    }

	//------------------------Triangulating the MPS/PS sampling------------------

	//here we re-implement the MPS for the uniform sampling, if you want the adaptive sampling, please use the method above
    void Delaunay::maximal_sampling(){
		clear() ;
		update_grid(radius_min_) ;
		if(sample_boundary_) {
			sample_domain_boundary() ;
			optimize_boundary_samples();
		}

		int B=10;
		int fixArrayLength=B*x_res_*y_res_;
		int_two *activeCells=(int_two *)malloc(sizeof(int_two)*fixArrayLength);
		for(int i=0;i<fixArrayLength;i++){
			activeCells[i].x=-1;
			activeCells[i].y=-1;
		}
		int_two removeFlag;
		removeFlag.x=-2;
		removeFlag.y=-2;

		int iterationLevel=0;
		int activeCellSize=0;
		for(int i=0;i<x_res_; i++){
			for(int j=0;j<y_res_;j++){
				if(!grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE){
					activeCells[activeCellSize].x=i;
					activeCells[activeCellSize].y=j;
					activeCellSize++;
				}
			}
		}
		//we begin iteration untill there are no active cells or the level bigger than the bits of precision (here is 23)
		double A=2;
		srand((unsigned int)time(0));
		while(activeCellSize>0 && iterationLevel<=23){
			int dartThrowTimes=(int)(A*activeCellSize);
			int pow=1;
			for(int i=1;i<=iterationLevel;i++){
				pow=pow*2;
			}
			double gridsize_cur=gridsize_/pow;
			for(int j=0;j<dartThrowTimes;j++){
				if(activeCellSize>0){
					int randomCell=(rand()*((int)activeCellSize/RAND_MAX+1)) % activeCellSize;
					//We locate the base cell from any subcell (x,y)
					int x_indices=activeCells[randomCell].x, y_indices=activeCells[randomCell].y;
					int base_x_indices=(int)(x_indices/pow), base_y_indices=(int)(y_indices/pow);
					if(grid_[base_x_indices][base_y_indices].sampled){
						activeCells[randomCell]=activeCells[activeCellSize-1];
						activeCells[activeCellSize-1].x=removeFlag.x;
						activeCells[activeCellSize-1].y=removeFlag.y;
						activeCellSize--;
					}
					else{
						vec2 randomPoint = gridsize_cur*(vec2(x_indices, y_indices)+random_v()) ;
						if(diskFreeCheck(randomPoint,base_x_indices,base_y_indices) && (in_boundary(randomPoint) || on_boundary(randomPoint))){
							grid_[base_x_indices][base_y_indices].point=randomPoint;
							grid_[base_x_indices][base_y_indices].sampled=true;
							grid_[base_x_indices][base_y_indices].radius=radius_min_;
							samples_.push_back(randomPoint) ;
							weights_.push_back(SQR(radius_min_)) ;
						}
					}									
				}
			}
 			int newSize=activeCellSize;
			int_two *activeCellsCopy=(int_two *)malloc(sizeof(int_two)*fixArrayLength);
			for(int i=0;i<fixArrayLength;i++){
				activeCellsCopy[i].x=-1;
				activeCellsCopy[i].y=-1;
			}
			int copySize=0;
			for(int k=0;k<newSize;k++){
				int x_indicesCur=activeCells[k].x, y_indicesCur=activeCells[k].y;
				int base_x_indicesCur=(int)(x_indicesCur/pow), base_y_indicesCur=(int)(y_indicesCur/pow);
				if(!grid_[base_x_indicesCur][base_y_indicesCur].sampled && grid_[base_x_indicesCur][base_y_indicesCur].type != CELL_OUTSIDE){
					for(int i=0;i<=1;i++){
						for(int j=0;j<=1;j++){
							vec2 currentIndex(2*activeCells[k].x+i,2*activeCells[k].y+j);
							if(!subcell_is_covered(currentIndex.x, currentIndex.y, iterationLevel+1, base_x_indicesCur, base_y_indicesCur)){
								activeCellsCopy[copySize].x=currentIndex.x;
								activeCellsCopy[copySize].y=currentIndex.y;
								copySize++;
							}												
						}
					}
				}
			}
			activeCellSize=0;
			for(int i=0;i<fixArrayLength;i++){
				if(activeCellsCopy[i].x!=-1 && activeCellsCopy[i].y!=-1){
					activeCells[i]=activeCellsCopy[i];
					activeCellSize++;
				}
				else{
					break;
				}
			}
			delete []activeCellsCopy;
			iterationLevel++;
		}
		unsigned int order=0;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE ){
					grid_[i][j].vertexOrder=order;
					order++;
					sample_number_++;
				}
			}
		}
		std::cerr << "Sample radius: " << radius_min_ <<std::endl ;
		std::cerr << "Number of sampled points: " << sample_number_ <<std::endl ;
	}

	bool Delaunay::diskFreeCheck(vec2 p, int u, int v){

		for(int i=u-2; i<=u+2; ++i) {
			if(i<0 || i>=x_res_) continue ;
			for(int j=v-2; j<=v+2; ++j) {
				if(j<0 || j>=y_res_) continue ;
				if(grid_[i][j].sampled && grid_[i][j].type!=CELL_OUTSIDE && (i!=u || j!=v)) {
					double distance = SQR(p.x-grid_[i][j].point.x)+SQR(p.y-grid_[i][j].point.y);
					if(distance <= SQR(radius_min_)){
						return false;
					}
				}
			}
		}
		return true ;
	}

	bool Delaunay::is_hit_free(int u, int v, vec2 p, double r2, double gridlen) {
		if(grid_[u][v].sampled) { return false ; }

		for(int i=u-2; i<=u+2; ++i) {
			if(i<0 || i>=x_res_) continue ;
			for(int j=v-2; j<=v+2; ++j) {
				if(j<0 || j>=y_res_) continue ;
				if(grid_[i][j].sampled && grid_[i][j].type!=CELL_OUTSIDE && (i!=u || j!=v)) {
					double distance = SQR(p.x-grid_[i][j].point.x)+SQR(p.y-grid_[i][j].point.y);
					if(distance <= SQR(radius_min_)){
						return false;
					}
				}
			}
		}
		return true ;
	}

	bool Delaunay::subcell_is_covered(int x_index, int y_index, int iterationLevel, int base_x_index, int base_y_index){
		int pow=1;
		for(int i=1;i<=iterationLevel;i++){
			pow=pow*2;
		}
		double new_gridlen=gridsize_/pow;
		double r2=SQR(radius_min_);
		vec2 left_low,cellCenter;
		left_low.x=x_index*new_gridlen;
		left_low.y=y_index*new_gridlen;
		cellCenter.x=left_low.x+(new_gridlen/2);
		cellCenter.y=left_low.y+(new_gridlen/2);
		//check if the cell is coveraged by any nearby cell sample
		for(int i=base_x_index-2;i<=base_x_index+2;i++){
			for(int j=base_y_index-2;j<=base_y_index+2;j++){
				if(i>=0 && j>=0 && i<=x_res_-1 && j<=y_res_-1){
					if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE){
						double dx=fabs(cellCenter.x-grid_[i][j].point.x)+(new_gridlen*0.5);
						double dy=fabs(cellCenter.y-grid_[i][j].point.y)+(new_gridlen*0.5);
						double distance = SQR(dx)+SQR(dy);
						if(distance <= r2){
							return true;
						}
					}
				}
			}
		}
		if(!in_boundary(left_low) && !on_boundary(left_low)){
			return true;
		}
		return false;

	}

	//cluster grid cells for uniform or adaptive
	void Delaunay::boundary_sample_cluster(){

		if(cluster_colors_==nil) {
			cluster_colors_ = new vec3*[x_res_] ;
			for(int i=0; i<x_res_; ++i) {
				cluster_colors_[i] = new vec3[y_res_] ;
			}

			for(int i=0; i<x_res_; ++i) {
				for(int j=0; j<y_res_; ++j) {
					cluster_colors_[i][j]=vec3(-1.0,-1.0,-1.0);
				}
			}
		}
		std::priority_queue<cellPair, std::vector<cellPair>, std::greater<cellPair> > PQ ;
		vec3 newcolor;
		int_two id, cen_id;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE){
					newcolor.x=Numeric::random_float64();
					newcolor.y=Numeric::random_float64();
					newcolor.z=Numeric::random_float64();
					cluster_colors_[i][j]=newcolor;

					id.x=i; id.y=j;
					cen_id.x=i; cen_id.y=j;
					PQ.push(cellPair(grid_[i][j].point, grid_[i][j].point, grid_[i][j].radius, id, cen_id));
				}
				if(!grid_[i][j].sampled || grid_[i][j].type == CELL_OUTSIDE){
					vec2 c(gridsize_*(i+0.5), gridsize_*(j+0.5)) ;
					grid_[i][j].point=c;

					newcolor.x=newcolor.y=newcolor.z=0.0;
					cluster_colors_[i][j]=newcolor;
				}
			}
		}
		while(!PQ.empty()){
			cellPair cell=PQ.top();
			PQ.pop();
			id.x=cell.id_.x; id.y=cell.id_.y;
			cen_id.x=cell.center_id_.x; cen_id.y=cell.center_id_.y;
			if(grid_[id.x][id.y].type == CELL_BOUNDARY && grid_[id.x][id.y].clusterIndex.x == -1 && grid_[id.x][id.y].clusterIndex.y == -1){
				grid_[id.x][id.y].clusterIndex.x=cen_id.x;
				grid_[id.x][id.y].clusterIndex.y=cen_id.y;
				grid_[id.x][id.y].locked=true;
	
				int_two cur_cell_id;
				if(id.x-1>=0){
					cur_cell_id.x=id.x-1; cur_cell_id.y=id.y;
					if(grid_[cur_cell_id.x][cur_cell_id.y].type == CELL_BOUNDARY && grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 && grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
				if(id.y-1>=0){
					cur_cell_id.x=id.x; cur_cell_id.y=id.y-1;
					if(grid_[cur_cell_id.x][cur_cell_id.y].type == CELL_BOUNDARY && grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 && grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
				if(id.x+1<=x_res_-1){
					cur_cell_id.x=id.x+1; cur_cell_id.y=id.y;
					if(grid_[cur_cell_id.x][cur_cell_id.y].type == CELL_BOUNDARY && grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 && grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
				if(id.y+1<=y_res_-1){
					cur_cell_id.x=id.x; cur_cell_id.y=id.y+1;
					if(grid_[cur_cell_id.x][cur_cell_id.y].type == CELL_BOUNDARY && grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 && grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
			}
		}
	}

	void Delaunay::proxy_clustering(){
		if(cluster_colors_==nil) {
			cluster_colors_ = new vec3*[x_res_] ;
			for(int i=0; i<x_res_; ++i) {
				cluster_colors_[i] = new vec3[y_res_] ;
			}

			for(int i=0; i<x_res_; ++i) {
				for(int j=0; j<y_res_; ++j) {
					cluster_colors_[i][j]=vec3(-1.0,-1.0,-1.0);
				}
			}
		}
		

		std::priority_queue<cellPair, std::vector<cellPair>, std::greater<cellPair> > PQ ;
		//initial seeds
		vec3 newcolor;
		int_two id, cen_id;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE){
					if(!grid_[i][j].locked){
						newcolor.x=Numeric::random_float64();
						newcolor.y=Numeric::random_float64();
						newcolor.z=Numeric::random_float64();
						cluster_colors_[i][j]=newcolor;
					}
					id.x=i; id.y=j;
					cen_id.x=i; cen_id.y=j;
					PQ.push(cellPair(grid_[i][j].point, grid_[i][j].point, grid_[i][j].radius, id, cen_id));
				}
				if(!grid_[i][j].sampled || grid_[i][j].type == CELL_OUTSIDE){
					if(!grid_[i][j].locked){
						vec2 c(gridsize_*(i+0.5), gridsize_*(j+0.5)) ;
						grid_[i][j].point=c;

						newcolor.x=newcolor.y=newcolor.z=0.0;
						cluster_colors_[i][j]=newcolor;
					}
				}
			}
		}
		//cluster
		while(!PQ.empty()){
			cellPair cell=PQ.top();
			PQ.pop();
			id.x=cell.id_.x; id.y=cell.id_.y;
			cen_id.x=cell.center_id_.x; cen_id.y=cell.center_id_.y;
			if(grid_[id.x][id.y].clusterIndex.x == -1 || grid_[id.x][id.y].clusterIndex.y == -1){
				grid_[id.x][id.y].clusterIndex.x=cen_id.x;
				grid_[id.x][id.y].clusterIndex.y=cen_id.y;
	
				int_two cur_cell_id;
				if(id.x-1>=0){
					cur_cell_id.x=id.x-1; cur_cell_id.y=id.y;
					//if(!grid_cell_->base_cell_[cur_cell_id].outdomain && grid_cell_->base_cell_[cur_cell_id].clusterIndex==-1){
					if((grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 || grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1) && grid_[cur_cell_id.x][cur_cell_id.y].type != CELL_OUTSIDE){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
				if(id.y-1>=0){
					cur_cell_id.x=id.x; cur_cell_id.y=id.y-1;
					if((grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 || grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1) && grid_[cur_cell_id.x][cur_cell_id.y].type != CELL_OUTSIDE){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
				if(id.x+1<=x_res_-1){
					cur_cell_id.x=id.x+1; cur_cell_id.y=id.y;
					if((grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 || grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1) && grid_[cur_cell_id.x][cur_cell_id.y].type != CELL_OUTSIDE){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
				if(id.y+1<=y_res_-1){
					cur_cell_id.x=id.x; cur_cell_id.y=id.y+1;
					if((grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 || grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1) && grid_[cur_cell_id.x][cur_cell_id.y].type != CELL_OUTSIDE){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
			}
		}
	}

	void Delaunay::proxy_clustering_vary(){
		if(cluster_colors_==nil) {
			cluster_colors_ = new vec3*[x_res_] ;
			for(int i=0; i<x_res_; ++i) {
				cluster_colors_[i] = new vec3[y_res_] ;
			}

			for(int i=0; i<x_res_; ++i) {
				for(int j=0; j<y_res_; ++j) {
					cluster_colors_[i][j]=vec3(-1.0,-1.0,-1.0);
				}
			}
		}
		vec3 newcolor;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE){
					if(!grid_[i][j].locked){
						newcolor.x=Numeric::random_float64();
						newcolor.y=Numeric::random_float64();
						newcolor.z=Numeric::random_float64();
						cluster_colors_[i][j]=newcolor;
					}
				}
				if(!grid_[i][j].sampled || grid_[i][j].type == CELL_OUTSIDE){
					if(!grid_[i][j].locked){
						vec2 c(gridsize_*(i+0.5), gridsize_*(j+0.5)) ;
						grid_[i][j].point=c;
						newcolor.x=newcolor.y=newcolor.z=0.0;
						cluster_colors_[i][j]=newcolor;
					}
				}
			}
		}
		vec2 sample1,sample2;
		double weight1,weight2;
		int u1,v1,u2,v2;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(!grid_[i][j].locked && !grid_[i][j].valid && grid_[i][j].type!=CELL_OUTSIDE) {
//					gx_assert(grid_[i][j].indices.size() > 0) ;
					std::vector<int>& indices = grid_[i][j].indices ;
					if(indices.size()==1){
						sample1=samples_[indices[0]];
						u1=sample1.x/gridsize_;
					    v1=sample1.y/gridsize_;
						grid_[i][j].clusterIndex.x=u1;
						grid_[i][j].clusterIndex.y=v1;
					}
					else if(indices.size()==2){
						sample1=samples_[indices[0]];
						sample2=samples_[indices[1]];
						weight1=weights_[indices[0]];
						weight2=weights_[indices[1]];
						u1=sample1.x/gridsize_; v1=sample1.y/gridsize_;
						u2=sample2.x/gridsize_; v2=sample2.y/gridsize_;
						vec2 c(gridsize_*(i+0.5), gridsize_*(j+0.5)) ;
						double dis1=distance2(grid_[i][j].point, grid_[u1][v1].point)- weight1;
						double dis2=distance2(grid_[i][j].point, grid_[u2][v2].point)- weight2;
						if(dis1<dis2){
							grid_[i][j].clusterIndex.x=u1;
						    grid_[i][j].clusterIndex.y=v1;
						}
						else{
							grid_[i][j].clusterIndex.x=u2;
						    grid_[i][j].clusterIndex.y=v2;
						}
					}
				}
			}
		}

		std::priority_queue<cellPair, std::vector<cellPair>, std::greater<cellPair> > PQ ;
		//initial seeds
		int_two id, cen_id, nei_cell;
		int_two nei[4];
	    nei[0].x=-1; nei[0].y=0;
	    nei[1].x=0;  nei[1].y=1;
	    nei[2].x=1;  nei[2].y=0;
	    nei[3].x=0;  nei[3].y=-1;
		std::vector<int_two> reset;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].clusterIndex.x!=-1 &&  grid_[i][j].clusterIndex.y!=-1 && grid_[i][j].type != CELL_OUTSIDE){
					bool find=false;
					for(int k=0;k<4;k++){
						nei_cell.x=i+nei[k].x;
						nei_cell.y=j+nei[k].y;
						if(nei_cell.x>=0 && nei_cell.x<x_res_ && nei_cell.y>=0 && nei_cell.y<y_res_){
							if(grid_[nei_cell.x][nei_cell.y].clusterIndex.x!=-1 &&  grid_[nei_cell.x][nei_cell.y].clusterIndex.y!=-1){
								continue;
							}
							find=true;
							break;
						}
					}
					if(find){
						id.x=i; id.y=j;
					    cen_id.x=grid_[i][j].clusterIndex.x; cen_id.y=grid_[i][j].clusterIndex.y;
					    PQ.push(cellPair(grid_[i][j].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, id, cen_id));
						reset.push_back(id);
					}
				}
			}
		}
		for(int i=0;i<reset.size();i++){
			grid_[reset[i].x][reset[i].y].clusterIndex.x=-1;
			grid_[reset[i].x][reset[i].y].clusterIndex.y=-1;
		}
		//cluster
		while(!PQ.empty()){
			cellPair cell=PQ.top();
			PQ.pop();
			id.x=cell.id_.x; id.y=cell.id_.y;
			cen_id.x=cell.center_id_.x; cen_id.y=cell.center_id_.y;
			if(grid_[id.x][id.y].clusterIndex.x == -1 || grid_[id.x][id.y].clusterIndex.y == -1){
				grid_[id.x][id.y].clusterIndex.x=cen_id.x;
				grid_[id.x][id.y].clusterIndex.y=cen_id.y;
	
				int_two cur_cell_id;
				if(id.x-1>=0){
					cur_cell_id.x=id.x-1; cur_cell_id.y=id.y;
					//if(!grid_cell_->base_cell_[cur_cell_id].outdomain && grid_cell_->base_cell_[cur_cell_id].clusterIndex==-1){
					if((grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 || grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1) && grid_[cur_cell_id.x][cur_cell_id.y].type!=CELL_OUTSIDE ){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
				if(id.y-1>=0){
					cur_cell_id.x=id.x; cur_cell_id.y=id.y-1;
					if((grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 || grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1)  && grid_[cur_cell_id.x][cur_cell_id.y].type!=CELL_OUTSIDE){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
				if(id.x+1<=x_res_-1){
					cur_cell_id.x=id.x+1; cur_cell_id.y=id.y;
					if((grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 || grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1) && grid_[cur_cell_id.x][cur_cell_id.y].type!=CELL_OUTSIDE){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
				if(id.y+1<=y_res_-1){
					cur_cell_id.x=id.x; cur_cell_id.y=id.y+1;
					if((grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.x == -1 || grid_[cur_cell_id.x][cur_cell_id.y].clusterIndex.y == -1) && grid_[cur_cell_id.x][cur_cell_id.y].type!=CELL_OUTSIDE){
						PQ.push(cellPair(grid_[cur_cell_id.x][cur_cell_id.y].point, grid_[cen_id.x][cen_id.y].point, grid_[cen_id.x][cen_id.y].radius, cur_cell_id, cen_id));
					}
				}
			}
		}
	}

	//extract triangles
	void Delaunay::generate_triangle(){
		std::vector<int> faceList;
		int faceNum=0;
		int u, v;
		for(int i=1;i<x_res_;i++){
			for(int j=1;j<y_res_;j++){
					vec2 nearbyCell[4];
					nearbyCell[0]=nearbyCell[1]=nearbyCell[2]=nearbyCell[3]=vec2(-1,-1);
					int activeSize=0;
					bool valid;
					vec2 nei[4] = {vec2(0,0), vec2(-1,0), vec2(-1,-1), vec2(0,-1)}; 
					for(int k=0;k<4;k++){
						int_two baseIndices;
						baseIndices.x = i + nei[k].x;
						baseIndices.y = j + nei[k].y;
						bool noRepeate=true;
						for(int l=0;l<activeSize;l++){
							if(grid_[baseIndices.x][baseIndices.y].clusterIndex.x == nearbyCell[l].x && grid_[baseIndices.x][baseIndices.y].clusterIndex.y == nearbyCell[l].y){
								noRepeate=false;
								break;
							}
						}
						if(noRepeate && grid_[baseIndices.x][baseIndices.y].clusterIndex.x!=-1 && grid_[baseIndices.x][baseIndices.y].clusterIndex.y!=-1){
							nearbyCell[activeSize].x=grid_[baseIndices.x][baseIndices.y].clusterIndex.x;
							nearbyCell[activeSize].y=grid_[baseIndices.x][baseIndices.y].clusterIndex.y;
							activeSize++;
						}
					}
					if(activeSize==3){
						for(int k=0;k<3;k++){
							u=(int)nearbyCell[k].x; v=(int)nearbyCell[k].y;
							faceList.push_back(grid_[u][v].vertexOrder);
						}
						faceNum++;
					}
					if(activeSize==4){
						vec2 triangle[3];
						vec3 curVer,preVer,nextVer,dualVer;
						u=(int)nearbyCell[1].x; v=(int)nearbyCell[1].y;
						curVer.x=grid_[u][v].point.x;
						curVer.y=grid_[u][v].point.y;
						curVer.z=curVer.x*curVer.x + curVer.y*curVer.y;

						u=(int)nearbyCell[0].x; v=(int)nearbyCell[0].y;
						preVer.x=grid_[u][v].point.x;
						preVer.y=grid_[u][v].point.y;
						preVer.z=preVer.x*preVer.x + preVer.y*preVer.y;

						u=(int)nearbyCell[2].x; v=(int)nearbyCell[2].y;
						nextVer.x=grid_[u][v].point.x;
						nextVer.y=grid_[u][v].point.y;
						nextVer.z=nextVer.x*nextVer.x + nextVer.y*nextVer.y;

						u=(int)nearbyCell[3].x; v=(int)nearbyCell[3].y;
						dualVer.x=grid_[u][v].point.x;
						dualVer.y=grid_[u][v].point.y;
						dualVer.z=dualVer.x*dualVer.x + dualVer.y*dualVer.y;

						vec3 vector1(curVer.x-preVer.x,curVer.y-preVer.y,curVer.z-preVer.z);
						vec3 vector2(curVer.x-nextVer.x,curVer.y-nextVer.y,curVer.z-nextVer.z);
						vec3 vector3(curVer.x-dualVer.x,curVer.y-dualVer.y,curVer.z-dualVer.z);

						vec3 crossVector=Geex::cross(vector1,vector2);
						double det=Geex::dot(crossVector,vector3);
						if(det<0){
							triangle[0]=nearbyCell[3];
							triangle[1]=nearbyCell[0];
							triangle[2]=nearbyCell[1];
							for(int k=1;k<4;k++){
								u=(int)nearbyCell[k].x; v=(int)nearbyCell[k].y;
								faceList.push_back(grid_[u][v].vertexOrder);
							}
							faceNum++;
							for(int k=0;k<3;k++){
								u=(int)triangle[k].x; v=(int)triangle[k].y;
								faceList.push_back(grid_[u][v].vertexOrder);
							}
							faceNum++;				
						}
						else{
							triangle[0]=nearbyCell[2];
							triangle[1]=nearbyCell[3];
							triangle[2]=nearbyCell[0];
							for(int k=0;k<3;k++){
								u=(int)nearbyCell[k].x; v=(int)nearbyCell[k].y;
								faceList.push_back(grid_[u][v].vertexOrder);
							}
							faceNum++;
							for(int k=0;k<3;k++){
								u=(int)triangle[k].x; v=(int)triangle[k].y;
								faceList.push_back(grid_[u][v].vertexOrder);
							}
							faceNum++;
						}
					}
				//}
			}
		}

		Geex::MapBuilder builder(triangleMesh_);
		builder.begin_surface();
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE ){
					builder.add_vertex(vec3(grid_[i][j].point.x,grid_[i][j].point.y,0.0));
				}
			}
		}
		for(int j=0;j<faceNum;j++){
			builder.begin_facet();
			builder.add_vertex_to_facet(faceList[3*j]);
			builder.add_vertex_to_facet(faceList[3*j+1]);
			builder.add_vertex_to_facet(faceList[3*j+2]);
			builder.end_facet();
		}
		builder.end_surface();

		MapFacetProperty<vec4> facet_color ;
		facet_color.bind(triangleMesh_, "facet_color") ;

		MapFacetProperty<bool> facet_visit ;
		facet_visit.bind(triangleMesh_, "facet_visit") ;

		MapVertexProperty<int> vertex_order ;
		vertex_order.bind(triangleMesh_, "vertex_order") ;

		FOR_EACH_FACET(Geex::Map,triangleMesh_,fit){
			facet_color[fit]=vec4(color1[0],color1[1],color1[2],color1[3]);
			facet_visit[fit]=false;
		}

		int count=0;
		FOR_EACH_VERTEX(Geex::Map,triangleMesh_,vit){
			vertex_order[vit]=count;
			count++;
		}
	}

    void Delaunay::generate_triangle_vary(){ //This function should be checked carefully, It may be not good written
		
		std::vector<int_two> tri_cells;
		int_two cell;
		int_two neib[2],nei_cell;
	    neib[0].x=-1; neib[0].y=0;
	    neib[1].x=0;  neib[1].y=-1;
		//std::vector<int_two> reset;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].clusterIndex.x!=-1 &&  grid_[i][j].clusterIndex.y!=-1 && grid_[i][j].type != CELL_OUTSIDE){
					bool find=false;
					for(int k=0;k<2;k++){
						nei_cell.x=i+neib[k].x;
						nei_cell.y=j+neib[k].y;
						if(nei_cell.x>=0 && nei_cell.x<x_res_ && nei_cell.y>=0 && nei_cell.y<y_res_){
							if(grid_[nei_cell.x][nei_cell.y].clusterIndex.x!=grid_[i][j].clusterIndex.x ||  grid_[nei_cell.x][nei_cell.y].clusterIndex.y!=grid_[i][j].clusterIndex.y){
								find=true;
							    break;
							}
						}
					}
					if(find){
						cell.x=i;cell.y=j;
					    tri_cells.push_back(cell);
					}
				}
			}
		}

		std::vector<int> faceList;
		int faceNum=0;
		int u, v, u1, v1;
		for(int m=0;m<tri_cells.size();m++){
			u1=tri_cells[m].x;
			v1=tri_cells[m].y;
			if(u1==0 || v1==0){
				continue;
			}
			vec2 nearbyCell[4];
			nearbyCell[0]=nearbyCell[1]=nearbyCell[2]=nearbyCell[3]=vec2(-1,-1);
			int activeSize=0;
			bool valid;
			vec2 nei[4] = {vec2(0,0), vec2(-1,0), vec2(-1,-1), vec2(0,-1)}; 
			for(int k=0;k<4;k++){
				int_two baseIndices;
				baseIndices.x = u1 + nei[k].x;
				baseIndices.y = v1 + nei[k].y;
				bool noRepeate=true;
				for(int l=0;l<activeSize;l++){
					if(grid_[baseIndices.x][baseIndices.y].clusterIndex.x == nearbyCell[l].x && grid_[baseIndices.x][baseIndices.y].clusterIndex.y == nearbyCell[l].y){
						noRepeate=false;
						break;
					}
				}
				if(noRepeate && grid_[baseIndices.x][baseIndices.y].clusterIndex.x!=-1 && grid_[baseIndices.x][baseIndices.y].clusterIndex.y!=-1){
					nearbyCell[activeSize].x=grid_[baseIndices.x][baseIndices.y].clusterIndex.x;
					nearbyCell[activeSize].y=grid_[baseIndices.x][baseIndices.y].clusterIndex.y;
					activeSize++;
				}
			}
			if(activeSize==3){
				for(int k=0;k<3;k++){
					u=(int)nearbyCell[k].x; v=(int)nearbyCell[k].y;
					faceList.push_back(grid_[u][v].vertexOrder);
				}
				faceNum++;
			}
			if(activeSize==4){
				vec2 triangle[3];
				vec3 curVer,preVer,nextVer,dualVer;
				u=(int)nearbyCell[1].x; v=(int)nearbyCell[1].y;
				curVer.x=grid_[u][v].point.x;
				curVer.y=grid_[u][v].point.y;
				curVer.z=curVer.x*curVer.x + curVer.y*curVer.y;

				u=(int)nearbyCell[0].x; v=(int)nearbyCell[0].y;
				preVer.x=grid_[u][v].point.x;
				preVer.y=grid_[u][v].point.y;
				preVer.z=preVer.x*preVer.x + preVer.y*preVer.y;

				u=(int)nearbyCell[2].x; v=(int)nearbyCell[2].y;
				nextVer.x=grid_[u][v].point.x;
				nextVer.y=grid_[u][v].point.y;
				nextVer.z=nextVer.x*nextVer.x + nextVer.y*nextVer.y;

				u=(int)nearbyCell[3].x; v=(int)nearbyCell[3].y;
				dualVer.x=grid_[u][v].point.x;
				dualVer.y=grid_[u][v].point.y;
				dualVer.z=dualVer.x*dualVer.x + dualVer.y*dualVer.y;

				vec3 vector1(curVer.x-preVer.x,curVer.y-preVer.y,curVer.z-preVer.z);
				vec3 vector2(curVer.x-nextVer.x,curVer.y-nextVer.y,curVer.z-nextVer.z);
				vec3 vector3(curVer.x-dualVer.x,curVer.y-dualVer.y,curVer.z-dualVer.z);

				vec3 crossVector=Geex::cross(vector1,vector2);
				double det=Geex::dot(crossVector,vector3);
				if(det<0){
					triangle[0]=nearbyCell[3];
					triangle[1]=nearbyCell[0];
					triangle[2]=nearbyCell[1];
					for(int k=1;k<4;k++){
						u=(int)nearbyCell[k].x; v=(int)nearbyCell[k].y;
						faceList.push_back(grid_[u][v].vertexOrder);
					}
					faceNum++;
					for(int k=0;k<3;k++){
						u=(int)triangle[k].x; v=(int)triangle[k].y;
						faceList.push_back(grid_[u][v].vertexOrder);
					}
					faceNum++;				
				}
				else{
					triangle[0]=nearbyCell[2];
					triangle[1]=nearbyCell[3];
					triangle[2]=nearbyCell[0];
					for(int k=0;k<3;k++){
						u=(int)nearbyCell[k].x; v=(int)nearbyCell[k].y;
						faceList.push_back(grid_[u][v].vertexOrder);
					}
					faceNum++;
					for(int k=0;k<3;k++){
						u=(int)triangle[k].x; v=(int)triangle[k].y;
						faceList.push_back(grid_[u][v].vertexOrder);
					}
					faceNum++;
				}
			}
		}

		Geex::MapBuilder builder(triangleMesh_);
		builder.begin_surface();
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE ){
					builder.add_vertex(vec3(grid_[i][j].point.x,grid_[i][j].point.y,0.0));
				}
			}
		}
		for(int j=0;j<faceNum;j++){
			builder.begin_facet();
			builder.add_vertex_to_facet(faceList[3*j]);
			builder.add_vertex_to_facet(faceList[3*j+1]);
			builder.add_vertex_to_facet(faceList[3*j+2]);
			builder.end_facet();
		}
		builder.end_surface();

		MapFacetProperty<vec4> facet_color ;
		facet_color.bind(triangleMesh_, "facet_color") ;

		MapFacetProperty<bool> facet_visit ;
		facet_visit.bind(triangleMesh_, "facet_visit") ;

		MapVertexProperty<int> vertex_order ;
		vertex_order.bind(triangleMesh_, "vertex_order") ;

		FOR_EACH_FACET(Geex::Map,triangleMesh_,fit){
			facet_color[fit]=vec4(color1[0],color1[1],color1[2],color1[3]);
			facet_visit[fit]=false;
		}

		int count=0;
		FOR_EACH_VERTEX(Geex::Map,triangleMesh_,vit){
			vertex_order[vit]=count;
			count++;
		}
	}

	//edge flipping
	void Delaunay::show_original_triangle(){
		MapFacetProperty<vec4> facet_color ;
		facet_color.bind(triangleMesh_, "facet_color") ;
		FOR_EACH_FACET(Geex::Map,triangleMesh_,fit){
			facet_color[fit]=vec4(color1[0],color1[1],color1[2],color1[3]);
		}
	}

	void Delaunay::show_non_delanuary(){
		int non_delanuary=0;

		MapFacetProperty<vec4> facet_color ;
		facet_color.bind(triangleMesh_, "facet_color") ;

		MapFacetProperty<vec2> circumcircle_center ;
		circumcircle_center.bind(triangleMesh_, "circumcircle_center") ;

		MapFacetProperty<double> circumcircle_rad ;
		circumcircle_rad.bind(triangleMesh_, "circumcircle_rad") ;

		MapFacetProperty<bool> non_delanuary_facet ;
		non_delanuary_facet.bind(triangleMesh_, "non_delanuary_facet") ;

		MapHalfedgeProperty<bool> flip_edge ;
		flip_edge.bind(triangleMesh_, "flip_edge") ;

		//for each edge, use vector cross product to detect its empty-circle property
		FOR_EACH_EDGE(Geex::Map,triangleMesh_,it){
			Geex::MapCombels::Halfedge* opposite=it->opposite();
			flip_edge[it]=false;
			flip_edge[opposite]=false;
			if(!it->is_border() && !opposite->is_border()){
				Geex::MapCombels::Facet* cur_facet=it->facet();
				Geex::MapCombels::Facet* adj_facet=opposite->facet();
				vec2 circle_center;
				double x1=it->vertex()->point().x, y1=it->vertex()->point().y,
					x2=it->next()->vertex()->point().x,y2=it->next()->vertex()->point().y,
					x3=it->prev()->vertex()->point().x,y3=it->prev()->vertex()->point().y,
					x4=opposite->next()->vertex()->point().x, y4=opposite->next()->vertex()->point().y;
				circle_center.x=((y2-y1)*(y3*y3-y1*y1+x3*x3-x1*x1)-(y3-y1)*(y2*y2-y1*y1+x2*x2-x1*x1))/(2*(x3-x1)*(y2-y1)-2*((x2-x1)*(y3-y1)));
				circle_center.y=((x2-x1)*(x3*x3-x1*x1+y3*y3-y1*y1)-(x3-x1)*(x2*x2-x1*x1+y2*y2-y1*y1))/(2*(y3-y1)*(x2-x1)-2*((y2-y1)*(x3-x1)));
				double circle_radius=(circle_center.x-x1)*(circle_center.x-x1)+(circle_center.y-y1)*(circle_center.y-y1);
				double dis=(circle_center.x-x4)*(circle_center.x-x4)+(circle_center.y-y4)*(circle_center.y-y4);
				circumcircle_center[cur_facet]=circle_center;
				circumcircle_rad[cur_facet]=sqrt(circle_radius);
				if(dis<circle_radius){
					facet_color[cur_facet]=vec4(color2[0],color2[1],color2[2],color2[3]);
					facet_color[adj_facet]=vec4(color2[0],color2[1],color2[2],color2[3]);

					non_delanuary+=2;
					flip_edge[it]=true;
					flip_edge[opposite]=true;
				}
			}
		}
		double ratio=(double)(triangleMesh_->nb_facets()-non_delanuary)/triangleMesh_->nb_facets();
		std::cerr<<"Ratio of delanuay triangles: "<<ratio<<std::endl;
	}

	void Delaunay::edge_flip(){
		
		MapFacetProperty<vec4> facet_color ;
		facet_color.bind(triangleMesh_, "facet_color") ;

		MapHalfedgeProperty<bool> flip_edge ;
		flip_edge.bind(triangleMesh_, "flip_edge") ;

		std::vector<Geex::MapCombels::Halfedge*> heList;
		std::vector<Geex::MapCombels::Vertex*> verList;
		int neighbor=0;

		FOR_EACH_EDGE(Geex::Map,triangleMesh_,it){
			if(flip_edge[it]){
				Geex::MapCombels::Halfedge* opposite=it->opposite();

				heList.push_back(it);
				heList.push_back(opposite);

				verList.push_back(it->next()->vertex());
				verList.push_back(it->prev()->vertex());
				verList.push_back(opposite->next()->vertex());
				verList.push_back(it->vertex());

				if(flip_edge[it->prev()] || flip_edge[it->next()]){
					flip_edge[it->prev()]=flip_edge[it->next()]=false;
					flip_edge[it->prev()->opposite()]=flip_edge[it->next()->opposite()]=false;
					neighbor++;
				}
				if(flip_edge[opposite->prev()] || flip_edge[opposite->next()]){
					flip_edge[opposite->prev()]=flip_edge[opposite->next()]=false;
					flip_edge[opposite->prev()->opposite()]=flip_edge[opposite->next()->opposite()]=false;
					neighbor++;
				}
			}
		}

		Geex::MapEditor editor(triangleMesh_);
		if(verList.size()>0){
			for(int i=0;i<=verList.size()-4;i+=4){
				Geex::MapCombels::Halfedge* h1=editor.make_triangle(verList[i]->point(),verList[i+1]->point(),verList[i+2]->point());
				Geex::MapCombels::Halfedge* h2=editor.make_triangle(verList[i+2]->point(),verList[i+3]->point(),verList[i]->point());
			}
		}
		
		for(int j=0;j<heList.size();j++){
			if(heList[j]){
				if(!heList[j]->is_border()){
					editor.erase_facet(heList[j]);
				}
				
			}
		}

		FOR_EACH_FACET(Geex::Map,triangleMesh_,fit){
			facet_color[fit]=vec4(color1[0],color1[1],color1[2],color1[3]);
		}
		FOR_EACH_HALFEDGE(Geex::Map,triangleMesh_,eit){
			flip_edge[eit]=false;
		}

		if(neighbor>0){
			show_non_delanuary();
			edge_flip();
		}
	}

	double Delaunay::determinant(double a[], int n)
	{
		if(n==9){
			return a[0]*a[4]*a[8]+a[1]*a[5]*a[6]+a[2]*a[3]*a[7]-a[2]*a[4]*a[6]-a[0]*a[5]*a[7]-a[1]*a[3]*a[8];
		}
	}

	void Delaunay::show_non_rt(){

		MapFacetProperty<vec4> facet_color ;
		facet_color.bind(triangleMesh_, "facet_color") ;

		MapFacetProperty<bool> non_delanuary_facet ;
		non_delanuary_facet.bind(triangleMesh_, "non_delanuary_facet") ;

		MapHalfedgeProperty<bool> flip_edge ;
		flip_edge.bind(triangleMesh_, "flip_edge") ;

		MapVertexProperty<double> ver_weight ;
		ver_weight.bind(triangleMesh_, "ver_weight") ;

		FOR_EACH_VERTEX(Geex::Map,triangleMesh_,vit){
			int u=vit->point().x/gridsize_;
			int v=vit->point().y/gridsize_;
			ver_weight[vit]=SQR(grid_[u][v].radius);
		}

		FOR_EACH_EDGE(Geex::Map,triangleMesh_,it){
			Geex::MapCombels::Halfedge* opposite=it->opposite();
			flip_edge[it]=false;
			flip_edge[opposite]=false;
			if(!it->is_border() && !opposite->is_border()){
				Geex::MapCombels::Facet* cur_facet=it->facet();
				Geex::MapCombels::Facet* adj_facet=opposite->facet();
				vec2 circle_center;
				double x1=it->vertex()->point().x, y1=it->vertex()->point().y, w1=ver_weight[it->vertex()],
					x2=it->next()->vertex()->point().x,y2=it->next()->vertex()->point().y, w2=ver_weight[it->next()->vertex()],
					x3=it->prev()->vertex()->point().x,y3=it->prev()->vertex()->point().y, w3=ver_weight[it->prev()->vertex()],
					x4=opposite->next()->vertex()->point().x, y4=opposite->next()->vertex()->point().y, w4=ver_weight[opposite->next()->vertex()];
				double a1[9]={x2,y2,SQR(x2)+SQR(y2)-w2, x3,y3,SQR(x3)+SQR(y3)-w3, x4,y4,SQR(x4)+SQR(y4)-w4};
				double a2[9]={x1,y1,SQR(x1)+SQR(y1)-w1, x3,y3,SQR(x3)+SQR(y3)-w3, x4,y4,SQR(x4)+SQR(y4)-w4};
				double a3[9]={x1,y1,SQR(x1)+SQR(y1)-w1, x2,y2,SQR(x2)+SQR(y2)-w2, x4,y4,SQR(x4)+SQR(y4)-w4};
				double a4[9]={x1,y1,SQR(x1)+SQR(y1)-w1, x2,y2,SQR(x2)+SQR(y2)-w2, x3,y3,SQR(x3)+SQR(y3)-w3};
				double det=determinant(a1,9) - determinant(a2,9) + determinant(a3,9) - determinant(a4,9);
				if(det<0){
					facet_color[cur_facet]=vec4(color2[0],color2[1],color2[2],color2[3]);
					facet_color[adj_facet]=vec4(color2[0],color2[1],color2[2],color2[3]);

					flip_edge[it]=true;
					flip_edge[opposite]=true;
				}
			}
		}
	}

	//save or load functions
	void Delaunay::save_ps_new(const std::string& filename){
		std::string domain_filename = filename ;
		std::string grid_filename = filename ;
		std::string pts_filename = filename ;
		std::string pow_filename = filename ;
		std::string pow_cir_filename = filename ;
		std::string cluster_filename = filename ;
		std::string mesh_filename = filename ;
		std::string model_filename = filename ;
		//std::string vertex_filename = filename ;

		domain_filename.insert(filename.size()-3, "_domain") ;
		grid_filename.insert(filename.size()-3, "_grid") ;
		pts_filename.insert(filename.size()-3, "_pts") ;
		pow_filename.insert(filename.size()-3, "_pow") ;
		pow_cir_filename.insert(filename.size()-3, "_pow_cir") ;
		cluster_filename.insert(filename.size()-3, "_cluster") ;
		mesh_filename.insert(filename.size()-3, "_mesh") ;
		model_filename.insert(filename.size()-3, "_model") ;
		//vertex_filename.insert(filename.size()-3, "_vertex") ;

		double scale = 1000 ;
		char* translate="2 2 translate";

		//save domain
		std::ofstream out(domain_filename.c_str()) ;
		out<<translate<<std::endl;
		for(unsigned int i=0; i<boundary_.size(); i++) {
			out << "newpath" << std::endl ;
			out << boundary_[i].vertex[0].x*scale << " " << boundary_[i].vertex[0].y*scale << " moveto" << std::endl ;
			out << boundary_[i].vertex[1].x*scale << " " << boundary_[i].vertex[1].y*scale << " lineto" << std::endl ;
			out << "closepath" << std::endl ;
		    out << "3.0 setlinewidth" << std::endl ;
		    out << "0.0 0.0 0.5 setrgbcolor" << std::endl ;
		    out << "stroke" << std::endl ;
        }
		out<<"showpage"<<std::endl;
		out.close() ;

		//save grid cells
		out.open(grid_filename.c_str()) ;
		out<<translate<<std::endl;
		for(int i=0; i<x_res_; ++i) {
			for(int j=0; j<y_res_; ++j) {
					out << "newpath" << std::endl ;
					out << i*gridsize_*scale << " " << j*gridsize_*scale << " moveto" << std::endl ;
					out << (i+1)*gridsize_*scale << " " << j*gridsize_*scale << " lineto" << std::endl ;
					out << (i+1)*gridsize_*scale << " " << (j+1)*gridsize_*scale << " lineto" << std::endl ;
					out << i*gridsize_*scale << " " << (j+1)*gridsize_*scale << " lineto" << std::endl ;
					out << "closepath" << std::endl ;
					/*if(grid_[i][j].valid && grid_[i][j].type!=CELL_OUTSIDE ){
						out << "gsave" << std::endl ;
						out << "0.95 0.95 0.012 setrgbcolor" << std::endl ;
						out << "fill" << std::endl ;
						out << "grestore" << std::endl ;
					}*/
					out << "0.5 setlinewidth" << std::endl ;
					out << "stroke" << std::endl ;
			}
		}
		out<<"showpage"<<std::endl;
		out.close() ;

		//save clusters
		out.open(cluster_filename.c_str()) ;
		out<<translate<<std::endl;
		for(int i=0; i<x_res_; ++i) {
			for(int j=0; j<y_res_; ++j) {
				    int_two index;
					index.x=grid_[i][j].clusterIndex.x;
					index.y=grid_[i][j].clusterIndex.y;
					if(index.x!=-1 && index.y!=-1){
						out << "newpath" << std::endl ;
						out << i*gridsize_*scale << " " << j*gridsize_*scale << " moveto" << std::endl ;
						out << (i+1)*gridsize_*scale << " " << j*gridsize_*scale << " lineto" << std::endl ;
						out << (i+1)*gridsize_*scale << " " << (j+1)*gridsize_*scale << " lineto" << std::endl ;
						out << i*gridsize_*scale << " " << (j+1)*gridsize_*scale << " lineto" << std::endl ;
						out << "closepath" << std::endl ;
						out << "gsave" << std::endl ;
						out << cluster_colors_[index.x][index.y].x<<" "<<cluster_colors_[index.x][index.y].y<<" "<<cluster_colors_[index.x][index.y].z <<" setrgbcolor" << std::endl ;
						out << "fill" << std::endl ;
						out << "grestore" << std::endl ;
						out << "0.7 setlinewidth" << std::endl ;
						out << "stroke" << std::endl ;
					}
			}
		}
		for(int i=0; i<x_res_; ++i) {
			for(int j=0; j<y_res_; ++j) {
				    int_two index;
					index.x=grid_[i][j].clusterIndex.x;
					index.y=grid_[i][j].clusterIndex.y;
					if(index.x!=-1 && index.y!=-1){
						out << "newpath" << std::endl ;
						out << (i+0.5)*gridsize_*scale << " " << (j+0.5)*gridsize_*scale << " moveto" << std::endl ;
						out << grid_[index.x][index.y].point.x*scale << " " << grid_[index.x][index.y].point.y*scale << " lineto" << std::endl ;
						out << "closepath" << std::endl ;
						out << "1.0 1.0 1.0 setrgbcolor" << std::endl ;
						out << "0.4 setlinewidth" << std::endl ;
						out << "stroke" << std::endl ;
					}
			}
		}
		out<<"showpage"<<std::endl;
		out.close() ;

		// save points
		out.open(pts_filename.c_str()) ;
		out<<translate<<std::endl;
		for(int i=0; i<samples_.size(); ++i) {
			out << samples_[i][0]*scale << " " << samples_[i][1]*scale << " " ;
			out << 0.7 << " 0 360 arc closepath" << std::endl ;
			out << "0.2 0.55 1.0 setrgbcolor" << std::endl ;
			out << "fill" << std::endl ;
		}
		//out<<"showpage"<<std::endl;
		out.close() ;

		// save power circle
		out.open(pow_cir_filename.c_str()) ;
		out<<translate<<std::endl;
		for(int i=0; i<samples_.size(); ++i) {
			out << samples_[i][0]*scale << " " << samples_[i][1]*scale << " " ;
			out << sqrt(weights_[i])*scale << " 0 360 arc closepath" << std::endl ;
			out << "0.5 setlinewidth" << std::endl ;
			out << "0.5 0.5 0.5 setrgbcolor" << std::endl ;
			out << "stroke" << std::endl ;
		}
		//out<<"showpage"<<std::endl;
		out.close() ;

		// save power fill
		out.open(pow_filename.c_str()) ;
		out<<"2 2 translate"<<std::endl;
		for(int i=0; i<samples_.size(); ++i) {
			out << samples_[i][0]*scale << " " << samples_[i][1]*scale << " " ;
			out << sqrt(weights_[i])*scale << " 0 360 arc closepath" << std::endl ;
			//out << "gsave" << std::endl ;
			out << "0.8 0.8 0.3 setrgbcolor" << std::endl ;
			out << "fill" << std::endl ;
			/*out << "grestore" << std::endl ;
			out << "0.8 0.8 0.8 setrgbcolor" << std::endl ;
			out << "0.5 setlinewidth" << std::endl ;
			out << "stroke" << std::endl ;*/
		}
		out.close() ;

		//save triangles
		MapFacetProperty<vec4> facet_color ;
		facet_color.bind(triangleMesh_, "facet_color") ;

		MapHalfedgeProperty<bool> flip_edge ;
		flip_edge.bind(triangleMesh_, "flip_edge") ;

		static float color1[4] = {
		0.0f, 0.5f, 0.0f, 0.3f
	    } ;
		out.open(mesh_filename.c_str()) ;
		out<<translate<<std::endl;
		FOR_EACH_FACET(Geex::Map,triangleMesh_,it){
			if(facet_color[it].x==0.0){
				Geex::MapCombels::Halfedge* thisHe=it->halfedge();
				Geex::MapCombels::Halfedge* prevHe=thisHe->prev();
				Geex::MapCombels::Halfedge* nextHe=thisHe->next();
				out << "newpath" << std::endl ;
				out << prevHe->vertex()->point().x*scale << " " << prevHe->vertex()->point().y*scale << " moveto" << std::endl ;
				out << thisHe->vertex()->point().x*scale << " " << thisHe->vertex()->point().y*scale << " lineto" << std::endl ;
				out << nextHe->vertex()->point().x*scale << " " << nextHe->vertex()->point().y*scale << " lineto" << std::endl ;
				out << "closepath" << std::endl ;
				//out << facet_color[it].x<<" "<< facet_color[it].y <<" "<< facet_color[it].z <<" setrgbcolor" << std::endl ;
				out << "0.013 0.600 0.028 setrgbcolor" << std::endl ;
				out << "0.6 setlinewidth" << std::endl ;
				out << "stroke" << std::endl ;
			}
		}
		FOR_EACH_FACET(Geex::Map,triangleMesh_,it){
			if(facet_color[it].x!=0.0){
				Geex::MapCombels::Halfedge* thisHe=it->halfedge();
				Geex::MapCombels::Halfedge* prevHe=thisHe->prev();
				Geex::MapCombels::Halfedge* nextHe=thisHe->next();
				out << "newpath" << std::endl ;
				out << prevHe->vertex()->point().x*scale << " " << prevHe->vertex()->point().y*scale << " moveto" << std::endl ;
				out << thisHe->vertex()->point().x*scale << " " << thisHe->vertex()->point().y*scale << " lineto" << std::endl ;
				out << nextHe->vertex()->point().x*scale << " " << nextHe->vertex()->point().y*scale << " lineto" << std::endl ;
				out << "closepath" << std::endl ;
				//out << facet_color[it].x<<" "<< facet_color[it].y <<" "<< facet_color[it].z <<" setrgbcolor" << std::endl ;
				out << "0.85 0.12 0.12 setrgbcolor" << std::endl ;
				out << "0.6 setlinewidth" << std::endl ;
				out << "stroke" << std::endl ;
			}
		}
		FOR_EACH_EDGE(Geex::Map,triangleMesh_,eit){
			if(flip_edge[eit]){
				out << "newpath" << std::endl ;
				out << eit->vertex()->point().x*scale << " " << eit->vertex()->point().y*scale << " moveto" << std::endl ;
				out << eit->opposite()->vertex()->point().x*scale << " " << eit->opposite()->vertex()->point().y*scale << " lineto" << std::endl ;
				out << "closepath" << std::endl ;
				out << "0.1 0.8 0.8 setrgbcolor" << std::endl ;
				out << "0.6 setlinewidth" << std::endl ;
				out << "stroke" << std::endl ;
			}
		}
		out<<"showpage"<<std::endl;
		out.close() ;
		//save model
		std::ofstream ver_out("vertex.txt") ;
		MapVertexProperty<int> vertex_order ;
		vertex_order.bind(triangleMesh_, "vertex_order") ;
		int count=1;
		out.open(model_filename.c_str()) ;
		FOR_EACH_VERTEX(Geex::Map,triangleMesh_,vit){
			vertex_order[vit]=count;
			out<<"v "<<vit->point().x<<" "<<vit->point().y<<" 0.0"<<std::endl;
			ver_out<<vit->point().x<<" "<<vit->point().y<<std::endl;
			count++;
		}
		FOR_EACH_FACET(Geex::Map,triangleMesh_,fit){
			Geex::MapCombels::Halfedge* thisHe=fit->halfedge();
			Geex::MapCombels::Halfedge* prevHe=thisHe->prev();
			Geex::MapCombels::Halfedge* nextHe=thisHe->next();
			out<<"f "<<vertex_order[thisHe->vertex()]<<" "<<vertex_order[nextHe->vertex()]<<" "<<vertex_order[prevHe->vertex()]<<std::endl;
		}
		ver_out.close();
		out.close();
		
	}

	void Delaunay::save_ps(const std::string& filename){
		std::string domain_filename = filename ;
		std::string mesh_filename = filename ;
		std::string model_filename = filename ;

		domain_filename.insert(filename.size()-3, "_domain") ;
		mesh_filename.insert(filename.size()-3, "_mesh") ;

		double scale = 300 ;
		char* translate="2 2 translate";

		//save domain
		std::ofstream out(domain_filename.c_str()) ;
		out<<translate<<std::endl;
		for(unsigned int i=0; i<boundary_.size(); i++) {
			out << "newpath" << std::endl ;
			out << boundary_[i].vertex[0].x*scale << " " << boundary_[i].vertex[0].y*scale << " moveto" << std::endl ;
			out << boundary_[i].vertex[1].x*scale << " " << boundary_[i].vertex[1].y*scale << " lineto" << std::endl ;
			out << "closepath" << std::endl ;
		    out << "3.0 setlinewidth" << std::endl ;
		    out << "0.0 0.0 0.5 setrgbcolor" << std::endl ;
		    out << "stroke" << std::endl ;
        }
		out<<"showpage"<<std::endl;
		out.close() ;

		//save triangles
		out.open(mesh_filename.c_str()) ;
		out<<translate<<std::endl;
		FOR_EACH_FACET(Geex::Map,triangleMesh_,it){
			Geex::MapCombels::Halfedge* thisHe=it->halfedge();
			Geex::MapCombels::Halfedge* prevHe=thisHe->prev();
			Geex::MapCombels::Halfedge* nextHe=thisHe->next();
			out << "newpath" << std::endl ;
			out << prevHe->vertex()->point().x*scale << " " << prevHe->vertex()->point().y*scale << " moveto" << std::endl ;
			out << thisHe->vertex()->point().x*scale << " " << thisHe->vertex()->point().y*scale << " lineto" << std::endl ;
			out << nextHe->vertex()->point().x*scale << " " << nextHe->vertex()->point().y*scale << " lineto" << std::endl ;
			out << "closepath" << std::endl ;
			out << "0.013 0.600 0.028 setrgbcolor" << std::endl ;
			out << "0.6 setlinewidth" << std::endl ;
			out << "stroke" << std::endl ;
		}
		out<<"showpage"<<std::endl;
		out.close() ;

		std::cerr << "Save ps to: " << filename << std::endl ;

	}

	void Delaunay::load(const std::string& filename) {
        std::cerr << "loading " << filename << std::endl ;
        clear() ;
        std::ifstream in(filename.c_str()) ;
        if(!in) {
            std::cerr << "could not open file" << std::endl ;
            return ;
        }

		vec2 p ;
		double w, min_rad ;
		char sample_mode = 'U';
		int count = 0;
		begin_insert() ;
        while(in) {
			if(count == 0){
				in >> min_rad >> sample_mode;
			}
			else{
				in >> p >> w ;
				if(in) {  // we need to do this additional check else we got the last point twice !
					insert(p, w) ; 
				}
			}
			count++;
        } ;
		end_insert(true) ;
		in.close() ;

		radius_min_ = min_rad;
		update_grid(radius_min_) ;
		if(sample_mode == 'U'){
			weight_scheme_ = WT_UNIFORM;
		}
		if(sample_mode == 'A'){
			weight_scheme_ = WT_LFS;
		}
		// update weight
		for(unsigned int i=0; i<samples_.size(); ++i) {
			int u = samples_[i].x/gridsize_ ;
			int v = samples_[i].y/gridsize_ ;
			grid_[u][v].sampled = true ;
			grid_[u][v].point  = samples_[i] ;
			grid_[u][v].radius=sqrt(weights_[i]);
			mark_hit_cells(u, v, gridsize_, samples_[i], weights_[i], i) ;
		}
		//new added
		unsigned int order=0;
		for(int i=0;i<x_res_;i++){
			for(int j=0;j<y_res_;j++){
				if(grid_[i][j].sampled && grid_[i][j].type != CELL_OUTSIDE ){
					grid_[i][j].vertexOrder=order;
					order++;
					sample_number_++;
				}
			}
		}
		std::cerr<<std::endl ;
		std::cerr << "Sample radius: " << radius_min_ <<std::endl ;
		std::cerr << "Number of sampled points: " << sample_number_ <<std::endl ;
    }

	void Delaunay::save(const std::string& filename){
		char sampling_mode;
		if(weight_scheme_ == WT_UNIFORM){
			sampling_mode = 'U';
		}
		if(weight_scheme_ == WT_LFS){
			sampling_mode = 'A';
		}

		std::ofstream file(filename.c_str());
		file<< radius_min_ << " " << sampling_mode << std::endl;
		for(int i=0; i<samples_.size(); ++i) {
			file <<samples_[i][0] << " " << samples_[i][1] << " " << weights_[i] << std::endl ;
		}
		file.close();
		std::cerr << "Save pts to: " << filename << std::endl ;
	}
}
