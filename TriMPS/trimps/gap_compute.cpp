#include "gap_compute.h"

namespace Geex {

#ifndef SQR(x)
#define SQR(x) ((x)*(x))
#endif 

	//
	void GapCompute::compute_face_gap(Delaunay::Face_handle f, std::vector<FaceGap>& fgaps) {
		vec2   v[3] = {to_geex(f->vertex(0)->point()), to_geex(f->vertex(1)->point()), to_geex(f->vertex(2)->point())} ;
		double w[3] = {f->vertex(0)->point().weight(), f->vertex(1)->point().weight(), f->vertex(2)->point().weight()} ;
		double r[3] = {f->vertex(0)->radius, f->vertex(1)->radius, f->vertex(2)->radius} ;
		vec2   e[3] = {v[2]-v[1], v[0]-v[2], v[1]-v[0]} ;
		double l[3] = {e[0].length(), e[1].length(), e[2].length()} ;
		FaceGap gap ;

		for(int i=0; i<3; ++i) { 
			e[i] /= l[i] ; 
		}		

		gx_assert(f->has_gap) ;

		for(int i=0; i<3; ++i) {
			Delaunay::Face_handle& fo = f->neighbor(i) ;
			int j = f->ccw(i) ;
			int k = f->cw(i) ;

			// we don't handle this case for face gap, compute vertex gap instead
			if(delaunay_->is_infinite(fo) && l[i]>r[j]+r[k])
				return ;

			if(!fo->has_gap) {
				double a = 0.5*(l[i]+(w[j]-w[k])/l[i]) ;
				double b = l[i]-a ;
				double h = sqrt(w[j]-a*a) ;
				double h2 =sqrt(w[k]-b*b) ; 
				gx_assert(fabs(h-h2)<1e-10) ;
				vec2 p = (a/l[i]*v[k] + b/l[i]*v[j]) + h*vec2(-e[i].y, e[i].x) ;
				gap.push_back(p) ;
			}
			else {
				Line<double> L(v[j], vec2(-e[i].y, e[i].x)) ;
				// the current edge is covered
				if(l[i]<r[j]+r[k]) {
					// same as case 1: two gaps can't see each other
					if(L.side(f->dual)>0) { // acute 
						if(L.side(fo->dual)<0) { // acute neighbor
							double a = 0.5*(l[i]+(w[j]-w[k])/l[i]) ;
							double b = l[i]-a ;
							double h = sqrt(w[j]-a*a) ;
							double h2 =sqrt(w[k]-b*b) ; 
							gx_assert(fabs(h-h2)<1e-10) ;
							vec2 p = (a/l[i]*v[k] + b/l[i]*v[j]) + h*vec2(-e[i].y, e[i].x) ;
							gap.push_back(p) ;
						} 
						else {
							vec2 dualj = fo->dual - v[j] ;
							vec2 dualk = fo->dual - v[k] ;
							gap.push_back(v[j] + r[j]*normalize(dualj)) ;
							gap.push_back(v[k] + r[k]*normalize(dualk)) ;
						}
					}
					// two gaps see each other, obtuse triangle edge
					else {
						vec2 dualj = f->dual - v[j] ;
						vec2 dualk = f->dual - v[k] ;
						gap.push_back(v[j] + r[j]*normalize(dualj)) ;
						gap.push_back(v[k] + r[k]*normalize(dualk)) ;
					}
				} 
				else {
					if(L.side(f->dual)>0) { // acute triangle
						if(L.side(fo->dual)<0) { // acute neighbor
							gap.push_back(v[j] + r[j]*e[i]) ;
							gap.push_back(v[k] - r[k]*e[i]) ;
						}
						else { // obtuse neighbor
							vec2 dualj = fo->dual - v[j] ;
							vec2 dualk = fo->dual - v[k] ;
							gap.push_back(v[j] + r[j]*normalize(dualj)) ;
							gap.push_back(v[k] + r[k]*normalize(dualk)) ;
						// gap.add_disk(vo) ;
						} 
					} 
					else { // obtuse triangle
						vec2 dualj = f->dual - v[j] ;
						vec2 dualk = f->dual - v[k] ;
						gap.push_back(v[j] + r[j]*normalize(dualj)) ;
						gap.push_back(v[k] + r[k]*normalize(dualk)) ;
					}
				}
			}
		}

		if(gap.size()>=3) {
			for(int i=0; i<3; ++i) {
				gap.add_disk(to_geex(f->vertex(i)->point()), f->vertex(i)->point().weight()) ;
			}
			fgaps.push_back(gap) ;
		}
	}

	//
	void GapCompute::compute_face_vertex_gap(Delaunay::Face_handle f, std::vector<FaceGap>& fgaps) {
	}

	//
	void GapCompute::compute_vertex_gaps(Delaunay::Vertex_handle v, std::vector<VertexGap>& vgaps) {
		if(v->dual_intersects_boundary) {
			compute_vertex_gaps_clipped(v, vgaps) ;
		} 
		else {
			compute_vertex_gaps_inner(v, vgaps) ;
		}
	}

	//
	void GapCompute::compute_vertex_gaps_inner(Delaunay::Vertex_handle v, std::vector<VertexGap>& vgaps) {
		Delaunay::Face_circulator f1 = delaunay_->incident_faces(v) ;
		Delaunay::Face_circulator f2 = f1 ; f2++ ;
		vec2 p0 = to_geex(v->point()) ;
		double area = 0 ;

		do {
			if(f1->has_gap && f2->has_gap) {
				vec2 p1 = f1->dual ;
				vec2 p2 = f2->dual ;
				Delaunay::Vertex_handle vopp = f1->vertex(f1->cw(f1->index(v))) ;
				vec2 p3 = to_geex(vopp->point()) ;
				vec2 dir = normalize(p3-p0) ;
				vec2 nor = vec2(-dir.y, dir.x) ;
				Line<double> L(p0, nor) ;

				// p1 p2 on same side
				if(L.side(p1) * L.side(p2) > 0) {
					VertexGap gap ;
					gap.set_disk(p0, v->radius) ;
					gap.push_back(p0+v->radius*normalize(p1-p0)) ;
					gap.push_back(p1) ;
					gap.push_back(p2) ;
					gap.push_back(p0+v->radius*normalize(p2-p0)) ;
					gap.gap_index() = f1->gap_index ;
					vgaps.push_back(gap) ;
//					double garea = compute_sector_gap_area(p0, p1, p2, radius) ;
//					area += garea ;
				}
				else {
					double d03 = distance(p0, p3) ;
					// TODO: this case can be split into two
					if(d03>v->radius+vopp->radius) {
						vec2 mp = delaunay_->get_mid_point(v->point(), vopp->point()) ; //0.5*(p0+p3) ;
						vec2 p03 = p0 + v->radius*dir ;
						vec2 p01 = p0 + v->radius*normalize(p1-p0) ;
						vec2 p02 = p0 + v->radius*normalize(p2-p0) ;

						VertexGap gap0, gap1 ;

						gap0.set_disk(p0, v->radius) ;
						gap0.push_back(p03) ;
						gap0.push_back(p0+v->radius*normalize(p1-p0)) ;
						gap0.push_back(p1) ;
						gap0.push_back(mp) ;
						gap0.gap_index() = f1->gap_index ;
						vgaps.push_back(gap0) ;

						gap1.set_disk(p0, v->radius) ;
						gap1.push_back(p0+v->radius*normalize(p2-p0)) ;
						gap1.push_back(p03) ;
						gap1.push_back(mp) ;
						gap1.push_back(p2) ;
						gap1.gap_index() = f2->gap_index ;
						vgaps.push_back(gap1) ;
						//double garea = compute_sector_gap_area(p0, p1, p2, exclude_radius) ;
						//area += garea ;
					} else {
						vec2 mp = delaunay_->get_mid_point(v->point(), vopp->point()) ; // 0.5*(p0+p3) ;
						double h = sqrt(SQR(v->radius)-distance2(p0, mp)) ; //SQR(d03*0.5)) ;
						vec2 p12 = mp - h*nor;
						vec2 p21 = mp + h*nor;

						VertexGap gap1 ;
						gap1.set_disk(p0, v->radius) ;
						gap1.push_back(p0+v->radius*normalize(p1-p0)) ;
						gap1.push_back(p1) ;
						gap1.push_back(p12) ;
						gap1.gap_index() = f1->gap_index ;
						vgaps.push_back(gap1) ;

						VertexGap gap2 ;
						gap2.set_disk(p0, v->radius) ;
						gap2.push_back(p0+v->radius*normalize(p2-p0)) ;
						gap2.push_back(p21) ;
						gap2.push_back(p2) ;
						gap2.gap_index() = f2->gap_index ;
						vgaps.push_back(gap2) ;
						//double garea1 = compute_sector_gap_area(p0, p1, p12, exclude_radius) ;
						//double garea2 = compute_sector_gap_area(p0, p21, p2, exclude_radius) ;
						//area += garea1 ;
						//area += garea2 ;
					}
				}								
			} 
			else if(f1->has_gap && !f2->has_gap) {
				vec2 p1 = f1->dual ;
				Delaunay::Vertex_handle vopp = f1->vertex(f1->cw(f1->index(v))) ;
				vec2 p3 = to_geex(vopp->point()) ;
				vec2 dir = normalize(p3-p0) ;
				vec2 nor = vec2(-dir.y, dir.x) ;
//				double d03 = distance(p0, p3) ;
				vec2 mp = delaunay_->get_mid_point(v->point(), vopp->point()) ; //0.5*(p0+p3) ;
				double h = sqrt(SQR(v->radius)-distance2(p0, mp)) ; //SQR(d03*0.5)) ;
				vec2 p12 = mp - h*nor;

				VertexGap gap ;
				gap.set_disk(p0, v->radius) ;
				gap.push_back(p0+v->radius*normalize(p1-p0)) ;
				gap.push_back(p1) ;
				gap.push_back(p12) ;
				gap.gap_index() = f1->gap_index ;
				vgaps.push_back(gap) ;
				//double garea = compute_sector_gap_area(p0, p1, p12, radius) ;
				//area += garea ;
			}
			else if(!f1->has_gap && f2->has_gap) {
				vec2 p2 = f2->dual ;
				Delaunay::Vertex_handle vopp = f1->vertex(f1->cw(f1->index(v))) ;
				vec2 p3 = to_geex(vopp->point()) ;
				vec2 dir = normalize(p3-p0) ;
				vec2 nor = vec2(-dir.y, dir.x) ;
//				double d03 = distance(p0, p3) ;
				vec2 mp = delaunay_->get_mid_point(v->point(), vopp->point()) ; //0.5*(p0+p3) ;
				double h = sqrt(SQR(v->radius)-distance2(p0, mp)) ; //SQR(d03*0.5)) ;
				vec2 p21 = mp + h*nor ;

				VertexGap gap ;
				gap.set_disk(p0, v->radius) ;
				gap.push_back(p0+v->radius*normalize(p2-p0)) ;
				gap.push_back(p21) ;
				gap.push_back(p2) ;
				gap.gap_index() = f2->gap_index ;
				vgaps.push_back(gap) ;
				//double garea = compute_sector_gap_area(p0, p21, p2, exclude_radius) ;
				//area += garea ;
			}
			else {
				// no gap
			}
			++f1 ;
			++f2 ;
		} while(f1!=delaunay_->incident_faces(v)) ;
	}

	void GapCompute::compute_vertex_gaps_clipped(Delaunay::Vertex_handle v, std::vector<VertexGap>& vgaps) {
		gx_assert(v->dual_intersects_boundary) ;

		Polygon2* P = delaunay_->dual_convex_clip(v) ;
		double r2 = v->point().weight() ;
		vec2   p0 = to_geex(v->point()) ;

		for(int i=0; i<P->size(); ++i) {
			bool intersect0 = distance2((*P)[i].vertex[0], p0) > r2 ;
			bool intersect1 = distance2((*P)[i].vertex[1], p0) > r2 ;
			if(!intersect0 && !intersect1) {
				continue ;
			}
			else if(intersect0 && !intersect1) {
				vec2 pv0 = p0+v->radius*normalize((*P)[i].vertex[0]-p0) ;
				vec2 v01 = normalize((*P)[i].vertex[1] - (*P)[i].vertex[0]) ;
				Line<double> L((*P)[i].vertex[0], vec2(-v01.y, v01.x)) ;
				double h = L.side(p0) ;
				double t = sqrt(distance2((*P)[i].vertex[0], p0)-SQR(h)) - sqrt(r2-SQR(h)) ;
				vec2 p01 = (*P)[i].vertex[0] + t*v01 ;

				VertexGap gap ;
				gap.push_back((*P)[i].vertex[0]) ;
				gap.push_back(p01) ;
				gap.push_back(pv0) ;
				gap.set_disk(p0, v->radius) ;
//				gap.gap_index() = -1 ;
				vgaps.push_back(gap) ;
			}
			else if(!intersect0 && intersect1) {
				vec2 pv1 = p0+v->radius*normalize((*P)[i].vertex[1]-p0) ;
				vec2 v01 = normalize((*P)[i].vertex[1] - (*P)[i].vertex[0]) ;
				Line<double> L((*P)[i].vertex[0], vec2(-v01.y, v01.x)) ;
				double h = L.side(p0) ;
				double t = sqrt(distance2((*P)[i].vertex[1], p0)-SQR(h)) - sqrt(r2-SQR(h)) ;
				vec2 p10 = (*P)[i].vertex[1] - t*v01 ;

				VertexGap gap ;
				gap.push_back((*P)[i].vertex[1]) ;
				gap.push_back(pv1) ;
				gap.push_back(p10) ;
				gap.set_disk(p0, v->radius) ;
				vgaps.push_back(gap) ;
			}
			else {
				vec2 pv0 = p0+v->radius*normalize((*P)[i].vertex[0]-p0) ;
				vec2 pv1 = p0+v->radius*normalize((*P)[i].vertex[1]-p0) ;
				vec2 v01 = normalize((*P)[i].vertex[1] - (*P)[i].vertex[0]) ;
				Line<double> L((*P)[i].vertex[0], vec2(-v01.y, v01.x)) ;
				double h = L.side(p0) ;
				bool t0 = dot(p0-(*P)[i].vertex[0], v01)>0 ;
				bool t1 = dot(p0-(*P)[i].vertex[1], v01)>0 ;

				if(t0 && t1 || !t0 && !t1) {
					VertexGap gap ;
					gap.push_back(pv0) ;
					gap.push_back((*P)[i].vertex[0]) ;
					gap.push_back((*P)[i].vertex[1]) ;
					gap.push_back(pv1) ;
					gap.set_disk(p0, v->radius) ;
					vgaps.push_back(gap) ;
				}
				else{ // to gaps
					vec2   mp = p0+h*vec2(v01.y, -v01.x);
					double t = sqrt(r2-SQR(h)) ;
					VertexGap gap0, gap1 ;
					if(h>v->radius) {
						vec2 mp0 = p0 + v->radius*vec2(v01.y, -v01.x);

						gap0.push_back(mp0) ;
						gap0.push_back(pv0) ;
						gap0.push_back((*P)[i].vertex[0]) ;
						gap0.push_back(mp) ;
						gap0.set_disk(p0, v->radius) ;
						vgaps.push_back(gap0) ;

						gap1.push_back(pv1) ;
						gap1.push_back(mp0) ;
						gap1.push_back(mp) ;
						gap1.push_back((*P)[i].vertex[1]) ;
						gap1.set_disk(p0, v->radius) ;
						vgaps.push_back(gap1) ;
					} else {
						gap0.push_back((*P)[i].vertex[0]) ;
						gap0.push_back(mp-t*v01) ;
						gap0.push_back(pv0) ;
						gap0.set_disk(p0, v->radius) ;
						vgaps.push_back(gap0) ;

						gap1.push_back((*P)[i].vertex[1]) ;
						gap1.push_back(pv1) ;
						gap1.push_back(mp+t*v01) ;
						gap1.set_disk(p0, v->radius) ;
						vgaps.push_back(gap1) ;
					}
				}
			}
		}	
	}

} // end of namespace Geex
