#ifndef __POLYGONS__
#define __POLYGONS__


#include <Geex/mathematics/glsl_linear.h>
#include "geometry.h"
#include <set>
#include <algorithm>
#include <iterator>

namespace Geex {

    class Convex : public std::vector< Line<real> > {
    public:
        inline bool contains(const vec2& p) {
            for(unsigned int i=0; i<size(); i++) {
                if((*this)[i].side(p) < 0) {
                    return false ;
                }
            }
            return true ;
        }
    } ;


    template <class T> inline void sets_intersect(
        const std::set<T>& S1, const std::set<T>& S2, std::set<T>& I
    ) {
        I.clear() ;
        std::set_intersection(S1.begin(), S1.end(), S2.begin(), S2.end(), std::inserter(I, I.begin())) ;
    }

    class PolygonVertex : public vec2 {
    public:
        PolygonVertex() : on_boundary(false), gap_index(-1) { }
        PolygonVertex(real x, real y) : vec2(x,y), on_boundary(false), gap_index(-1) { }
        PolygonVertex(const vec2& V) : vec2(V), on_boundary(false), gap_index(-1) { }
        
        std::set<int> bisectors ;
        std::set<int> boundary_edges ;
        bool on_boundary ;
		int  gap_index ; // used for classify gaps incident to boundary.
    } ;


    class Polygon : public std::vector<PolygonVertex> {
    public:
        void load(const std::string& filename) ;
        void normalize() ;
        bool contains(const vec2& p) ;
        Line<real> edge_line(unsigned int i) const ;
        void convex_clip(Polygon& to, const Line<real>& L, int E) ;
    } ;

    class PolygonEdge {
    public:
        PolygonEdge() { }
        PolygonEdge(const PolygonVertex& v1, const PolygonVertex& v2) {
            vertex[0] = v1 ;
            vertex[1] = v2 ;
        }
        Line<real> line() const {
            vec2 V = vertex[1] - vertex[0] ;
            return Line<real>(vertex[0], vec2(-V.y, V.x)) ;
        }
		bool on_boundary() {
			std::set<int> e ;
			sets_intersect(vertex[0].boundary_edges, vertex[1].boundary_edges, e) ;
			return !e.empty() ;
		} 
		int boundary_index() {
			std::set<int> e ;
			sets_intersect(vertex[0].boundary_edges, vertex[1].boundary_edges, e) ;
			if(e.empty())
				return -1 ;
			else {
				gx_assert(e.size()==1) ;
				return *(e.begin()) ;
			}
		}
		double length() {
			return (vertex[1]-vertex[0]).length() ;
		}
		vec2 mid_point() {
			return 0.5*(vertex[1]+vertex[0]) ;
		}
	public:
        PolygonVertex vertex[2] ;
    } ;

    class Polygon2 : public std::vector<PolygonEdge> {
    public:
        void load(const std::string& filename) ;
        void normalize() ;
        bool contains(const vec2& p) ;
        void convex_clip(Polygon2& to, const Line<real>& L, int E, bool close = true) ;
		double area() ;
		//New Add
		void getBoundingBox();
		double bbox[4];//the order is x_max, x_min,y_max,y_min
		bool onBounary(const vec2& p);
		bool segement_intersect(vec2 p,vec2 q);
		double cross(vec2 a,vec2 b,vec2 c);
		double max(double a,double b);
		double min(double a,double b);
		//End Add
    } ;
}

#endif
