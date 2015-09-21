#ifndef __GAP_FUNCS_H__
#define __GAP_FUNCS_H__

#include "gap_primitive.h"
#include "delaunay.h"

namespace Geex {

	class GapCompute {
	public:
		GapCompute(Delaunay* dt) : delaunay_(dt) {
		}
	public:
		//
		void compute_face_gap(Delaunay::Face_handle f, std::vector<FaceGap>& fgaps) ;
		//
		void compute_face_vertex_gap(Delaunay::Face_handle f, std::vector<FaceGap>& fgaps) ;
		//
		void compute_vertex_gaps(Delaunay::Vertex_handle v, std::vector<VertexGap>& vgaps) ;

	protected:
		void compute_vertex_gaps_inner(Delaunay::Vertex_handle v, std::vector<VertexGap>& vgaps) ;
		void compute_vertex_gaps_clipped(Delaunay::Vertex_handle v, std::vector<VertexGap>& vgaps) ;

	private:
		Delaunay* delaunay_ ;

	} ;
}

#endif
