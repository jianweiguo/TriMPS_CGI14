#ifndef __GRID_H__
#define __GRID_H__

#include "geometry.h"

namespace Geex {

	enum GridCellType {CELL_INSIDE, CELL_OUTSIDE, CELL_BOUNDARY} ;

	class GridCell ;

	struct int_two{
		int x;
		int y;
	};
	/**
	  * Each grid cell indicate whether valid/invalid; records covered by which disk
	*/
	class GridCell {
	public:
		GridCell() : valid(true), type(CELL_INSIDE) 
		{ u = -1 ; v = -1 ; }

	public:
		bool valid ;
		GridCellType type ;
		std::vector<int> indices ; // stores sample indices that fully covers the cell
		unsigned int u ;
		unsigned int v ;
		//---------New added
		bool sampled;
		vec2 point;
		unsigned int vertexOrder;
		int_two clusterIndex;
		double radius;
		double min_dis;
		bool locked;

		bool visited;
		vec2 star_neighbor[25];
		int neighbor_number;
		vec2 edges[15];
		int edge_number;
        //---------End added
	} ;
	/**
	  * Uniform grid for Poisson disk sampling
	*/
	class UniformGrid {
	public:
		UniformGrid() {	}
		~UniformGrid() {
			if(cells_!=nil) {
				for(int i=0; i<res_; ++i) {
					delete [] cells_[i] ;
				}
				delete [] cells_ ;
			}
		}


		GridCell*& operator [] (const int i) { 
			gx_assert(i>=0 && i<res_) ;
			return cells_[i] ; 
		}

	protected:
		int        res_ ;
		GridCell** cells_ ;
	} ;


}

#endif
