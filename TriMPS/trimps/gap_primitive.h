#ifndef __GAP_PRIMITIVES_H__
#define __GAP_PRIMITIVES_H__

#include "geometry.h"

namespace Geex {


// The basic representation of a  gap element related to a vertex
// A gap unit is the sector of the Voronoi cell of a vertex 
// excluded by the circle. The polygon can only be a triangle or a quad
// Each gap unit is represented by a set of vertex and the incident disk.

   class VertexGap : public std::vector<vec2> {
   public:
	   VertexGap():gap_index_(-1) {}

	   vec2 random_point() {
		   gx_assert(size()==3 || size()==4) ;
		   if(size()==3) {
			   return random_point_tri((*this)[0], (*this)[1], (*this)[2]) ;
		   }
		   else {
			   return random_point_quad((*this)[0], (*this)[1], (*this)[2], (*this)[3]) ;
		   }
	   }

	   void set_disk(vec2& c, double r) {
		   disk_.first = c ;
		   disk_.second = r ;
	   }

	   double area() {
		   if(size()==3) {
			   return triangle_area((*this)[0], (*this)[1], (*this)[2]) ;
		   }
		   else {
			   return triangle_area((*this)[0], (*this)[1], (*this)[2]) 
				     + triangle_area((*this)[0], (*this)[2], (*this)[3]) ;
		   }
	   }

	   vec2 center() {
		   vec2 c(0 ,0) ;
		   for(unsigned int i=0; i<this->size(); ++i) {
			   c += (*this)[i] ;
		   }
		   return 1.0/(this->size())*c ;
	   }

	   std::pair<vec2, double>& disk() { return disk_ ; }
	   int& gap_index() { return gap_index_ ; }

   protected:
	   std::pair<vec2, double> disk_ ;
	   int gap_index_ ;
   } ;


   // gap unit related to a cell in the triangulation
   class FaceGap : public std::vector<vec2> {
   public: 
	   FaceGap() { } ;

	   void add_disk(vec2& c, double w) {
		   disks_.push_back(c) ;
		   weights_.push_back(w) ;
	   } 

	   void clear_all() {
		   clear() ;
		   disks_.clear() ;
		   weights_.clear() ;
	   }

	   std::vector<vec2>&   disks()   { return disks_ ; }
	   std::vector<double>& weights() { return weights_ ; }

	   double area() {
		   double area = 0 ;
		   for(unsigned int i=0; i<size()-2; ++i) {
			   area += triangle_area((*this)[0], (*this)[i+1], (*this)[i+2]) ;
		   }
		   return area ;
	   }

	   vec2 center() {
		   vec2   g(0, 0) ;
		   double w = 0 ;
		   for(int i=0; i<disks_.size(); ++i) {
			   g = g + weights_[i]*disks_[i] ;
			   w += weights_[i] ;
		   }

		   g = 1.0/w * g ;
		   return g ;
	   }

	   bool random_point(vec2& p) {
//		   vec2 p(0, 0) ;
		   bool suc = false ; 
		   double totalw = 0 ;
		   int ntry = 0 ;


		   while(!suc) {
			   suc = true ;
			   
			   switch(size()) {
			   case 3: {
				   p = random_point_tri((*this)[0], (*this)[1], (*this)[2]) ; 
				   break ;
			   }
			   case 4: {
				   p = random_point_quad((*this)[0], (*this)[1], (*this)[2], (*this)[3]) ;
				   break ;
			   }
			   default:
				   gx_assert(size()>4) ;
				   p = random_point_poly(*this) ; 
				   break ;
			   //else {
				  // p = vec2(0, 0) ;
				  // for(unsigned int i=0; i<size()-1; ++i) {
					 //  //double w = Numeric::random_float64()*(1.0-totalw) ;
					 //  double w = Random::UniformRandom()*(1.0-totalw) ;
					 //  totalw += w ;
					 //  p += w*(*this)[i] ;
				  // }
				  // p += (1-totalw)*(*this)[size()-1] ;
			   //}
			   }

			   for(int j=0; j<disks_.size(); ++j) {
				   if(distance2(disks_[j], p) < weights_[j]) {
					   suc = false ;
					   break ;
				   }
			   }
			   ntry ++ ;
			   if(ntry > 50) {
				   std::cerr << "can't generate point in gap" << std::endl ;
				   for(unsigned int i=0; i<size(); ++i) {
					   std::cerr << (*this)[i] << std::endl ;
				   }
				   std::cerr << "disks " << disks_.size() << std::endl ;
				   for(unsigned int i=0; i<disks_.size(); ++i) {
					   std::cerr << disks_[i] << std::endl ;
				   }
				   return false ;
//				   break ;
			   }
		   }
//		   return p ;
		   return true ;
	   }

   protected:
	   std::vector<vec2>   disks_ ;
	   std::vector<double> weights_ ;
   } ;
}

#endif
