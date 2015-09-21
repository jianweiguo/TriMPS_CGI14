/*
 * The code for the paper "Efficient triangulation of Poisson-disk sampled point sets".
 * If you have any probelms, don't hesitate to contact the authors.
 */

#include "delaunay.h"
#include <Geex/basics/file_system.h>
#include <stdlib.h>
#include <fstream>

namespace Geex {

    class DelaunayApp {
    public:
        DelaunayApp(){ 
			delaunay_instance_ = new Delaunay;
			non_convex_ = true ;
			min_radius_ = 0.01 ;
			delaunay_instance_->set_non_convex_mode(non_convex_);
			delaunay_instance_->radius_min() = min_radius_;
        }
		void load_boundary(const std::string& filename, double min_r){
			reset();

			min_radius_ = min_r ;
			delaunay_instance_->radius_min() = min_radius_;

			boundary_filename_ = Geex::FileSystem::get_project_root() + "/data/" + filename;
			if(boundary_filename_.length() > 0) {
                delaunay_instance_->load_boundary(boundary_filename_) ;
            }
		}
		
		void save_ps(const std::string& filename) {
			std::string fname = Geex::FileSystem::get_project_root() + "/data/" + filename ;
            delaunay_instance_->save_ps(fname) ;
		}

		void load(const std::string& filename) {
			std::string fname = filename ; 
			if(!Geex::FileSystem::is_file(fname)) {
				fname = Geex::FileSystem::get_project_root() + "/data/" + filename ;
			}
           delaunay_instance_->load(fname) ;
		}

        void reset() {
			delaunay_instance_->resetDelaunay();
        }
		void set_sample_mode(WeightScheme scheme){
			delaunay_instance_->sample_mode() = scheme;
		}

		void set_lfs_paras(double max_ratio, double lfs_w){
			delaunay_instance_->ratio_max_radius() = max_ratio;
			delaunay_instance_->lfs_weight() = lfs_w;
		}
		void generate_poisson_disk() {
			reset();
			delaunay_instance_->generate_poisson_disk() ;
			delaunay_instance_->fill_gaps() ;
		}

		void boundary_sample_cluster(){
			delaunay_instance_->boundary_sample_cluster();
		}

		void maximal_sampling(){
			delaunay_instance_->maximal_sampling();
		}

		void show_original_triangle(){
			delaunay_instance_->show_original_triangle() ;
		}

		void cell_cluster(){
			WeightScheme scheme = delaunay_instance_->sample_mode();
			if(scheme==WT_UNIFORM){
				delaunay_instance_->proxy_clustering();
			}
			if(scheme==WT_LFS){
				//delaunay_instance_->proxy_clustering_vary();
				delaunay_instance_->proxy_clustering();
			}
		}

		void triangle(){
			delaunay_instance_->generate_triangle();
			//delaunay_instance_->generate_triangle_vary();
		}

		void show_non_delanuary(){
			WeightScheme scheme=delaunay_instance_->sample_mode();
			if(scheme==WT_UNIFORM){
				delaunay_instance_->show_non_delanuary() ;
			}
			if(scheme==WT_LFS){
				delaunay_instance_->show_non_rt();
			}
		}

		void edge_flip(){
			delaunay_instance_->edge_flip() ;
		}

		void tri_poisson_disk() {
			boundary_sample_cluster() ;
			cell_cluster() ;
			triangle();
			show_non_delanuary();
			edge_flip();
			save_ps("out.ps");
		}

    private:
		Delaunay *delaunay_instance_;
        std::string boundary_filename_ ;
		double min_radius_ ;
        bool non_convex_ ;
		
    } ;
}

int main(int argc, char** argv) {
	Geex::initialize() ;
    if(argc != 4 && argc!=6) {
        std::cerr << "usage : " << std::endl ;
		std::cerr << "    Uniform: " << "boundary_fname sample_mode min_radius"<< std::endl ;
		std::cerr << "    Adaptive: " << "boundary_fname sample_mode min_radius max_radius_ratio lfs_weight"<< std::endl ;
        std::cerr << "example: " << std::endl ;
		std::cerr << "    Uniform: " << "circle.line U 0.025"<< std::endl ;
		std::cerr << "    Adaptive: " << "circle.line A 0.007 12.0 0.8"<< std::endl ;
		std::cerr << std::endl ;
		std::cerr << "Note that the boundary_fname should be in the ./data folder"<< std::endl ;
        return -1 ;
    }
	std::string boundary_fname = argv[1];
	char* sample_mode = argv[2]; 
	double min_r = atof(argv[3]);
	
	Geex::DelaunayApp *delaunayApp = new Geex::DelaunayApp;
    std::cerr << "============= Load boundary ==========" << std::endl ;
	delaunayApp->load_boundary(boundary_fname, min_r);
	std::cerr << std::endl ;

    std::cerr << "============= Poisson-disk sampling  ==========" << std::endl ;
	if(sample_mode[0] == 'U'){
		delaunayApp->set_sample_mode(Geex::WT_UNIFORM);
		delaunayApp->maximal_sampling();
	}
	else if(sample_mode[0] == 'A'){
		delaunayApp->set_sample_mode(Geex::WT_LFS);
		double max_ratio = atof(argv[4]);
		double lfs_weight = atof(argv[5]);
		delaunayApp->set_lfs_paras(max_ratio, lfs_weight);
		delaunayApp->generate_poisson_disk();
	}
	std::cerr << std::endl ;

	std::cerr << "============= Triangulate Poisson-disk sampling  ==========" << std::endl ;
	delaunayApp->tri_poisson_disk();
    return 0 ;
}