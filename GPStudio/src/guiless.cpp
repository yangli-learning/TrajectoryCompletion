#include "guiless.h"
#include <QDebug>
#include "osg_utility.h"
#include "trajectory_completion.h"
#include "gps_trajectory.pb.h"
#include <boost/timer.hpp>
Guiless::Guiless(const std::string &trace_fname,
				 const std::string &bbox_fname){

	//	gflags::ParseCommandLineFlags(&argc,&argv,true);
	
	gps_trajectories_= new GPSTrajectories();
	
	valid_= false;
	qDebug() << "Start loading to file " << trace_fname.c_str() << ", " <<bbox_fname.c_str();
	if (gps_trajectories_->load(trace_fname)){
		qDebug()<<"finished loading traces...";
		if (importBBox(bbox_fname)){
			valid_ =true;
		}
	}
	qDebug() <<"Valid = " <<valid_;
}

Guiless::~Guiless(){
	//delete gps_trajectories_;
	// delete skeleton_;
}
void Guiless::reportTime(double dt, const char* label){
	qDebug() <<"[Time elapsed for " << label<<"]:" << dt<< endl;
}
int Guiless::exec(){
	// create voting image
	// initialize skeleton
	// extract skeleton
	// compute skeleton graph
	// cluster junction trajectory
	// complete trajectories
	boost::timer t; 
	double et0 = t.elapsed();
	skeleton_ = std::make_shared<L1Skeleton>(voting_image_,
											 heading_images_);

	skeleton_->initialize();

	int nBranches = skeleton_->nBranches();
	float hMax=2*skeleton_->getH0();
	int iter=0;
	while(skeleton_->getH0() < hMax){
		qDebug() <<"============== iter "<<iter<<"=================";
		qDebug() <<"nBranch = " << nBranches;
		bool converged;
		if (iter==0){
			converged = skeleton_->regularize(8);
		}else{
			converged = skeleton_->regularize();
		}
		skeleton_->identifyBranchPoints();
		int newNBranches = skeleton_->nBranches();
		if (converged || nBranches == newNBranches ){

			skeleton_->reinitialize(true);
		}else{
			skeleton_->reinitialize(false);
		}
		nBranches = newNBranches;
		iter++;
	}
	//	skeleton_->visualizeDensity();
	
	double et1 = t.elapsed();
	reportTime(et1-et0,"skeleton-construction");
	std::cout <<"Skeleton extraction finished in " << iter << " steps for "
			  << et1 <<" seconds";
	

	if (skeleton_->getBranches().size()>0){
		double et2 = t.elapsed();
		L1SkeletonGraph prototype(skeleton_->getBranches(),
                              skeleton_->getBranchGraph(),
                              skeleton_->getBridgePoints(),  voting_image_.rows
                              );
		
		double et3= t.elapsed();
		reportTime(et3-et0,"build-skeleton-graph");
		prototype.writeSkeletonImage("skeleton0.png",voting_image_.rows, voting_image_.cols);

		double et4= t.elapsed();
		
		junction_network_=std::make_shared<JunctionNetwork>(&prototype,skeleton_.get(),
                                          voting_image_.rows,
                                          voting_image_.cols );
		prototype.writeSkeletonImage("skeleton1.png",voting_image_.rows, voting_image_.cols);
		
		junction_network_->initJunctionGraph();
		junction_network_->validateTopologyByGPS(gps_trajectories_,min_corner_, resolution_,-1);
		
		double et5= t.elapsed();
		reportTime(et5-et4,"build-junction-network");
		prototype.writeSkeletonImage("skeleton2.png",voting_image_.rows, voting_image_.cols);
		
		exportClusters( "clusters.dat","clusters.index");


		SubTrajClusters  clusters = junction_network_->getTrajectoryClusters();
		GPSTrajectories::GpsTrajs gps_trajs = gps_trajectories_->getGpsTrajs();

		double et6= t.elapsed();
		TrajectoryCompletion trajCompl(&clusters, &gps_trajs,
									   cv::Point2f(min_corner_.x(), min_corner_.y()),
									   max_corner_.y()-min_corner_.y(), resolution_);
	
		trajCompl.completeAll2("completionResults.csv");
		double et7 = t.elapsed();
		reportTime(et7-et6,"trajectory-completion");
	   
		trajCompl.exportDenseTrajectories("dense_traj"); //export trajectories in txt format
		
	}
	return 0;
}

bool Guiless::importBBox(const std::string &fname){
	
	qDebug() <<"read bounding box from "<<fname.c_str() << endl;
	if (fname.length()== 0){
		qDebug() <<"No file selected " << endl;
		return false;
	}
	cv::FileStorage fs(fname, cv::FileStorage::READ);
	if (!fs.isOpened()){
		qDebug() <<"Error open " <<fname.c_str() << endl;
		return false;
	}
	//QWriteLocker lock(&read_write_lock_);
	//expired_ = true;

	cv::FileNodeIterator it1 = fs["min_corner"].begin();
	min_corner_[0] = *(it1++);
	min_corner_[1] = *it1;
	min_corner_[2] = 0;
	cv::FileNodeIterator it2 = fs["max_corner"].begin();
	max_corner_[0] = *(it2++);
	max_corner_[1] = *it2; max_corner_[2] = 0;
   
    fs["map_image"] >> voting_image_;
    resolution_ = fs["resolution"];
	cv::FileNode dn = fs["directions"];
	cv::FileNodeIterator it3;

    directions_.clear();
	for (it3 = dn.begin();it3!=dn.end();++it3){
		cv::Point2f f;
		(*it3) >> f;
        directions_.push_back(f);
	}
	cv::FileNode hn = fs["heading_images"];
    heading_images_.clear();
	for (cv::FileNodeIterator it4=hn.begin();it4!= hn.end();it4++){
		cv::Mat m;
		(*it4) >> m;
        heading_images_.push_back(m);
	}
    fs.release();

	qDebug()<<"resolution: " << resolution_;

	qDebug() <<"mincorner: " << min_corner_[0]<<"," << min_corner_[1]<< endl;
	qDebug() <<"maxcorner: " << max_corner_[0]<<"," << max_corner_[1]<< endl;
    //map_color_image_ = cv::Mat(voting_map_.map_image_.rows,
    //                           voting_map_.map_image_.cols,
	//							   CV_8UC3, cv::Scalar(255, 255, 255));

	double min, max;
    cv::minMaxLoc(voting_image_, &min, &max);
    qDebug() << "max_intensity: " << max << endl;
	qDebug() << "min_intensity: " << min << endl;

	
	return true;
}

void Guiless::exportClusters(//shared_ptr<JunctionNetwork> junction_network,
							 const std::string &fname,
							 const std::string &index_fname){
    if (junction_network_==nullptr){
        return;
    }
   
    std::ofstream file(fname.c_str(),std::ofstream::out);
    std::ofstream index_file(index_fname.c_str(), std::ofstream::out);
    index_file <<"TrajectoryId, ClusterId, StartPos, EndPos\n";
    SubTrajClusters clusters = junction_network_->getTrajectoryClusters();
    file << clusters.size() << std::endl;
	


    GPSTrajectories::GpsTrajs gps_trajs = gps_trajectories_->getGpsTrajs();
     int ct =0;
    for (auto it = clusters.begin();it!=clusters.end();++it){

		//int head  = it->first[0], tail = it->second[1];
		file << it->second.size() <<" " << it->first[0] << " "
			 << it->first[1] << std::endl;
		for (size_t j=0; j<it->second.size();++j){

			subtraj_descriptor subtraj = it->second[j];
			int trajId = std::get<0>(subtraj);
			int pStart = std::get<1>(subtraj);
			int pEnd = std::get<2>(subtraj);
			std::shared_ptr<GpsTraj> &traj = gps_trajs[trajId];
				
			file << pEnd- pStart;
			index_file << trajId<<"," << ct <<"," << pStart <<","
					   << pEnd //<<  ","
				//   << (subtraj.pEnd-subtraj.pStart)<<","<<traj->point_size()
					   << "\n";
			pEnd = std::min(pEnd,traj->point_size()-1);
			for (size_t k=pStart; k<pEnd; ++k){
				cv::Point2f p1 = junction_network_->toCvPoint2f(&(traj->point(k) ),
															  min_corner_,resolution_);
				cv::Point2f p((p1.x+0.5)*resolution_,
							  (voting_image_.rows-p1.y+0.5)*resolution_);
				file << std::setprecision(8)<< " "<<  p.x<<" " << p.y;
			}
			file <<std::endl;
		}
        
    }
    file.close();
    index_file.close();
   
}
