#include <cassert>
#include <QDebug>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>

#include <opencv2/core/core_c.h>
#include <opencv2/core/version.hpp>
#if (CV_VERSION_EPOCH  > 2)
#include <opencv2/imgcodecs.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#else
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
//#include "l1skeleton_graph.h"
//#include <opencv2/highgui/highgui_c.h>
#endif
#include "common.h"
#include "parameters.h"
#include "osg_utility.h"
#include "main_window.h"
#include "scene_widget.h"
//#include "gps_trajectories.h"
#include "gps_trajectory.pb.h"

#include "bbox_2d.h"
#include "trajectory_completion.h"

BBox2D::BBox2D(GPSTrajectories* scene_image) :
    gps_trajectories_(scene_image), corner_number_(0), built_(false),
    show_map_image_(true),
    voting_map_(this){

	resolution_= Parameters::getInstance()->map_image_resolution;
	skeleton_ = 0;
   // skeleton_heading_=0;
	skeleton_converged_ = false;
	skeleton_nodes_= 0 ;
    branch_nodes_ = 0;
	cluster_nodes_=0;
	heading_nodes_=0;
    junction_network_ = 0;
}

BBox2D::~BBox2D(void) {
	delete skeleton_;
    if (junction_network_ != nullptr){
        delete junction_network_;
	}
}

void BBox2D::setBuilt(bool built) {
	qDebug()<<"note: set built is called!\n";
  if (built_ == built)
    return;

  QWriteLocker lock(&read_write_lock_);
  expired_ = true;

  built_ = built;

  if (built_) {
    float min_x = std::min(min_corner_.x(), max_corner_.x());
    float min_y = std::min(min_corner_.y(), max_corner_.y());
    float max_x = std::max(min_corner_.x(), max_corner_.x());
    float max_y = std::max(min_corner_.y(), max_corner_.y());
    min_corner_ = osg::Vec3(min_x, min_y, 0.0f);
    max_corner_ = osg::Vec3(max_x, max_y, 0.0f);

    lock.unlock();
    voting_map_.init();
    //voting_map_.voteMapImage();
  }
  return;
}
bool BBox2D::isDegenerated(void) {
  if (corner_number_ != 2)
    return false;
  return getWidth() * getHeight() > 400 * resolution_ * resolution_;
}

void BBox2D::clearCorners(void) {
  if (corner_number_ == 0)
    return;

  QWriteLocker lock(&read_write_lock_);
  expired_ = true;

  corner_number_ = 0;

  return;
}

void BBox2D::addCorner(const osg::Vec3& corner) {
  QWriteLocker lock(&read_write_lock_);
  expired_ = true;
  if (corner_number_ == 0) {
    min_corner_ = corner;
    corner_number_ = 1;
  } else if (corner_number_ == 1 || corner_number_ == 2) {
    max_corner_ = corner;
    corner_number_ = 2;
  } else {
    assert(false);
  }

  return;
}

void BBox2D::pickEvent(PickMode pick_mode, float x, float y) {
  switch (pick_mode) {
  case PickMode::CTRL:
    break;
  case PickMode::SHIFT:
    MainWindow::getInstance()->getSceneWidget()->removeBBox(this);
    break;
  case PickMode::CTRLSHIFT:
	  {
		  int idx = getPickedSkeletonPoint(x,y);

		  if (idx >= 0 ){
              std::vector<SkeletonPoint> pts = skeleton_->getSkeleton();
			  cv::Point2f coor = pts[idx].pos;
              qDebug()<<"=======\nlocation: ("<<coor.x <<"," << coor.y<<")";
			  //visualize direction 
              emit scene_widget_->pointSelected(idx,coor.x,coor.y);
              highlightNeighbors(idx);
              qDebug() <<"index " << idx;
              qDebug() << "id: " << pts[idx].id;
              switch (pts[idx].type){
                case SkeletonPointType::BRANCH:
                    qDebug() << "type: branch"; break;
                case SkeletonPointType::NOISE:
                    qDebug() << "type: noise"; break;
                case SkeletonPointType::BRIDGE:
                    qDebug() << "type: bridge"; break;
                default:
                    qDebug() << "type: non-branch point"; break;
              }

              qDebug() << "sigma: " << skeleton_->getSkeletonLinearity(idx);
              qDebug() <<"===============";
		  }
		  break;
	  }
  case PickMode::ALTSHIFT:
	  {
		  int idx = getPickedSkeletonPoint(x,y);
		  if (idx >= 0 ){
			  qDebug() <<"Trace a branch from " << idx<< "th skeleton point ";
			  // std::vector<SkeletonPoint> pts = skeleton_->getSkeleton();
			  traceBranchFromPoint(idx);
		  }
		  break;
	  }
  default:
	  break;
  }
}
void BBox2D::traceBranchFromPoint(int skeleton_id){
	// identify if selected branch is a candidate point
	// if so, find the branch traced from this point 
	// draw the trace result
    if (L1Skeleton* sk = dynamic_cast<L1Skeleton*>(skeleton_)){
        auto branchNoise = sk->getBranchFromSeed(skeleton_id);
        auto branch = branchNoise.first;
        auto noise = branchNoise.second;
        std::vector<cv::Point2f> points , noisePoints;
        qDebug() <<"Traced branch of length " << branch.size()<< " with "
				 << noisePoints.size() << " noise pts";
        for (size_t i=0;i<branch.size();++i){
            points.push_back(branch[i].pos);
        }
        for (size_t i=0;i<noise.size();++i){
            noisePoints.push_back(noise[i].pos);
        }
        if (branch_nodes_ == NULL){
            branch_nodes_ = new osg::Geode();
            content_root_->addChild(branch_nodes_);
        }
        drawPolyLine(branch_nodes_,points,osg::Vec4(1.0f, 1.0f, 0.8f, 0.3f));
        drawPoints(branch_nodes_,points,osg::Vec4(1.0f, 1.0f, 0.8f, 0.3f));
        drawPoints(branch_nodes_, noisePoints,osg::Vec4(0.f,0.f,1.f,1.f));
    }
}

int BBox2D::getPickedSkeletonPoint(float x, float y){
    if (skeleton_ && skeleton_->size() == 0) {
		qDebug() << "No unidentified skeleton points. No point selected.";
		return -1;
    }
    cv::Point2f mapCoor = toMapCoordinates(x,y);
    std::vector<IndexVal<float> > result =
		skeleton_->findNearestKSkeletonPoints(mapCoor.x,mapCoor.y,1);
    /*
	  if (!Parameters::getInstance()->skeleton_use_heading){
        result = dynamic_cast<L1SkeletonWithHeading*>(skeleton_)->findNearestKSkeletonPoints(mapCoor.x,mapCoor.y,1);
    }else{
        result = dynamic_cast<L1Skeleton*>(skeleton_)->findNearestKSkeletonPoints(mapCoor.x,mapCoor.y,1);
    }*/
	if (result.size()>0){
		IndexVal<float> nearest =result.front(); 

		if (nearest.val < resolution_){
			return nearest.index;
		}
	}
    return -1;
}
void BBox2D::toggleRenderMapImage(bool show_map_image) {
	if (show_map_image == show_map_image_){
		return;
	}

  QWriteLocker locker(&read_write_lock_);
  expired_ = true;

  show_map_image_ = !show_map_image_;

  return;
}



void BBox2D::renderMapImage(void) {
    if (!show_map_image_ || corner_number_ != 2
        || voting_map_.map_image_.data == nullptr){

        return;
    }
    qDebug()<<"compute map color image with "<< corner_number_<<" corners";
    voting_map_.updateMapColorImage(VotingMapRenderer::GAMMA_COMPRESSION);

    int width = voting_map_.map_image_.cols;
    int height = voting_map_.map_image_.rows;
    osg::ref_ptr<osg::Image> image = new osg::Image;
    image->setImage(width, height, 3, GL_RGB, GL_BGR, GL_UNSIGNED_BYTE,
                    (unsigned char*) (voting_map_.map_color_image_.ptr()),
                    osg::Image::AllocationMode::NO_DELETE);

    osg::ref_ptr<osg::Geometry> quad =
		osg::createTexturedQuadGeometry(osg::Vec3(min_corner_.x(),
												  min_corner_.y(), -0.5f),
										osg::Vec3(getWidth(), 0.0f, 0.0f),
										osg::Vec3(0.0f, getHeight(), 0.0f),
										0.0f, 1.0f, 1.0f, 0.0f);

    osg::ref_ptr<osg::Texture2D> texture = new osg::Texture2D;
    texture->setResizeNonPowerOfTwoHint(false);
    texture->setImage(image);
    quad->getOrCreateStateSet()->setTextureAttributeAndModes(0, texture.get());
    quad->getOrCreateStateSet()->setMode(GL_DEPTH_TEST, osg::StateAttribute::ON);

    osg::ref_ptr<osg::Geode> texture_holder = new osg::Geode();
    texture_holder->addDrawable(quad);
    content_root_->addChild(texture_holder);
}

void BBox2D::updateImpl(void) {
  if (corner_number_ == 0) {
    return;
  }
  osg::Vec4 color(1.0f, (built_ ? (1.0f) : (0.0f)), 0.0f, 1.0f);
  float corner_radius = 1.5f * resolution_;
  float edge_radius = 1.0f * resolution_;
  if (corner_number_ == 1) {
    content_root_->addChild(OSGUtility::drawSphere(min_corner_, corner_radius, color));
  } else if (corner_number_ == 2) {
    osg::Vec3 corner_1(min_corner_);
    osg::Vec3 corner_4(max_corner_);
    corner_1.x() = max_corner_.x();
    corner_4.x() = min_corner_.x();

    content_root_->addChild(OSGUtility::drawSphere(min_corner_, corner_radius, color));
    content_root_->addChild(OSGUtility::drawSphere(corner_1, corner_radius, color));
    content_root_->addChild(OSGUtility::drawSphere(max_corner_, corner_radius, color));
    content_root_->addChild(OSGUtility::drawSphere(corner_4, corner_radius, color));

    content_root_->addChild(OSGUtility::drawCylinder(min_corner_, corner_1, edge_radius, color));
    content_root_->addChild(OSGUtility::drawCylinder(corner_1, max_corner_, edge_radius, color));
    content_root_->addChild(OSGUtility::drawCylinder(max_corner_, corner_4, edge_radius, color));
    content_root_->addChild(OSGUtility::drawCylinder(corner_4, min_corner_, edge_radius, color));

  }

  renderMapImage();
  return;
}

cv::Point BBox2D::toCvPoint(const TrajPoint* traj_point,
                            bool convert) {
    cv::Point2f pf = toCvPoint2f(traj_point,convert);
    return cv::Point(round(pf.x), round(pf.y));
}

cv::Point2f BBox2D::toCvPoint2f(const TrajPoint* traj_point,
                                bool convert) {
  float x = traj_point->x();
  float y = traj_point->y();
  return (convert)?toMapCoordinates(x,y):cv::Point2f(x,y);
}

#ifdef _MSC_VER
static float round(float x) {
  return x >= 0.0f ? floorf(x + 0.5f) : ceilf(x - 0.5f);
}
#endif



void BBox2D::exportBBox(const std::string &fname){
		qDebug()<<"write bounding box info to "<<fname.c_str() << endl;
	cv::FileStorage fs(fname,cv::FileStorage::WRITE);
	fs <<"min_corner"<< "[:" << min_corner_[0]<<min_corner_[1]
	   <<"]";
	fs <<"max_corner"<< "[:"<< max_corner_[0] << max_corner_[1]
	   <<"]";
    fs <<"map_image"<< voting_map_.map_image_;
	fs <<"resolution" << resolution_;
	fs<<"directions"<< "[";
    for (size_t d=0;d<voting_map_.directions_.size();++d){
        fs << voting_map_.directions_[d];
	}
	fs<<"]";
	fs <<"heading_images"<< "[";
    for (size_t d=0;d<voting_map_.heading_images_.size();++d){
        fs <<voting_map_.heading_images_[d];
	}
	fs <<"]";
	fs.release();
	return ;

}

bool BBox2D::importBBox(const std::string& fname){

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
    voting_map_.init();
    fs["map_image"] >> voting_map_.map_image_;
    resolution_ = fs["resolution"];
	cv::FileNode dn = fs["directions"];
	cv::FileNodeIterator it3;

    voting_map_.directions_.clear();
	for (it3 = dn.begin();it3!=dn.end();++it3){
		cv::Point2f f;
		(*it3) >> f;
        voting_map_.directions_.push_back(f);
	}
	cv::FileNode hn = fs["heading_images"];
    voting_map_.heading_images_.clear();
	for (cv::FileNodeIterator it4=hn.begin();it4!= hn.end();it4++){
		cv::Mat m;
		(*it4) >> m;
        voting_map_.heading_images_.push_back(m);
	}
    fs.release();
	corner_number_=2;
	built_ = true;

	qDebug() <<"mincorner: " << min_corner_[0]<<"," << min_corner_[1]<< endl;
	qDebug() <<"maxcorner: " << max_corner_[0]<<"," << max_corner_[1]<< endl;
    voting_map_.map_color_image_ = cv::Mat(voting_map_.map_image_.rows,
                               voting_map_.map_image_.cols,
							   CV_8UC3, cv::Scalar(255, 255, 255));
   voting_map_.setValid(true);
   // lock.unlock();
	updateImpl();

	return true;
}
//>>>>>>>>>>>>> EXTRACT SKELETON >>>>>>>>>>>>>>>>>>>>>>>>>>>>> 
void BBox2D::initializeSkeleton(void){
    if (!voting_map_.isValid()){
        qDebug() <<"Error: voting map has not been initialized!";
        return;
    } 
    skeleton_ = new L1Skeleton(voting_map_.map_image_,
                               voting_map_.heading_images_);
    skeleton_->initialize();
	double min, max;
    cv::minMaxLoc(voting_map_.map_image_, &min, &max);
    qDebug() << "max_intensity: " << max << endl;
	qDebug() << "min_intensity: " << min << endl;

	std::vector<SkeletonPoint> S = skeleton_->getSamples();


	osg::Vec4 color(0.0f, 0.0f, 1.0f, 0.3f);
	//	float s = Parameters::getInstance()->skeleton_h0;
	qDebug() <<"Size of supporting neighborhood h0: " << skeleton_->getH0(); 


	double corner_radius = 1.0f * resolution_;
	osg::ref_ptr<osg::Vec3Array> sample_locs = new osg::Vec3Array;
	osg::ref_ptr<osg::Vec4Array> sample_colors = new osg::Vec4Array;

    for (size_t i=0; i<S.size(); ++i) {
        sample_locs->push_back(toVec3(S[i].pos));
		sample_colors->push_back(color);
	}

	skeleton_nodes_ =OSGUtility::drawSpheres(sample_locs, 
											 corner_radius,
											 sample_colors); 
	content_root_->addChild(skeleton_nodes_);
	if (EXPORT_SAMPLES){
		std::ostringstream oss;
		oss <<"skeletonPoints_iter"<< std::setfill('0') << std::setw(2)
			<<skeleton_->getIter()<<".csv"; 
		exportSamples(oss.str().c_str());
		
	}
}
void BBox2D::incrementH0(void){
	if (! skeleton_ ){
		qDebug() << "skeleton has not been initialized!\n";
		return;
	}
	skeleton_->reinitialize(true);

}
void BBox2D::extractSkeletonStep(void){	
	
    if (! skeleton_ ){
		qDebug() << "skeleton has not been initialized!\n";
		return;
	}
	//float s = Parameters::getInstance()->skeleton_scale;
	skeleton_converged_ = skeleton_->iterate_step(); //skeleton_->getH0()*s); 

    if (skeleton_converged_ ){
		skeleton_->identifyBranchPoints();
		qDebug() << "Found " << skeleton_->nBranches() << " branches..";
		drawBranches();
		skeleton_->reinitialize();

	}
	highlightNeighbors(-1); //draw skeleton without highlight

	if (EXPORT_SAMPLES){
		std::ostringstream oss;
		oss <<"skeletonPoints_iter"<< std::setfill('0') << std::setw(2)
			<<skeleton_->getIter()<<".csv"; 
		exportSamples(oss.str().c_str());	
	}
	
	return;
}

void BBox2D::extractSkeleton(void){	
	
	if (! skeleton_ ){
		qDebug() << "skeleton has not been initialized!\n";
		return;
	}
    if ( graph_nodes_!=nullptr){
        content_root_->removeChild(graph_nodes_);
    }
	//float s = Parameters::getInstance()->skeleton_scale;
	skeleton_converged_ = skeleton_->regularize(); //skeleton_->getH0()*s); 

	//	if (skeleton_converged_){
	skeleton_->identifyBranchPoints();
	qDebug() << "Found " << skeleton_->nBranches() << " branches..";
	drawBranches();
	qDebug() << "reinitialize skeleton...";
	skeleton_->reinitialize();
	highlightNeighbors(-1); //draw skeleton without highlight
	//}
	qDebug()<<"finish extracting skeleton..";
	return;

}

void BBox2D::drawSkeletonGraph(L1SkeletonGraph &graph){

    content_root_->removeChild(skeleton_nodes_);
    content_root_->removeChild(branch_nodes_);

    graph_nodes_ = new osg::Geode();
    //typedef graph_traits<SGraph>::out_edge_iterator edge_iterator;

    //std::map<EdgeMapKey, std::vector<int> >::const_iterator it;

    std::vector<std::vector<cv::Point2f> > edgeSegments = graph.getEdgeSegments();
    // std::vector<(edgeIndexSegments.size());

    //draw graph edges
    //std::vector<cv::Point2f> allPts;
    for (size_t ct = 0; ct < edgeSegments.size(); ++ct){

        int colorIdx = (ct) % Common::palette_size;

        osg::Vec4 color( (Common::palette[colorIdx][0])/255,
                         (Common::palette[colorIdx][1])/255,
                         (Common::palette[colorIdx][2])/255, 1.f);

       // allPts.insert(allPts.end(), edgeSegments[ct].begin(), edgeSegments[ct].end());
         drawPolyLine(graph_nodes_, edgeSegments[ct] ,color,0.4*resolution_ );

    }

    //draw vertices
     std::vector<cv::Point2f> vertexPoints = graph.getVertexPositions();

    drawPoints(graph_nodes_, vertexPoints, osg::Vec4(0.35,0.35,0.35,0.5f));

   // drawPoints(graph_nodes_, allPts, osg::Vec4(1.f,1.f,1.f,0.6f),0.6*resolution_);
    content_root_->addChild(graph_nodes_);

    return;
}

void BBox2D::drawSkeletonGraph(CompactSkeletonGraph &graph ){
	// remove branches and sample points from view
	//draw skeleton edges as 

	content_root_->removeChild(skeleton_nodes_);
	content_root_->removeChild(branch_nodes_);

	graph_nodes_ = new osg::Geode();
	
	std::map<EdgeMapKey, std::vector<int> >::const_iterator it;
	int ct=0;

	//draw graph edges
	for (it = graph.edge_map_.begin();it!= graph.edge_map_.end();
		 ++it){
		if (it->first.first < it->first.second)
			continue;

        int colorIdx = (ct++) % Common::palette_size;
		
		osg::Vec4 color( (Common::palette[colorIdx][0])/255, 
						 (Common::palette[colorIdx][1])/255,
						 (Common::palette[colorIdx][2])/255, 1.f);

		std::vector<cv::Point2f> edgePoints;
		for (size_t i=0;i<it->second.size();++i){
			edgePoints.push_back( graph.points_[it->second[i]]);
			/*
			qDebug()<<"draw edge " << ct-1 <<" [" << it->first.first<<"," 
					<< it->first.second
					<<"] " << graph.points_[it->second[i]].x <<"," 
					<< graph.points_[it->second[i]].y;*/
		}
		drawPolyLine(graph_nodes_,edgePoints,color );

	}
	//draw vertices
	std::vector<cv::Point2f> vertexPoints;
	for (size_t i=0;i<graph.vertex_map_.size();++i){
	   if (graph.vertices[i].outSize >=1){
	
		   vertexPoints.push_back(graph.points_[graph.vertex_map_[i]]);		
		   //  qDebug()<<"Keep: " ;
		   //graph.printVertex(i);
	   }
	}
    drawPoints(graph_nodes_, vertexPoints,osg::Vec4(0.35,0.35,0.35,0.5f) );
	content_root_->addChild(graph_nodes_);

	return;

}
void BBox2D::computeSkeletonGraph(float minFlow){
    if ( skeleton_ == nullptr || skeleton_->getBranches().size()==0){
        qDebug() << "Warning: can not construct skeleton graph!\n";
        return ;
    }
    qDebug()<<"Computing skeleton graph...";
    L1SkeletonGraph prototype(skeleton_->getBranches(),
                              skeleton_->getBranchGraph(),
                              skeleton_->getBridgePoints(),  voting_map_.map_image_.rows
                              );

    junction_network_=new JunctionNetwork(&prototype,skeleton_,
                                          voting_map_.map_image_.rows,
                                          voting_map_.map_image_.cols );
     drawSkeletonGraph(prototype);
     if (graph_nodes_!=nullptr){
         std::vector<cv::Point2f> intersectPoints = junction_network_->getIntersectionPoints();
         drawPoints(graph_nodes_, intersectPoints, osg::Vec4(1.f,0.f,0.f,1.0f), 1.4f*resolution_);
     }
     junction_network_->initJunctionGraph();
     junction_network_->validateTopologyByGPS(gps_trajectories_,min_corner_, resolution_, minFlow);
	 prototype.writeSkeletonImage("skeleton.png",voting_map_.map_image_.rows,
								  voting_map_.map_image_.cols);
}
/*
void BBox2D::computeSkeletonGraph(){
	if ( skeleton_ == nullptr){
		qDebug() << "skeleton has not been initialized!\n";
		return ;
	}
	//	qDebug() <<"Recentering branches ... ";
	//	skeleton_->recenterBranches();
	if (skeleton_->getBranches().size()==0){
		qDebug() << "branches have not been detected...";
		return;
	}
	SkeletonGraph prototype(skeleton_->getBranches(),
							skeleton_->getBranchGraph(),
							skeleton_->getBridgePoints(),
							skeleton_->size());
	qDebug()<< "skeleton graph has "<< prototype.nVertices ;
	CompactSkeletonGraph compactGraph(prototype.nCompactVertices);
	prototype.buildCompactGraph(&compactGraph);
	qDebug()<<"finished computing compact graph ";
	drawSkeletonGraph(compactGraph);	
	process_graph_ = new ProcessGraph(&compactGraph, 
                                      voting_map_.map_image_.rows,
                                      voting_map_.map_image_.cols);

    drawSkeletonGraphJunction();
	process_graph_->validateTopologyByGPS(gps_trajectories_,min_corner_, resolution_);
	std::string dotstr=	process_graph_->toDotString();
	qDebug() << dotstr.c_str();
	process_graph_->exportGraph("junctionGraph.gv");
	qDebug() <<"num of clusters (after): " << process_graph_->getNumClusters();
	//	delete compactGraph;

}
*/


void BBox2D::exportClusters(std::vector<bool> selected,
                            const std::string& fname, const std::string &index_fname){
    if (cluster_nodes_==nullptr){
        return;
    }   
    int numSelected=0;
    for (size_t i=0; i<selected.size();++i){
        if (selected[i]) numSelected++;
    }
    std::ofstream file(fname.c_str(),std::ofstream::out);
    std::ofstream index_file(index_fname.c_str(), std::ofstream::out);
    index_file <<"TrajectoryId, ClusterId, StartPos, EndPos\n";
    file << numSelected << std::endl;
    SubTrajClusters clusters = junction_network_->getTrajectoryClusters();

    GPSTrajectories::GpsTrajs gps_trajs = gps_trajectories_->getGpsTrajs();
     int ct =0;
    for (auto it = clusters.begin();it!=clusters.end();++it){
        if (selected[ct]){
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
                    cv::Point2f p1 = toCvPoint2f(&(traj->point(k) ),true);
					cv::Point2f p((p1.x+0.5)*resolution_,
                                  (voting_map_.map_image_.rows-p1.y+0.5)*resolution_);
                    file << std::setprecision(8)<< " "<<  p.x<<" " << p.y;
                }
                file <<std::endl;
            }
        }
        ct++;
    }
    file.close();
    index_file.close();
   
}

void BBox2D::clusterJunctionTrajectory(){
	
	if (graph_nodes_!=NULL){
		content_root_->removeChild(graph_nodes_);
	}
    if (cluster_nodes_ != NULL){
        content_root_->removeChild(cluster_nodes_);
    }
    cluster_nodes_ = new osg::Group();

	qDebug() <<"visualize junction trajectories...";
	SubTrajClusters  clusters = junction_network_->getTrajectoryClusters();
	int i=0;
	osg::Vec4 white(0.8f, 0.8f,0.8f,0.2f);
	GPSTrajectories::GpsTrajs gps_trajs = gps_trajectories_->getGpsTrajs();

    for (auto it=clusters.begin(); it!=clusters.end();++it){
		float tint = 0.25;
        int colorIdx = i% Common::palette_size;

        float tr=Common::palette[colorIdx][0]*(1-tint)+tint*255,
            tg = Common::palette[colorIdx][1]*(1-tint)+tint*255,
            tb = Common::palette[colorIdx][2]*(1-tint)+tint*255;
        osg::Vec4 color(tr/255.f, tg/255.f,tb/255.f,0.2f);
		qDebug() << "cluster of  " << it->second.size() << " trajectories. C: "
				 << color[0] <<","<<color[1] <<","<<color[2] <<","<<color[3];

        osg::ref_ptr<osg::Geode> clusterNode = new osg::Geode();
        cluster_nodes_->addChild(clusterNode);

        for (size_t j=0; j<it->second.size();++j){
			subtraj_descriptor subtraj = it->second[j];
			int trajId = std::get<0>(subtraj);
			int pStart = std::get<1>(subtraj);
			int pEnd = std::get<2>(subtraj);			
			//			SubTrajectory subtraj =it->second[j];
			std::shared_ptr<GpsTraj>& traj =gps_trajs[trajId];
            std::vector<cv::Point2f> cvpts,cvptsSub;
			
			for (size_t k= 0;k < traj->point_size();++k){
                cv::Point2f p = toCvPoint2f(&(traj->point(k) ));
				cvpts.push_back(p);
				if ( k >= pStart && k< pEnd){
					cvptsSub.push_back(p);
				}
            }
		   
            //drawPolyLine(clusterNode, cvpts, white,0.15f*resolution_);
			drawPoints(clusterNode, cvptsSub, color, 0.6f*resolution_);
            clusterNode->setNodeMask(0x0);
		}
        i++;
	}

    content_root_->addChild(cluster_nodes_);
   
}
void BBox2D::completeTrajectories(){
	if (junction_network_==nullptr){
		return;
	}
	SubTrajClusters  clusters = junction_network_->getTrajectoryClusters();
	GPSTrajectories::GpsTrajs gps_trajs = gps_trajectories_->getGpsTrajs();
	TrajectoryCompletion trajCompl(&clusters, &gps_trajs,
					 cv::Point2f(min_corner_.x(), min_corner_.y()),
								   getHeight(), resolution_);
    trajCompl.completeAll2("completionResults.csv");
	trajCompl.exportDenseTrajectories("dense_traj_g");
}

void BBox2D::completeSingleTrajectory(int trajId){
    if (junction_network_==nullptr){
        return;
    }
    SubTrajClusters  clusters = junction_network_->getTrajectoryClusters();
    GPSTrajectories::GpsTrajs gps_trajs = gps_trajectories_->getGpsTrajs();
    TrajectoryCompletion trajCompl(&clusters, &gps_trajs,
                     cv::Point2f(min_corner_.x(), min_corner_.y()),
                                  voting_map_.map_image_.rows, resolution_);
    std::vector<cv::Point2f> denseTraj;bool success;
	tie(denseTraj ,success)= trajCompl.completeTrajectory(trajId);
	if (!success)
		qDebug() << "completion failed...";
	if (denseTraj.size()==0)
		qDebug()<<"Error: empty dense trajectory!";
	
    if (dense_traj_nodes_!=NULL){
        content_root_->removeChild(dense_traj_nodes_);
    }

    dense_traj_nodes_ = new osg::Group();
	osg::Vec4 color(1.0f,0.0f,0.0f,0.3f);
    for (size_t i=0; i<denseTraj.size();++i){
        cv::Point2f dir = (i==denseTraj.size()-1)? denseTraj[i]-denseTraj[i-1] :
                denseTraj[i+1]-denseTraj[i];
        cv::Point2f dirN;
        Common::normalize(dir, dirN);
		if (dirN.x==0 && dirN.y==0){
            qDebug()<<"Error: duplicated points at " << i
					<< " th points in dense trajectory " << trajId
					<<" of length " << denseTraj.size();
		}
        osg::Vec3 base = toVec3(denseTraj[i]),
               top = toVec3(denseTraj[i] + 2*dirN);
        osg::ref_ptr<osg::Geode> arrow = OSGUtility::drawCone(base,top,0.6*resolution_,color);
         dense_traj_nodes_->addChild(arrow);
    }

    std::shared_ptr<GpsTraj>& traj =gps_trajs[trajId];
    //std::vector<cv::Point2f> sparseTraj;
	osg::Vec4 color2(1.0f,1.0f,0.0f,0.3f);
    for (size_t k= 0;k < traj->point_size();++k){
        cv::Point2f pt = toCvPoint2f(&(traj->point(k) ));
        //sparseTraj.push_back(p);

		try{
			double angle = (int32_t)( traj->point(k).head());

			osg::Vec3 base = toVec3(pt),
				top = toVec3(pt+2.2*LShapeIterator::UnitDirection(angle));
			osg::ref_ptr<osg::Geode> arrow = OSGUtility::drawCone(base,top,0.8*resolution_,color2);
			dense_traj_nodes_->addChild(arrow);
		}catch(int e){
		}

		
    }
	
	//osg::ref_ptr<osg::Geode> sparseTrajNode = new osg::Geode();
    //drawPoints(sparseTrajNode, sparseTraj, color2, 1*resolution_);
	//dense_traj_nodes_->addChild(sparseTrajNode);
		
    content_root_->addChild(dense_traj_nodes_);
}

void BBox2D::drawSkeletonGraphJunction(){
    //============ DEBUG=====================================
    /*
	CompactSkeletonGraph *graph = process_graph_->getGraph();

	//	graph_nodes_ = new osg::Geode();
	std::vector<ExtendedEdge> edges = process_graph_->getExtendedEdges();
	qDebug() << edges.size() <<" extended edges";
	osg::Vec4 color(158.f/255.f,202.f/255.f, 225.f/255.f,0.3f);
	//	std::vector<cv::Point2f> ext_points;
	std::vector<cv::Point2f> nodes;

	for (int i=0; i< edges.size(); ++i){
		//		if (edges[i].final){
		int head = edges[i].head, tail =edges[i].tail;
		std::vector<int> points = graph->edge_map_[std::make_pair(head,tail)];
		std::vector<cv::Point2f> edge_points;
		if (points.size() >=2){
			int id1 = points[points.size()-1];
			int id2 = points[points.size()-2];
			int id0 = points[0];
			cv::Point2f p1 = graph->points_[id1],
				p2 = graph->points_[id2],p0=graph->points_[id0];
			cv::Point2f pp = (p1+p2)*0.5f;
			qDebug() << i<<": p0="<<p0.x<<"," << p0.y 	<< " p1="<<p1.x<<"," <<p1.y 
					 <<" order=" <<edges[i].nNeighbors;
			nodes.push_back(p1);
			nodes.push_back(p2);
            osg::Vec3 locTXT=toVec3(pp);
            locTXT.z()=3.0f;

			osg::ref_ptr<osgText::Text> text=new osgText::Text;
			text->setPosition(locTXT);
			char buf[33];
			sprintf(buf,"%d",i);
			text->setText(buf);
			text->setColor( osg::Vec4(0.0f,0.0f,0.0f,1.0f));
			text->setCharacterSize(20.f);
			graph_nodes_->addDrawable(text);

			edge_points.push_back(p1);
			edge_points.push_back(p2);
			drawPolyLine(graph_nodes_,edge_points,color,0.25f*resolution_);
			//	}
			//				qDebug()<<"draw node at " << 
		}else{
			qDebug()<<"edges have 0 or 1 point";
		}
	}
	drawPoints(graph_nodes_,nodes,color,0.5f* resolution_);
*/
}

void BBox2D::drawPoints( osg::Geode* nodes, const std::vector<cv::Point2f> &points,
                         osg::ref_ptr<osg::Vec4Array> colors,float r){
    osg::ref_ptr<osg::Vec3Array> bases = new osg::Vec3Array;

    double radius =(r>0)? r: 1.2f * resolution_;
    for (size_t i=0; i< points.size() ; ++i) {
        osg::Vec3 loc1=toVec3(points[i]);
        bases->push_back(loc1);
    }

    if (bases->size() > 0) {
        osg::ref_ptr<osg::Geode> geode =OSGUtility::drawSpheres(bases,radius,colors);
        addDrawablesFromGeode(geode,nodes);
    }
}

void BBox2D::drawPoints(osg::Geode* nodes,const std::vector<cv::Point2f> &points,
						osg::Vec4 color, float r){

    osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
    for (size_t i=0;i< points.size();++i){
        colors->push_back(color);
    }
    drawPoints(nodes,points,colors,r);
}

void BBox2D::drawPolyLine(osg::Geode* nodes,const std::vector<cv::Point2f> &points,	
						  osg::Vec4 color, float r){
	osg::ref_ptr<osg::Vec3Array> bases = new osg::Vec3Array;
	osg::ref_ptr<osg::Vec3Array> tops = new osg::Vec3Array;
	osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
	std::vector<double> radii;
	double corner_radius = (r>0)? r:0.5f * resolution_;

    for (int i=0; i< (int)points.size() - 1; ++i) {
        osg::Vec3 loc1=toVec3(points[i]);
        osg::Vec3 loc2=toVec3(points[i+1]);

		bases->push_back(loc1);
		tops->push_back(loc2);
		radii.push_back(corner_radius);
			
		colors->push_back(color);
	}
	
			
	if (bases->size() >= 1) {//geode->getNumDrawables()>0 ){
		osg::ref_ptr<osg::Geode> geode =  OSGUtility::drawCylinders(bases,
																	tops, radii, 
																	colors); 
		addDrawablesFromGeode(geode,nodes);
	}
}
 osg::Vec3 BBox2D::toVec3(cv::Point2f p1, bool transform){
     float x,y;
     x = (transform)?min_corner_.x()+(p1.x+0.5)*resolution_:p1.x;
     y = (transform)?min_corner_.y()
                     +(voting_map_.map_image_.rows-p1.y+0.5)*resolution_:p1.y;
     osg::Vec3 loc(x,y, 0.0f );
     return loc;
 }

void BBox2D::exportSamples(const char* fname){//std::ofstream &file){
	std::ofstream file(fname, std::ofstream::out);
	if (!file.is_open()){
		qDebug()<<"Error: can not write samples to " << fname << endl;
	}
	qDebug()<<"Writing samples to " << fname <<"..."<< endl;
	file<<"id,type,x,y,x_scaled,y_scaled,sigma" << endl;
	std::vector<SkeletonPoint> S = skeleton_->getSamples();
	cv::Mat sm = skeleton_->getSigma();
	for (size_t i=0; i< S.size(); ++i){
		osg::Vec3 loc = toVec3(S[i].pos);
		
		
		file << S[i].id<<","<< S[i].type <<"," <<S[i].pos.x <<"," << S[i].pos.y<< ","
		  <<  loc[0] <<","<<loc[1]
		  << sm.at<float>(i,0)<< endl;
		
	}
	file.close();
}





 cv::Point2f BBox2D::toMapCoordinates(float x,float y){
	 return cv::Point2f((x-min_corner_[0])/(float)resolution_,//-0.5
						voting_map_.map_image_.rows - (y-min_corner_[1])/(float)resolution_);// + 0.5 );
 }

void BBox2D::addDrawablesFromGeode(osg::Geode* src, osg::Geode* target){
	if (src==NULL || target==NULL){
		qDebug() << "Error: Can not add drawable";
		return;
	}
	int n = src->getNumDrawables();
	for (int i=0;i<n;++i){
		target->addDrawable(src->getDrawable(i));
	}
	target->computeBound();
}
void BBox2D::drawBranches(void){
	//std::vector<std::vector<Candidate> > branches = skeleton_->getBranches();
	if (branch_nodes_ == 0){
		branch_nodes_ = new osg::Geode();
		content_root_->addChild(branch_nodes_);
	}
	int nDrawn = branch_nodes_->getNumDrawables();
	qDebug() << "# of existing drawables: " << nDrawn;
	qDebug() << "# of total branches " << skeleton_->nBranches() ;


	for (int b = 0; b < skeleton_->nBranches();++b){
		std::vector<Candidate> branch = skeleton_->getBranchById(b);
        int colorIdx = b % Common::palette_size;
		osg::Vec4 color( (Common::palette[colorIdx][0])/255, 
						 (Common::palette[colorIdx][1])/255,
						 (Common::palette[colorIdx][2])/255, 1.f);
		std::vector<cv::Point2f> points;
		for (size_t i=0; i<branch.size();++i){
			points.push_back(branch[i].pos);
		}
		drawPolyLine(branch_nodes_, points, color);
		drawPoints(branch_nodes_, points, color);
	}

    drawBridges();
}

void BBox2D::drawBridges(){
	// draw bridges
    std::map<int, std::set<int> > branchGraph = skeleton_->getBranchGraph();
	osg::ref_ptr<osg::Vec3Array> bases = new osg::Vec3Array;
	osg::ref_ptr<osg::Vec3Array> tops = new osg::Vec3Array;
	osg::ref_ptr<osg::Vec4Array> colors = new osg::Vec4Array;
	std::vector<double> radii;

    std::map<int, std::set<int> >::iterator bit;
	for (bit = branchGraph.begin(); bit!= branchGraph.end();  ++bit){

        int branchId0 = *(bit->second.begin()), branchId1;
		std::vector<Candidate> branch0 = skeleton_->getBranchById(branchId0);
		cv::Point2f e1 = branch0[0].pos, e2 = branch0[branch0.size()-1].pos; 
		osg::Vec4 color(0.f ,0.3,0.f,1.f);
        for (auto vit = bit->second.begin();vit !=  bit->second.end(); ++vit){
            // draw line from one end-point of b0 to the point in bk
            // that contains bit->first
            branchId1 = *vit; // (bit->second)[k];
			Candidate cand = skeleton_->getCandidateInBranch(branchId1, bit->first);

			if (cand.id==-1){//x ==-1 && bridgeEnd.y ==-1){
				continue;
			}
			std::ostringstream oss;
			std::vector<Candidate> bpts = skeleton_->getBranchById(branchId1);

			for (int i=0;i<bpts.size();++i){
				oss << bpts[i].id <<" ";
			}
			qDebug() << "Found bridge end " << cand.id << " in branch " << branchId1
					 << (oss.str()).c_str();

			cv::Point2f bridgeEnd =cand.pos;
			cv::Point2f bridgeStart = (cv::norm(e1-bridgeEnd) < cv::norm(e2-bridgeEnd))?
                e1:e2;
            bases->push_back(toVec3(bridgeStart));
            tops->push_back(toVec3(bridgeEnd));
			radii.push_back(0.5f*resolution_);
            colors->push_back(color);
		}
	}

	if (bases->size() > 0){
		osg::ref_ptr<osg::Geode> geode =  OSGUtility::drawCylinders(bases,
																		  tops, 
																		  radii,
																		  colors);
		addDrawablesFromGeode(geode,branch_nodes_);
	} 
	
	return;
	

}

void BBox2D::highlightNeighbors(int skeleton_pid){
    if ( skeleton_==nullptr )
		return;
    // collect all points
    std::vector<SkeletonPoint> K= skeleton_->getSkeleton();	
    osg::ref_ptr<osg::Vec3Array> skeleton_locs = new osg::Vec3Array;
	osg::ref_ptr<osg::Vec4Array> skeleton_colors = new osg::Vec4Array;
	std::vector<double> radii;

	double corner_radius =1.0f * resolution_;

	osg::Vec4 colorNormal =(skeleton_converged_)? 
		osg::Vec4(1.0f, 1.0f, 0.0f, 0.3f): osg::Vec4(1.0f, 0.0f, 1.0f, 0.3f);
    osg::Vec4 colorBridge = osg::Vec4(0.2f, 0.2f, 0.2f, 0.3f);

    for (size_t i=0; i<K.size(); ++i) {
        skeleton_locs->push_back(toVec3(K[i].pos ));
        radii.push_back(corner_radius);
		if (K[i].type == SkeletonPointType::BRIDGE){
            skeleton_colors->push_back(colorBridge);
		}else{
            skeleton_colors->push_back(colorNormal);
		}
	}
    // highlight neighbors
    float r=0.f;
	if (skeleton_pid >= 0){
        osg::Vec4 colorSelected(1.0f,0.0f,0.0f,1.0f),//red
            colorNeighbors(0.0f,1.0f,1.0f,1.0f);//yellow
		skeleton_colors->at(skeleton_pid) = colorSelected;
        r = skeleton_->getAdaptiveH(K[skeleton_pid].pos.x,K[skeleton_pid].pos.y);
        qDebug() << "neighborhood search within range " << r << "..." ;
        L1Skeleton_::SearchResult result = skeleton_->findSkeletonNeighbors(K[skeleton_pid].pos);
        if (result.indices.rows >0 && result.indices.cols >0 ){
            for (int i=0;i<result.indices.cols;++i){
                if (result.indices.at<int>(0,i)==-1)
                    continue;
                int neighborId = result.indices.at<int>(0,i);
               // qDebug() << "found neighbor point " << neighborId <<" dist = "
                //         << cv::norm(K[neighborId].pos - K[skeleton_pid].pos);
                skeleton_colors->at(neighborId) = colorNeighbors;
            }
        }else{
            qDebug() <<"Can not highlight neighbors!";
        }
	}
    // draw all points
	content_root_->removeChild(skeleton_nodes_);
	skeleton_nodes_=  new osg::Geode(*OSGUtility::drawSpheres(skeleton_locs,
															  radii,
															  skeleton_colors));
    // draw disk and arrow
	if (skeleton_pid>=0){

		osg::ref_ptr<osg::Geode> disk = 
			OSGUtility::drawDisk(skeleton_locs->at(skeleton_pid),
								 r*resolution_,
								 osg::Vec3(0.f,0.f,1.f), osg::Vec4(1.f,1.f,0.8f,0.2f));
		addDrawablesFromGeode(disk,skeleton_nodes_);
        cv::Point2f d = skeleton_->getPrincipalDirection(skeleton_pid);
        qDebug() <<"pca direction: "<< d.x<<","<< d.y;
        osg::Vec3 arrowDir(6*resolution_* d.x,
                        -6*resolution_* d.y, 0.f);
        osg::ref_ptr<osg::Geode> arrow = OSGUtility::drawCone(skeleton_locs->at(skeleton_pid),
                                                              skeleton_locs->at(skeleton_pid) + arrowDir,
                                                              resolution_, colorBridge);
        skeleton_nodes_->addDrawable(arrow->getDrawable(0));

        qDebug() <<"Check: Skeleton index= " << skeleton_pid
                 << " id of the point="
                 << K[skeleton_pid].id;
	}
    content_root_->addChild(skeleton_nodes_);

	return;
}



void BBox2D::toggleTrajectoryCluster(int clusterId, bool visible){

	GPSTrajectories::GpsTrajs gps_trajs = gps_trajectories_->getGpsTrajs();

    if (cluster_nodes_!=nullptr && cluster_nodes_->getNumChildren() > clusterId){
        osg::Geode *n = dynamic_cast<osg::Geode*>(cluster_nodes_->getChild(clusterId));
        n->setNodeMask( (visible)? 0xffff :0x0);//  (n->getNodeMask()));
		/*if (visible){
            process_graph_->printClusterInfo(clusterId);
			}*/
    }else{
        if (cluster_nodes_==nullptr){
            qDebug() <<"cluster nodes not exists";
        }else{
            qDebug() <<"cannot get cluster " <<clusterId <<", total # clusters="
                       <<cluster_nodes_->getNumChildren();
        }
    }

}

void BBox2D::displayTrajectoryClusters(bool visible){
	if (cluster_nodes_ == nullptr)
		return;
	for (size_t t = 0; t < cluster_nodes_->getNumChildren();++t){
		osg::Geode *n = dynamic_cast<osg::Geode*>(cluster_nodes_->getChild(t));
        n->setNodeMask( (visible)? 0xffff :0x0);
	}
}

void BBox2D::displayTrajectoryHeading(std::vector<bool> checkedIds, bool visible){
    
	if (heading_nodes_ != NULL){
        content_root_->removeChild(heading_nodes_);
    }
	if (!visible) return;
    qDebug() << "Display trajectory heading! ";
    heading_nodes_ = new osg::Group();
	
    GPSTrajectories::GpsTrajs gps_trajs = gps_trajectories_->getGpsTrajs();
    SubTrajClusters  clusters = junction_network_->getTrajectoryClusters();
   // if (cluster_nodes_->getNumChildren() > clusterId){
		int i=0;
		for (auto it=clusters.begin();it!=clusters.end();++it){
			//        for (size_t i=0;i<checkedIds.size();++i){
            //gps_trajs[i]
			if (checkedIds[i] ){
				//auto c=it->second;
				
				int trajId =std::get<0>( it->second[0]);
                std::shared_ptr<GpsTraj> gps_traj = gps_trajs[trajId];
				// draw arrows
                for (size_t j = 0, j_end =  gps_traj->point_size()-1; j < j_end; ++j){
                    cv::Point2f pt = toCvPoint2f(  & (gps_traj->point(j)));
                    //cv::Point2f pt2 = toCvPoint2f(  & (gps_traj->point(j+1)));
					//toVec
					try{
                        double angle = (int32_t)( gps_traj->point(j).head());
                        osg::Vec4 color(0.2f,0.2f,0.2f,0.3f);
                        osg::Vec3 base = toVec3(pt),
                               top = toVec3(pt+2*LShapeIterator::UnitDirection(angle));
                        osg::ref_ptr<osg::Geode> arrow = OSGUtility::drawCone(base,top,0.6*resolution_,color);
                        heading_nodes_->addChild(arrow);
					}catch(int e){
						qDebug() << "Error: cannot read angle of point " << j <<" in traj" << trajId;
					}
					
				}
			}
			
			i++;
        }
    //}
    content_root_->addChild(heading_nodes_);

}

std::vector<bool> BBox2D::filterClustersBySize(int minSize){
    
    if (junction_network_ == nullptr)
        return std::vector<bool>();

    SubTrajClusters  clusters = junction_network_->getTrajectoryClusters();
      std::vector<bool> mask(clusters.size(),false);
    if (cluster_nodes_ == nullptr)
        return mask;

    int i=0;
   for (auto it=clusters.begin();it!=clusters.end();++it){
       osg::Geode *n = dynamic_cast<osg::Geode*>(cluster_nodes_->getChild(i));
        if (it->second.size()<minSize){
            n->setNodeMask(0x0);
            mask[i] = false;
        }else{
            n->setNodeMask(0xffff);
            mask[i] = true;
        }
        i++;

    }
   return mask;


}
L1Skeleton_*  BBox2D::getSkeleton(void){
      /*if (Parameters::getInstance()->skeleton_use_heading){
           return skeleton_heading_;
      }else{
          return skeleton_;
      }*/
      return skeleton_;
  }
