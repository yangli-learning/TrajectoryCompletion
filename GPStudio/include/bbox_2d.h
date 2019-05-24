#pragma once
#ifndef BBOX_2D_H
#define BBOX_2D_H

#include <opencv2/core/core.hpp>
#include <fstream>
#include <vector>
#include "renderable.h"
#include "l1skeleton.h" 
#include "skeleton_graph.h"
#include "l1skeleton_graph.h"
#include "junction_network.h"
//#include "process_graph.h"
#include "voting_map_renderer.h"
#include "gps_trajectories.h"

#define EXPORT_SAMPLES true

class TrajPoint;

class BBox2D: public Renderable {
public:
  BBox2D(GPSTrajectories* scene_image);

  ~BBox2D(void);

  //============= bounding box creation =============
  void setBuilt(bool built);

  bool initialized(void) {
    return corner_number_ != 0;
  }
  void addCorner(const osg::Vec3& corner);
  void clearCorners(void);
  bool isDegenerated(void);

  osg::Vec3 getCenter(void) const {
    return (min_corner_ + max_corner_) / 2;
  }
  const osg::Vec3& getMinCorner(void) const {
    return min_corner_;
  }

  const osg::Vec3& getMaxCorner(void) const {
    return max_corner_;
  }

  float getWidth(void) const {
    return max_corner_.x() - min_corner_.x();
  }

  float getHeight(void) const {
    return max_corner_.y() - min_corner_.y();
  }

  void pickEvent(PickMode pick_mode,float x, float y);


  const cv::Mat& getMapImage(void) const {
    return voting_map_.map_image_;
  }
  const cv::Mat& getMapColorImage(void) const {
    return voting_map_.map_color_image_;
  }

  VotingMapRenderer* getVotingMap(){
      return &voting_map_;
  }
  GPSTrajectories::GpsTrajs  getTrajectories(){
      return gps_trajectories_->getGpsTrajs();
  }
  //============= export and import BBox data with YAML file =
  void exportBBox(const std::string& fname);
  bool importBBox(const std::string& fname);

  //============= build L1 medial skeleton ==================
  L1Skeleton_ *getSkeleton(void);
  float getResolution(){return resolution_;}
  void initializeSkeleton(void);
  void extractSkeleton(void);
  void extractSkeletonStep(void);
  void incrementH0(void);

  void initializeSkeletonWithHeading();

  //============= compute junction and traj clusters ========
  void exportClusters(std::vector<bool> selected, const std::string& fname,
                      const std::string &index_fname);
  void computeSkeletonGraph(float minFlow=-1);
  void clusterJunctionTrajectory();
  inline JunctionNetwork* getJunctionNetwork(){
      return junction_network_;
  }
  //============= trajectory completion ====================
  void completeTrajectories();
  
  //============= visualization =============================
  void toggleRenderMapImage(bool show_map_image);
  void toggleTrajectoryCluster(int clusterId,bool visible);
  void displayTrajectoryClusters(bool visible);
  void displayTrajectoryHeading(std::vector<bool> checkedIds, bool visible);
  void completeSingleTrajectory(int trajId);
  std::vector<bool> filterClustersBySize(int minSize);

  //=============== utility functions ===================
  cv::Point toCvPoint(const TrajPoint* traj_point,bool convert=true);
  cv::Point2f toCvPoint2f(const TrajPoint* traj_point,bool convert=true);
  osg::Vec3 toVec3(cv::Point2f p, bool transform=true);
  cv::Point2f toMapCoordinates(float x, float y);
  void exportSamples(const char* fname);
protected:
  void updateImpl(void);
  void renderMapImage(void);

private:
  //===== private member (bounding box) ====================
  osg::Vec3 min_corner_;
  osg::Vec3 max_corner_;
  int corner_number_;
  bool built_;
  bool skeleton_converged_;
  
  osg::observer_ptr<GPSTrajectories> gps_trajectories_;

  VotingMapRenderer voting_map_;
  bool show_map_image_;

  /**
   * @brief resolution of voting image defined in parameters.ini
   *  can be overwritten by resolution in saved bbox file
   */
  float resolution_;

  //===== private member (l1 skeleton) ====================
  L1Skeleton_  *skeleton_;
  //L1SkeletonWithHeading *skeleton_heading_;
  JunctionNetwork *junction_network_;

  //===== private member (graphics objects) ===============
  osg::ref_ptr<osg::Geode> skeleton_nodes_;
  osg::ref_ptr<osg::Geode> branch_nodes_;
  osg::ref_ptr<osg::Geode> graph_nodes_;
  osg::ref_ptr<osg::Group> cluster_nodes_;
  osg::ref_ptr<osg::Group> heading_nodes_;
  osg::ref_ptr<osg::Group> dense_traj_nodes_;

  //===== private member functions ========================
  void drawBranches();
  void drawPolyLine(osg::Geode* nodes, const std::vector<cv::Point2f> &points, 
					osg::Vec4 color,float r=-1);
  void drawPoints(osg::Geode* nodes, const std::vector<cv::Point2f> &points,
				  osg::Vec4 color,float r=-1);

  void drawPoints(osg::Geode* nodes, const std::vector<cv::Point2f> &points,
                  osg::ref_ptr<osg::Vec4Array> colors,float r=-1);
  void drawBridges();
  void drawSkeletonGraph(CompactSkeletonGraph &graph);
  void drawSkeletonGraph(L1SkeletonGraph &graph);

  void drawSkeletonGraphJunction();
  void traceBranchFromPoint(int skeleton_id);  
  int getPickedSkeletonPoint(float x, float y);
  void addDrawablesFromGeode(osg::Geode* src, osg::Geode* target);  
  void highlightNeighbors(int skeleton_pid);


};

#endif // !BBOX_2D_H
