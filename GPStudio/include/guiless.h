#pragma once
#ifndef GUILESS_H
#define GUILESS_H

#include <opencv2/core/core.hpp>
#include <vector>
#include "l1skeleton.h"
#include "skeleton_graph.h"
#include "l1skeleton_graph.h"
#include "junction_network.h"
#include "voting_map_renderer.h"
#include "gps_trajectories.h"
#include <string>
#include <memory>
class TrajPoint;

class Guiless{
 public:
	Guiless(const std::string &trace_fname,
			const std::string &bbox_fname);
	~Guiless();
	int exec();
	void reportTime(double dt, const char* label);
	bool importBBox(const std::string& fname);
	void exportClusters(const std::string &fname,
						const std::string &index_fname);
 private:
	osg::observer_ptr<GPSTrajectories>  gps_trajectories_ ;
	//VotingMap voting_map_;
	osg::Vec3 min_corner_;
	osg::Vec3 max_corner_;
	
	//VotingMapRenderer voting_map_;
	
	float resolution_;
	std::shared_ptr<L1Skeleton>  skeleton_;
	bool valid_;

	std::shared_ptr<JunctionNetwork> junction_network_;

	/**
	 * @brief voting image
	 */
	cv::Mat voting_image_;
	
	/**
	 * @brief heading distribution 
	 */
	std::vector<cv::Mat> heading_images_;

	/**
	 * @brief max heading directions at each pixel
	 */
	std::vector<cv::Point2f> directions_;


	void loadParameters();
};
#endif
