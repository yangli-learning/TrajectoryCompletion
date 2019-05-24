#pragma once
#ifndef TRAJ_COMPLETION_H
#define TRAJ_COMPLETION_H

#include "gps_trajectories.h"
#include "junction_graph.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/index/rtree.hpp>
//class JunctionGraph;
//typedef JunctionGraph::SubTrajClusters SubTrajClusters;
typedef GPSTrajectories::GpsTrajs GpsTrajs;
class TrajPoint;

//--------------------------------------------------------------
// nearest neighbor search typenames and data structure
using namespace boost;
namespace bg = geometry;
namespace bgm = geometry::model;
namespace bgi = geometry::index;
typedef bgm::point<double, 2, bg::cs::cartesian> point;
typedef std::pair<point,size_t> pvalue; //point rtree type {point, pid}
typedef bgi::rtree<point,bgi::rstar<16> > Rtree; // point rtree
typedef bgi::rtree<pvalue,bgi::rstar<16> > prtree;
//--------------------------------------------------------------
// cluster graph typenames and data structure
struct ClusterGraphEdge{
	float length;
};
struct ClusterGraphVertex{
	int id;
	int trajId;
	int pId;
};
//---------------------------------------------------------------
// misc structs
struct CompareClusterSize{
CompareClusterSize(SubTrajClusters* c,
				   std::vector<std::pair<int,int> > cindex):
	clusters(c),clusterIndex(cindex){}
	bool operator() (int c1, int c2){
		std::array<int,2> key1={(clusterIndex)[c1].first,
								(clusterIndex)[c1].second},
			key2={(clusterIndex)[c2].first,
				  (clusterIndex)[c2].second};
        return (*clusters)[key1].size()>= (*clusters)[key2].size();
	}
	SubTrajClusters *clusters;
	std::vector<std::pair<int,int> > clusterIndex;
};
//--------------------------------------------------------------
// nearest neighbors trajectory completion
class TrajectoryCompletion{

 public:
	// cluster graph
	typedef adjacency_list<vecS, vecS, directedS, ClusterGraphVertex,
		ClusterGraphEdge> CGraph;
	typedef graph_traits<CGraph>::vertex_descriptor vertex_descriptor;
	typedef graph_traits<CGraph>::adjacency_iterator adjacency_iterator;
	typedef graph_traits<CGraph>::edge_descriptor edge_descriptor;
	typedef graph_traits<CGraph>::vertex_iterator vertex_iterator;

	TrajectoryCompletion(SubTrajClusters *clusters, GpsTrajs * trajs ,
						 cv::Point2f minCoor, float height, float resolution);
	// complete all trajectories
    void completeAll();//{}

	// complete all trajectories with output to file
	void completeAll2(const char* fname);
	
	std::pair<std::vector<cv::Point2f>,bool> completeTrajectory(int trajId);
    void exportDenseTrajectories(const std::string & basename);

 protected:
	float step_;
	float height_;
	cv::Point2f min_corner_;
	float resolution_;
	
	SubTrajClusters *clusters_;
	GpsTrajs *trajs_;
	
	std::map<int,CGraph> cluster_graphs_;
	
	/**
	 * each index I[c][i][j] is the graph id of jth point on ith subtraj of 
	 * cth cluster in cth cluster graph
	 */
	std::map<int,std::vector<std::vector<int> > > cluster_graph_indices_;

	/**
	 * map trajectory id to 
	 */
	std::map<int, std::vector<int> > trajClusterMap_;
	
	std::vector<std::pair<int,int> > jEdges_;
	std::map<int, std::vector<cv::Point2f> > denseTrajs_;
	void debugTrajClusterMap();
	subtraj_descriptor getSubTrajFromCluster(int clusterId,
											 int trajId);

	/**
	 * @brief fitting trajectory to cluster using greedy method
	 **/
	void fitSubTrajToCluster(subtraj_descriptor st, int cid,
							 std::vector<std::vector<cv::Point2f> > &denseTraj);//, std::vector<subtraj_descriptor> cluster);

	/**
	 * @brief fitting trajectory to cluster using shortest path search
	 */
	void fitSubTrajToClusterSP(subtraj_descriptor st,int cid,
							   std::vector<std::vector<cv::Point2f> > &denseTraj,
							   std::vector<bool> &success);
	void makeClusterGraph(int cid);//,CGraph &g,
	//std::vector<std::vector<int> > &vIndex );
	void debugDenseTrajectories();

	/**
	 * @brief find the stable heading vector of every trajectory point in the cluster
	 */
	//std::vector<cv::Point2f>  findStableHeadings(int cid);
	//	void getDijkstrasPath(const CGraph &g);
	cv::Mat headingHistogram(std::vector<cv::Point2f> h);
	cv::Mat headingHistogram(cv::Point2f v);
	
	std::vector<cv::Point2f > unit_directions_;

    // convert utm coordinates to scene/bbox coordinates
	cv::Point2f toCvPoint2f(const TrajPoint *p);

    //inverse transformation of toCvPoint2f
    cv::Point2f fromCvPoint2f(const cv::Point2f p);
    void exportClusterGraph( );
	void exportClusterGraphCSV( );
};



#endif
