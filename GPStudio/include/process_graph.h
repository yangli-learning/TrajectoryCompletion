#pragma once
#ifndef PROCESS_GRAPH_H
#define PROCESS_GRAPH_H

#include <opencv2/core/core.hpp>
#include "skeleton_graph.h"
#include <vector>
#include <unordered_set>
#include <QStringList>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>

#include "osg_utility.h"
//-----------auxilliary data structure -----------------
/**
 * extended version of edge (i,j)->{v0,...,vk} 
 * contains mapping to its twin (j,i)
 */
struct ExtendedEdge{
	int head, tail; // vertec id of edge head and tail 
	//  ExtendedEdge *twin;
	cv::Point2f dir; // direction of extension
	bool final;
	int nNeighbors;
};

/**
 * descrition of a subtrajectory of trajId,
 * with point range pStart .. pEnd;
 */
struct SubTrajectory{
    int trajId;
    int pStart;
    int pEnd;
	//	int keyHead,keyTail;
};
/**
 * trajectory clusters map type definition
 * key: array<int,2> = {srcBranchId, dstBranchId} cluster identifier
 * value: vector<SubTrajectory> = trajectories belong to this cluster
 */
typedef std::vector<SubTrajectory> TrajectoryCluster;
typedef std::map<std::array<int,2>, TrajectoryCluster > TrajectoryClusters;

class TrajPoint;
class GPSTrajectories;

//----------- Boost::Geometry namespaces and typedefs ------
namespace bg = boost::geometry;
namespace bgm = boost::geometry::model;
namespace bgi = boost::geometry::index;
typedef bgm::point<double, 2, bg::cs::cartesian> point;
typedef bgm::segment<point> segment;
typedef std::pair<segment, size_t> value;
typedef bgi::rtree<value, bgi::rstar<16> > rtree;
typedef std::pair<point,size_t> pvalue;
typedef bgi::rtree<pvalue,bgi::rstar<16> > prtree;

//------------ satic utility functions ----------------------
static bool pointEqual(point p1, point p2){
	return (p1.get<0>() == p2.get<0>())&& (p1.get<1>()==p2.get<1>());
}
static cv::Point2f boostToCVPoint(point p){
	return cv::Point2f(p.get<0>(), p.get<1>());
}

/**
 * class: ProcessGraph
 * construct trajectory junction graph. cluster subtrajectories
 */
class ProcessGraph{

 public:


	/**
	 *  constructor
	 */
	ProcessGraph(CompactSkeletonGraph *g, float height, float width);
    
	//~ProcessGraph();

	/**
	 * process skeleton graph and find endpoint clusters
	 * iteratively extend all edges and check new intersections
	 */
	void process();

	/**
	 * initialize edges in the junction graph
	 */
	void initJunctionEdges();

	/**
	 * extend the second endpoint of each edge by len
	 */
	bool extendJunctionEdges(float len);

	/**
	 * check intersection of all 
	 */
	void checkIntersections();

	/**
	 * cluster endpoints of branches
	 */
	void clusterEndpoints();

	/**
	 * project trajectory to skeleton edge and validate junction
	 * topology
	 */
	void validateTopologyByGPS(osg::observer_ptr<GPSTrajectories> trajectories,
							   osg::Vec3  min_corner,  float resolution);

	/**
	 * get number of subtrajectory clusters
	 */
    int getNumClusters(){return trajectory_clusters_.size();}
    int getMaxClusterSize();
	/**
	 * get all extended edges
	 */
	inline std::vector<ExtendedEdge> getExtendedEdges(){
		return extended_edges_;
	}
	/**
	 * get skeleton graph
	 */
	inline	CompactSkeletonGraph * getGraph(){return graph_;}

	/**
	 * get trajectory clusters
	 */
	TrajectoryClusters  getTrajectoryClusters(){return trajectory_clusters_; }

	/**
	 * convert junction graph to dot string (graphviz) 
	 */
	std::string toDotString();

	/**
	 * export junction graph to dot file
	 */
	void exportGraph(const char* fname);

	/**
	 * print information of a cluster
	 */
	void printClusterInfo(int clusterId);


	/**
	 * get id of trajectories in a given cluster
	 */
	std::vector<int> getClusterTrajectoryIds(int clusterId);

	/**
	 * convert cluster info to QStringList (to be displayed in 
	 * the GUI selection list)
	 */
    QStringList clusterToQStringList();
	//clusterS
    /**
     * @brief getClusterKeys
     *
     */
	std::vector<std::array<int,2> >  getSortedClusterKeys();
	
 protected:
	//--- protected fields ------------------------

	/** junction extention step **/
	float step_;   

	/** compact skeleton graph **/
	CompactSkeletonGraph *graph_;

	/** bounding box width and height **/
	float width_, height_;

	/** list of extended edges **/
	std::vector<ExtendedEdge> extended_edges_;

	/** endpoint clusters 
	 * clusters_[i][j]: id of jth extended edge at junction i
	 */
	std::vector<std::vector<size_t> > clusters_;

	/** adjacency matrix at each junction **/
	std::vector<cv::Mat> cluster_adj_;

	/** trajectory that belong to each edge in junction graph
	 * (b1,b2),{[t,pstart,pend]}
	 */
    TrajectoryClusters trajectory_clusters_;

	//--- projected methods -----------------------
	/**
	 * compute projection of a point on a line segment
	 *  in: point to be projected
	 *  p1: first endpoint of segment
	 *  p2: 2nd endpoint of segment
	 *  out: projected point (output) 
	 */
	bool projectPointOnSegment(cv::Point2f in, cv::Point2f p1,
							   cv::Point2f p2, cv::Point2f &out);
	
	/**
	 * make r tree of all extended edges for knn search
	 */
	rtree makeRTreeFromEdges();

	//	void splitExtenededEdge(int v1, int v2, int v3, cv::Point2f insertPt);
	//	void addJunctionEdges();

	/**
	 * convert trajectory point to Boost::geometry point format
	 *  tp: pointer to TrajPoint pbject
	 *  min_corner: lower left corner of the bounding box
	 *  resolution: resolution of the voting image
	 */
	point toBoostPoint(const TrajPoint* tp,
					   osg::Vec3 min_corner,float resolution);

	/**
	 * vote edge from edge i to j
	 */
	bool vote(int i, int j);

	/**
	 * get the opposite edge (twin) of edge i
	 */
	std::vector<ExtendedEdge>::iterator getTwinExtendedEdge(int i);

	/**
	 * compute distance ratio between the nearest few search results
	 * used for choosing projected edge
	 */
	float computeDistRatio(point p, std::vector<value> result_n);

	/**
	 * cluster subtrajectories
	 */
	void clusterSubTrajectories(std::vector<std::array<int,3> > breakpt,
								std::vector<int> project_idx, int i);
	/**
	 * postprocess trajectory clusters to eliminate incorrect
	 * mappings
	 */
	void processTrajClusters();
	/**
	 * comparator for segment searchz: return true if the id of 
	 * segment v is not i 
	 */
	inline bool different_id(value const& v, size_t i){
		return v.second != i;
	}

	/**
	 * test if edge i and j are twin edges
	 */
	inline bool twinEdge(int i, int j){
		int twinId = (getTwinExtendedEdge(i)-extended_edges_.begin());
		return twinId == j;
	}
	
};

#endif
