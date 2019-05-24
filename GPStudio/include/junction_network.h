#pragma once
#ifndef JUNCTION_NETWORK_H
#define JUNCTION_NETWORK_H

#include <opencv2/core/core.hpp>

#include <memory>
#include <algorithm>
#include "l1skeleton_graph.h"
#include "junction_graph.h"
#include <vector>
#include <QDebug>
#include <unordered_set>
#include <QStringList>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point.hpp>
#include <boost/geometry/geometries/segment.hpp>
#include <boost/geometry/index/rtree.hpp>
#include "osg_utility.h"



class L1Skeleton_;
class GpsTraj;
class TrajPoint;
class GPSTrajectories;
//class JunctionGraph;
//-----------auxilliary data structure -----------------

enum EndpointType{
    HEAD,
    TAIL,
    UNKNOWN
};
/**
 * @brief the ExtendedEdge class
 * meta information to keep track of direction of edge extension
 */
class ExtendedEdge{
public:
    int id;
    cv::Point2f headDir;
    cv::Point2f tailDir;
    bool headFinal;
    bool tailFinal;
    bool invalid;
    std::vector<cv::Point2f> pts;
    int nNeighbors;

    ExtendedEdge();
    ExtendedEdge(int _id, cv::Point2f _headDir, cv::Point2f _tailDir,
                 bool _headFinal, bool _tailFinal, std::vector<cv::Point2f> _pts);

    ~ExtendedEdge();

    /** copy constructor **/
    ExtendedEdge(const ExtendedEdge &ee);
    /** assignment operator **/
    ExtendedEdge & operator=(const ExtendedEdge &ee);
};

/**
 * @brief descrition of a subtrajectory of trajId,
 * with point range pStart .. pEnd;
 */
struct SubTrajectory{
    int trajId;
    int pStart;
    int pEnd;
};


typedef L1SkeletonGraph::edge_descriptor edge_descriptor;
typedef L1SkeletonGraph::vertex_descriptor vertex_descriptor;

struct TrajProjection{
	int pid; //index of point on trajectory
	edge_descriptor edge; // projected edge on skeleton graph
	cv::Point2f pos; //projection point
	bool isForward; //projection points from tail to head of the edge 
	
};



/**
 * trajectory clusters map type definition
 * key: array<int,2> = {srcBranchId, dstBranchId} cluster identifier
 * value: vector<SubTrajectory> = trajectories belong to this cluster
 
typedef std::vector<SubTrajectory> TrajectoryCluster;
typedef std::map<std::array<int,2>, TrajectoryCluster > TrajectoryClusters;
**/
class TrajPoint;
class GPSTrajectories;

//----------- Boost::Geometry namespaces and typedefs ------
namespace bg = boost::geometry;
namespace bgm = boost::geometry::model;
namespace bgi = boost::geometry::index;
typedef bgm::point<double, 2, bg::cs::cartesian> point;
typedef bgm::segment<point> segment;
typedef std::tuple<segment, size_t, size_t> value; //segment rtree type {lineSegment,edgeId,locOnEdge}
typedef bgi::rtree<value, bgi::rstar<16> > rtree; //segment rtree
typedef std::pair<point,size_t> pvalue; //point rtree type {point, id}
typedef bgi::rtree<pvalue,bgi::rstar<16> > prtree; // point rtree

/**
 * @brief The EdgeIntersection struct
 */
struct EdgeIntersection{
    value qSeg, nSeg;
    EndpointType qType;
    cv::Point2f loc;
    // defer resolution to edge
    //edge_descriptor qEdge; // query edge
    //edge_descriptor nEdge; // neighbor edge
    //int nPos; // neighbor edge position (intersection is between
              // nEdge.pts[nPos] and nEdge.pts[nPos+1]
    //EndpointType qEndpointType; //HEAD or TAIL

    EdgeIntersection(value qSeg, value nSeg,
                     EndpointType type, cv::Point2f loc):
        qSeg(qSeg),nSeg(nSeg),qType(type), loc(loc){

    }
    EdgeIntersection():
         qType(UNKNOWN),loc(cv::Point2f(0.f,0.f)){

    }
};
//------------ satic utility functions ----------------------
static bool pointEqual(point p1, point p2){
	return (p1.get<0>() == p2.get<0>())&& (p1.get<1>()==p2.get<1>());
}
static cv::Point2f boostToCVPoint(point p){
	return cv::Point2f(p.get<0>(), p.get<1>());
}

/**
 * @brief The JunctionNetwork class
 */
class JunctionNetwork{
 public:

	/**
	 *  constructor
	 */
    JunctionNetwork(L1SkeletonGraph *g, L1Skeleton_ *sk, float height, float width);

    /**
     * get all extended edges
     */
    inline std::vector<shared_ptr<ExtendedEdge> > getExtendedEdges(){
        return extended_edges_;
    }
    /**
     * get skeleton graph
     */
    inline	L1SkeletonGraph * getGraph(){return graph_;}

	/**
	 * process skeleton graph and find endpoint clusters
	 * iteratively extend all edges and check new intersections
	 */
	void process();

	/**
	 * initialize edges in the junction graph
     */
    void initExtendedEdges();

    /**
     * initialize junction graph from skeleton graph
     */
    void initJunctionGraph();

	/**
	 * extend the second endpoint of each edge by len
	 */
    bool growExtendedEdges(float len);

	/**
	 * check intersection of all 
	 */
	void checkIntersections();

	/**
	 * cluster endpoints of branches
	 */
	void clusterEndpoints(float r);

	/**
	 * project trajectory to skeleton edge and validate junction
	 * topology
	 */
	void validateTopologyByGPS(osg::observer_ptr<GPSTrajectories> trajectories,
                               osg::Vec3  min_corner,  float resolution,
							   float minFlow);

	/**
	 * get number of subtrajectory clusters
	 */
    int getNumClusters(){return clusters_.size();}

	/**
	 * get max size of subtrajectory cluster
	 */
    int getMaxClusterSize();

	/**
	 * get constant reference to subtrajectory clusters
	 */
     const SubTrajClusters  & getTrajectoryClusters(){return clusters_;}

	/**
	 * convert cluster info to QStringList (to be displayed in 
	 * the GUI selection list)
	 */
    QStringList clusterToQStringList();

	/**
	 * get position of current skeletongraph intersections
	 */
    std::vector<cv::Point2f> getIntersectionPoints();
	static cv::Point2f toCvPoint2f(const TrajPoint* tp, cv::Point2f min_corner,
								   float height, float resolution);
	/**
	 * convert trajectory point to Boost::geometry point format
	 *  tp: pointer to TrajPoint pbject
	 *  min_corner: lower left corner of the bounding box
	 *  resolution: resolution of the voting image
	 */
	point toBoostPoint(const TrajPoint* tp,
					   osg::Vec3 min_corner,float resolution);
	
    cv::Point2f toCvPoint2f(const TrajPoint* tp, osg::Vec3 min_corner, float resolution);
	
 protected:
	//--- protected fields ------------------------

	/** junction extention step **/
	float step_;   
    float min_vertex_density_;
	float min_edge_length_;
	float max_missing_edge_length_;
	float min_edge_angle_;
    /** skeleton graph **/
    L1SkeletonGraph *graph_;

    /** junction graph -- directed dual of skeleton graph**/
    JunctionGraph jgraph_;

	/** bounding box width and height **/
	float width_, height_;

	/** list of extended edges **/
    std::vector<shared_ptr<ExtendedEdge> > extended_edges_; //indexed by edge Id in sgraph
    std::list<shared_ptr<EdgeIntersection> > intersections_; //list of edge intersections

    /** pointer to l1skeleton **/
    L1Skeleton_ *skeleton_;

	/** subtrajectory cluster **/
	SubTrajClusters clusters_;
	
	/** debugging variables **/
	int n_direct_votes_;
	int n_indirect_votes_;
	
	/** endpoint clusters 
	 * clusters_[i][j]: id of jth extended edge at junction i
	 */
    //std::vector<std::vector<size_t> > clusters_;

	/** adjacency matrix at each junction **/
    //std::vector<cv::Mat> cluster_adj_;


	//--- projected methods -----------------------
    bool extendEdge(shared_ptr<ExtendedEdge> ee, EndpointType type, float len);

    cv::Point2f updateEndpointDirection(cv::Point2f p0, cv::Point2f oldDir) ;

    bool isValidExtension(cv::Point2f newHead);
	/**
	 * compute projection of a point on a line segment
	 *  in: point to be projected
	 *  p1: first endpoint of segment
	 *  p2: 2nd endpoint of segment
	 *  out: projected point (output) 
	 */
    bool projectPointOnSegment(cv::Point2f in, cv::Point2f p1,
                           cv::Point2f p2, cv::Point2f &out);

    bool projectPointOnLine(cv::Point2f in, cv::Point2f p1,
                           cv::Point2f p2, cv::Point2f &out);
	/**
	 * make r tree of all extended edges for knn search
	 */
    rtree makeRTreeFromEdges();

	/**
	 * vote edge from edge i to j
	 */
    //bool vote(int i, int j){ return false;}

	//    std::vector<std::pair<int,edge_descriptor> >
	std::vector<TrajProjection> projectGPSToGraph(
				   std::shared_ptr<GpsTraj> traj,
                   const cv::Mat &candidateSegs,  const cv::Mat &candidatePos,
                   const cv::Mat &distances, const osg::Vec3 &min_corner);
	
	void  voteJunctionGraphNaive(int trajId,
				   std::shared_ptr<GpsTraj> traj,
                   const cv::Mat &candidateSegs,  const cv::Mat &candidatePos,
								 const cv::Mat &distances, const osg::Vec3 &min_corner,
								 int maxR);													   
	void vote(int trajId,std::vector<TrajProjection> projections,
			  std::vector<int> breakpts,int nPoints);									   
	void voteAlt(int trajId,std::vector<TrajProjection> projections,
			  std::vector<int> breakpts,int nPoints);

	/**
	 * postprocess trajectory clusters to eliminate incorrect
	 * mappings
	 */
    void processTrajClusters(){}
	/**
	 * comparator for segment searchz: return true if the id of 
	 * segment v is not i 
	 */
    inline bool different_id(value const& v, value const & u){
        return std::get<1>(v)!= std::get<1>(u);
    }
    void initializeSegmentRTree(std::vector<value> &queries, rtree &rt);
    const std::string valueToStr(const value &v) &;

    //inline bool different_id(value const& v, value const &u){
    //	return v.second != i;
    //}
    /**
     * @brief processIntersections
     */
    void processIntersections();

    /**
     * @brief get the vertex descriptor of the endpoint in the query segment
     * in an intersection
     */
    vertex_descriptor qEndpoint(const EdgeIntersection &intersect);
    bool inBBox(cv::Point2f p);

    /*=============== handle trajectory projection ==================*/
    float computeEmissionCost(float dist, float mean, float var);
    float computeTransitionCost(float pathlen, float d2, float lambda=100.f);

    /**
     * @brief backtrace get shortest path from v_from to v_to
     * @param parents
     * @param from
     * @param to
     * @return
     */
    std::list<edge_descriptor> backtrace(std::vector<vertex_descriptor> parents,
                               vertex_descriptor v_from, vertex_descriptor v_to);

    std::vector<std::vector<cv::Point2f> > computeCandidateProjections();
    float computePathLen(const cv::Point2f  &proj_j, const cv::Point2f &proj_k,
                         const std::vector<float> &distances1,
                         const std::vector<float> &distances2,
                         const std::vector<vertex_descriptor> &parents1,
                         const std::vector<vertex_descriptor> &parents2,
                         int edgeId, int edgeId_next, int segId, int segId_next,
                          std::list<edge_descriptor> & shortestpath );

	void processVotes(int trajId,std::shared_ptr<GpsTraj> traj,
					  const cv::Mat &candidateSegs,const cv::Mat &candidatePos,
					  const std::vector<TrajProjection> &edges);

	/**
	 * @brief construct TrajProjection object
	 * @param pid trajectory point id
	 * @param edge descriptor of projected edge on skeletong graph
	 * @param pos projection coordinates
	 * @param isForward whether projected path is in the same direction as the 
	 *        edge polyline
	 */
	TrajProjection makeTrajProjection(int pid,	edge_descriptor edge,
									  cv::Point2f pos, bool isForward);

	void computeSubTrajClusters();

	//*============== debugging functions =====================*/
	void debugDijkstras(const std::vector<float> &distances,
				   const std::vector<vertex_descriptor> &parents);
	void debugEdgeVotes(int id, std::vector<std::vector<int> > votes);	
};

#endif
