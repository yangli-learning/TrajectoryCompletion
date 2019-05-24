#pragma once
#ifndef SKELETON_GRAPH_H
#define SKELETON_GRAPH_H

#include <opencv2/core/core.hpp>

#include "l1skeleton.h"
#include "dgraph.h"
#include <vector>


enum class SkeletonVertexType{
	UNKNOWN=0,
		INTERNAL=1,
		TERMINAL=2,
		INTERSECTION=3
		};

typedef std::pair<int,int> EdgeMapKey ;
typedef std::map<int, std::set<int> > BridgeMap;
typedef std::vector<std::vector<Candidate> > BranchList ;

/**
 * class SkeletonGraph
 * stores the topology and geometry of map skeleton 
 */
class CompactSkeletonGraph: public DGraph{


 public:
	/**
	 * public memebers
	 */
	std::vector<int> vertex_map_;  //vertex_map_[i] is the "unique" id of ith node
	std::map<EdgeMapKey, std::vector<int> > edge_map_; //map edge to polyline
	std::vector<cv::Point2f> points_; //coordinates of all points (vertices and edges)

	/**
	 * default constructor
	 */
 CompactSkeletonGraph(int n):DGraph(n){}

	/**
	 * destructor
	 */
	//	~CompactSkeletonGraph();

	/**
	 *  copy constructor
	 */
	//	CompactSkeletonGraph(const CompactSkeletonGraph &graph);

	/**
	 * assignment operator
	 */ 
	//	CompactSkeletonGraph & operator = (CompactSkeletonGraph rhs);

	/** 
	 * utility functions
	 */
	int	getVertexId(int pointId);
	inline std::vector<int> getEdgePointIds(size_t head, size_t tail){
		return edge_map_[std::make_pair(head,tail)];
	}
	void printVertexMap();
	void printEdgeMap();
	void printVertex(int i);
	// protected:
	// edgeCurves: polyline representation of each curve;
};


/**
 * extract skeleton topology from redundant branches
 */
class SkeletonGraph: public DGraph{
 public:
	/**
	 * constructor with input collection of branches
	 */
	SkeletonGraph(BranchList branches,				  
                  BridgeMap bridges,
				  std::vector<SkeletonPoint> bridgePoints,
				  int n);

	void initialize( BranchList branches,
					 BridgeMap bridges, 
					 std::vector<SkeletonPoint> bridgePoints);

	void buildCompactGraph(CompactSkeletonGraph *compact);
	//	void removeRedundantEdges(CompactSkeletonGraph *compact);
	void printVertexMap();

	/**
	 * trace along internal edges 
	 * parameters
	 *   - edge: starting edge to trace
	 *   - eddgePoints: list of vertex indices on the traced path
	 */
	int getVertexId(int pointId);

	/**
	 * trace simple path
	 */
	void traceSimplePath(DGraphEdge *edge, std::vector<int> &edgePoints);
	int nCompactVertices;

 protected:
	
	//	vector<vector<Candidate> > branches_;
	std::unordered_map<int,int> reverseVertexIndex_;
	std::vector<SkeletonPoint> vertices_;
	std::vector<SkeletonVertexType> vertex_types_;
	void findDanglingEdges(CompactSkeletonGraph* compact, 
							 std::vector<std::pair<int,int> > &edgeToRemove);

	bool isTurningPoint(int i);

	inline int vIndex(int id){
		if (reverseVertexIndex_.find(id)==reverseVertexIndex_.end())
			return -1;// reverseVertexIndex_;
		return reverseVertexIndex_[id];
	}
	inline int vIndex(SkeletonPoint p){
		return vIndex(p.id);
	}
	/**
	 * return incremented id
	 */
    int addVertex(int i,SkeletonPoint p);

	/**
	 * return the total number of actual vertices after processing branches
	 * actual vertices are added using addVertex(int,SkeletonPoint)
	 */
	int addEdgesFromBranches(const BranchList &branches);

	/**
	 * parameters:
	 *   ct: number of actual vertices before 
	 * return  number of actual vertices after processing bridges
	 */
	int addEdgesFromBridges(BranchList branches, 
							 BridgeMap bridges,
							 std::vector<SkeletonPoint> bridgePoints,
							 int ct);

	void assignVertexType();
	//----debug functions ---
	bool containsDuplicates(std::vector<int> endPoints);
	//	CompactSkeletonGraph *compact_graph_;

}; 

/**
 * 
 */


#endif
