#pragma once
#ifndef L1_SKELETON_GRAPH_H
#define L1_SKELETON_GRAPH_H

#include <opencv2/core/core.hpp>
#include "l1skeleton.h"
#include <vector>
#include <string>
#include "parameters.h"
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include <boost/graph/breadth_first_search.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/visitors.hpp>
#include <exception>

using namespace boost;
enum class SkeletonGraphVertexType{
    UNKNOWN=0,
        INTERNAL=1,
        TERMINAL=2,
        INTERSECTION=3
        };
enum class SkeletonGraphEdgeType{
    BRANCH=0,
    BRIDGE=1
};



class BfsException: public std::exception
{ 
public:
  BfsException(bool success): m_success(success){ 
  }
  virtual bool success() const throw()
  {
    return m_success;
  }
  bool m_success;
} ;

/**
 * extract skeleton topology from redundant branches
 * as an undirected graph
 */
class L1SkeletonGraph{
 public:
    struct vertex_info{
        int id;
        SkeletonGraphVertexType type;
        SkeletonPoint p;
    };

    struct edge_info{
        int id;
        bool visited;
        float length;
        std::vector<cv::Point2f> internalPts;
    };
    /**
     *filter graph predicate for removing dangling edges
     */
    template <typename Graph>
    struct edge_predicate{
        edge_predicate() { }
         edge_predicate(Graph g) : g(g) { }
         template <typename Edge>
         bool operator()(const Edge& e) const {
           typename Graph::vertex_descriptor s = source(e,g), t = target(e,g);
           return ( g[e].internalPts.size() >3 ||
                    ( g[s].type != SkeletonGraphVertexType::TERMINAL
                   && g[t].type != SkeletonGraphVertexType::TERMINAL) );
         }
         Graph g;
    };

    //template <class Name>

    typedef std::map<int, std::set<int> > BridgeMap;
    typedef std::vector<std::vector<Candidate> > BranchList ;
    typedef adjacency_list < vecS, vecS, undirectedS, vertex_info, edge_info> SGraph;
    typedef graph_traits<SGraph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<SGraph>::edge_descriptor edge_descriptor;
    typedef graph_traits<SGraph>::adjacency_iterator adjacency_iterator;
    typedef graph_traits<SGraph>::out_edge_iterator edge_iterator;
    typedef graph_traits<SGraph>::vertex_iterator vertex_iterator;
    typedef graph_traits<SGraph>::edge_iterator un_edge_iterator;


	/**
	 * visitor class for breadth-first search 
	 * stores distance in m_distance, throws BfsException when
	 * search terminates
	 */
	template <class DistanceMap, class Tag> 
	struct SearchVisitor:public base_visitor<SearchVisitor<DistanceMap,Tag> >{
		typedef Tag event_filter;
	SearchVisitor(DistanceMap pa,int t,int maxDepth) :
		m_distance(pa),m_target(t),m_maxDepth(maxDepth) {}
		template <class Edge, class Graph>
			void operator()(Edge e, const Graph& g) {
			typename graph_traits<Graph>::vertex_descriptor u = source(e, g),
				v = target(e, g);
			if ( get(m_distance,u) +1 > m_maxDepth){
				throw BfsException(false);
			}else{
				put(m_distance, v, get(m_distance, u) + 1);
				if (g[v].id ==m_target){
					throw BfsException(true);
				}
			}
		}


		DistanceMap m_distance;
		int m_target;
		int m_maxDepth ;
	};
	/**
	 * utility function that creates a SearchVisitor object
	 */
	template <class DistanceMap, class Tag>
		inline SearchVisitor<DistanceMap, Tag> searchVisitor(DistanceMap m,
												  int t, int maxDepth,Tag ) {
		return SearchVisitor<DistanceMap, Tag>(m,t,maxDepth);
	}

    class EdgeWriter {
    public:
         EdgeWriter(SGraph _g) : g(_g) {}
         template <class VertexOrEdge>
         void operator()(std::ostream& out, const VertexOrEdge& v) const {
                out << "[label=\"" << g[v].id <<" (" << g[v].length << ") \"]";
         }
    private:
         SGraph g;
    };

    class VertexWriter {
    public:
         VertexWriter(SGraph _g, float _h) : g(_g),height(_h) {}
         template <class VertexOrEdge>
         void operator()(std::ostream& out, const VertexOrEdge& v) const {
             std::string colorStr;
             switch(g[v].type){
             case SkeletonGraphVertexType::TERMINAL:
                 colorStr ="firebrick1";break;
             case SkeletonGraphVertexType::INTERSECTION:
                 colorStr = "gold";break;
             case SkeletonGraphVertexType::INTERNAL:
                 colorStr = "deepskyblue";break;
             default:
                 colorStr = "azure";break;
                }
             float graphScale =  1.f/( 2.0f*Parameters::getInstance()->skeleton_sample_cell_size);

             cv::Point2f pos =  (g[v].p.pos)*graphScale;
                out << "[label=\"" << g[v].id
                    << "\" color=\"" <<  colorStr << "\" pos=\"" <<
                    pos.x <<"," << height *graphScale - pos.y
                    << "!\" ]";
         }
    private:
         SGraph g;
         float height;
    };

   const std::string vertexTypeToString(SkeletonGraphVertexType vtype) {
        switch(vtype){
        case SkeletonGraphVertexType::TERMINAL:
            return  "terminal" ;
        case SkeletonGraphVertexType::INTERNAL:
            return "internal" ;
        case SkeletonGraphVertexType::INTERSECTION:
            return "intersection" ;
        default:
            return "unknown" ;
        }
    }
    float height_;

   /**
     * @brief L1SkeletonGraph
     * @param branches skeleton branch list
     * @param bridges map between bridge points and lists of branches
     * @param bridgePoints positions of bridge points
     */
    L1SkeletonGraph(BranchList branches, BridgeMap bridges,
                    std::vector<SkeletonPoint> bridgePoints,float h  );
    /**
     * @brief getEdgeSegments
     * @return list of polyline coordinates of all edges
     */
    std::vector<std::vector<cv::Point2f> > getEdgeSegments();
    std::vector<cv::Point2f> getVertexPositions();
    SGraph getGraph(){return graph_;}
    void printGraph();
	void writeSkeletonImage(const char* fname,int height,int width);
    void writeDotGraph(const char* fname);
	void updateEdgeIndex();
    void updateEdgeGeometry(int edgeId, std::vector<cv::Point2f> pts);
							//	bool smooth=false);
	void updateEdgeGeometry(edge_descriptor e, std::vector<cv::Point2f> pts);
							//	bool smooth=false);
    std::pair< edge_descriptor,bool> getEdgeById(size_t edgeId);//edgeId: index in edge list
    int getCurrentEdgeId(edge_descriptor e); //get index of a edge_descriptor, -1 if not found

    /**
     * @brief moveVertex move vertex v to new location
     * and update the geometry of its incident edges
     */
    void moveVertex(vertex_descriptor v,  cv::Point2f newLoc);

    /**
     * @brief moveVertex move vertex v to new location
     * do not modify edge

    void updateVertexPos(vertex_descriptor v,  cv::Point2f newLoc);*/
	std::pair<edge_descriptor,bool> addEdge(vertex_descriptor src, vertex_descriptor dst,
                 std::vector<cv::Point2f> pts);
	// void addEdge(int src, int dst, std::vector<cv::Point2f> pts);	
    void removeEdge( vertex_descriptor src, vertex_descriptor dst);
    /**
     * @brief splitEdgeAt split polyline edge e at segPos'th segment by vertex v
     */
    void splitEdgeAt(edge_descriptor e, vertex_descriptor v, int segPos );
    void splitEdgeAtIntersection(edge_descriptor e, vertex_descriptor v, int segPos ,
								 cv::Point2f intersectPt,
								 bool keepFirstHalf, bool keepLastHalf);
    bool isValidEdge(edge_descriptor e);
    bool incidentEdges(edge_descriptor e1, edge_descriptor e2);
    bool incidentEdges(int eid1, int eid2);
    void assignVertexType();
    void removeDanglingEdges( );
    static float cosAngleBetween(const cv::Point2f &v1, const cv::Point2f &v2);
    void updateVertexPos();
    void removeInternalNodes();
	void mergeSimilarBranches(float minLen, float minAngle);
    void computeEdgeWeight();
	void collapseEdgesShorterThan(float minLen);
    float shortestpath(vertex_descriptor src, vertex_descriptor target);
    float shortestpath(vertex_descriptor src, vertex_descriptor target, 	
    				   std::vector<vertex_descriptor> &path);
     void dijkstras(vertex_descriptor vertex_from, std::vector<float>  &distances,
                   std::vector<vertex_descriptor>   &parents );//(vertex_descriptor vertex_from,

               //  std::shared_ptr< std::vector<float> > distances,
                 //  std::shared_ptr<std::vector<vertex_descriptor> > parents);
    //(vertex_descriptor vertex_from, std::vector<float> &distances,
                   //std::vector<vertex_descriptor> &parents);
	 bool bfs(vertex_descriptor src, vertex_descriptor target, int max_depth);
	 bool bfs(int src, int target, int max_depth);
	 
	 std::pair<bool,std::vector<vertex_descriptor> > nearAdjacent(vertex_descriptor v1,
	 															  vertex_descriptor v2);
	//float graphDistance(vertex_descriptor src, vertex_descriptor dst);
protected:
    /**
     * @brief graph_ graph data structure
     */
    SGraph graph_;
    /**
     * @brief vIndex_ maps skeleton_id (key) to index in graph (value)
     */
    std::map<int,int> vIndex_;

    /**
     * @brief addEdgesFromBranches add edges from all branches
     * @param branches
     * @return # of vertices added.
     */
    int addEdgesFromBranches(const BranchList &branches);

    /**
     * @brief addEdgesFromBridges
     * @param branches
     * @param bridges
     * @param bridgePoints
     * @param edgeCount
     * @return
     */
    int addEdgesFromBridges(BranchList branches,
                            BridgeMap bridges,
                            std::vector<SkeletonPoint> bridgePoints,
                            int edgeCount);

    /**
     * @brief addVertex create vertex from skeleton point p and
     * define vertex property.
     * @param i vertex index
     * @param p skeleton point
     * @return vertex descriptor of selected vertex
     */
    vertex_descriptor getVertex(SkeletonPoint p);
    vertex_descriptor getVertex(SkeletonPoint p, SGraph &g,
                                std::map<int,int> &vertexIndex);

    bool vIndexExists(int skeletonId);
    bool vIndexExists(int skeletonId,  std::map<int,int> vertexIndex);

    bool isTurningPoint(SGraph::vertex_descriptor v);

    //void traceSimplePath(vertex_descriptor vertex, edge_descriptor edge,
    //                     std::vector<vertex_descriptor> &edgePoints);
    void traceSimplePath(vertex_descriptor vertex, edge_descriptor edge,
                         std::vector<int> &edgePoints);
    void removePath(std::vector<vertex_descriptor> &edgePoints);
    int nCompactVertices ;
    void initializeEdgeStatus();
	void collapseEdge(edge_descriptor edgeToCollapse);
    //std::vector<int> mergeInternalPts(std::vector<int>  edgePoints);
    std::vector<cv::Point2f> mergeInternalPts(cv::Point2f start, std::vector<int>  edgePoints);
	float computeBranchAngle(vertex_descriptor v, vertex_descriptor u,
							 vertex_descriptor w);

}; 



#endif
