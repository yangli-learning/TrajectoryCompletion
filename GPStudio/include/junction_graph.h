#pragma once
#ifndef JUNCTION_GRAPH_H
#define JUNCTION_GRAPH_H

#include <QObject>

#include <vector>
#include <opencv2/core/core.hpp>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/filtered_graph.hpp>
#include <boost/graph/graph_utility.hpp>
#include "osg_utility.h"
#include "l1skeleton_graph.h"

using namespace boost;

typedef std::tuple<int,int,int> subtraj_descriptor;
typedef std::vector<subtraj_descriptor> SubTrajCluster;
typedef std::map<std::array<int,2>, SubTrajCluster> SubTrajClusters;

/**
 * obtain the junction graph by computing the pseudo dual graph
 */
class JunctionGraph{
public:

    /**
     * @brief The vertex_info struct
     * Represent a skeleton branch segment in a given direction. Graphically,
     * it contains two endpoints head, tail. Head connects to outgoing traffic
     * Tail connects to incoming traffic. It is linked with another
     * vertex pairId, which share the same branch in opposite direction
     */
    struct vertex_info{
        int id; // consistent with index in graph
        int src;
        int target;
        int pairId;
        std::vector<cv::Point2f> shape;
    };

    /**
     * @brief The edge_info struct
     *  Represents traffic between two road segments, in_seg
     */
    struct edge_info{
        //FlowDirection in_dir;
        //FlowDirection out_dir;
        int flow;
        std::vector<subtraj_descriptor> sub_trajs;
    };

    typedef adjacency_list < vecS, vecS, directedS, vertex_info, edge_info> JGraph;
    typedef graph_traits<JGraph>::vertex_descriptor vertex_descriptor;
    typedef graph_traits<JGraph>::vertex_iterator vertex_iterator;
    typedef graph_traits<JGraph>::edge_descriptor edge_descriptor;
    typedef graph_traits<JGraph>::out_edge_iterator o_edge_iterator;
    typedef graph_traits<JGraph>::in_edge_iterator i_edge_iterator;
	typedef graph_traits<JGraph>::edge_iterator edge_iterator;

	/**
	 * @brief The edge_predicate struct
	 * Whether edge weight (flow) is above a given threshold
	 */
	template <typename Graph>
	struct filter_edge_predicate{
		filter_edge_predicate(){}	
	filter_edge_predicate(Graph _g, float _thresh):g(_g),thresh(_thresh){}
		template <typename Edge>
		bool operator()(const Edge& e) const{
			return g[e].flow > thresh;
		}
		Graph g;
		float thresh;
	};
	

    class EdgeWriter {
    public:
         EdgeWriter(JGraph _g) : g(_g) {}
         template <class VertexOrEdge>
			 void operator()(std::ostream& out, const VertexOrEdge& e) const;
    private:
         JGraph g;
    };

    class VertexWriter {
    public:
         VertexWriter(JGraph _g, float _h) : g(_g),height(_h) {}
         template <class VertexOrEdge>
         void operator()(std::ostream& out, const VertexOrEdge& v) const;
    private:
         JGraph g;
         float height;
    };
    class GraphWriter{
    public:
        GraphWriter( JGraph _g ):g(_g){}

        void operator()(std::ostream& out) const ;
    private:
         JGraph g;

    };

    JunctionGraph();
    JunctionGraph(L1SkeletonGraph* g);
    void init(L1SkeletonGraph::SGraph g);

    bool pairVertices(vertex_descriptor v1, vertex_descriptor v2);
    //vertex_descriptor getVertexById(int i);
    bool vertexExists(int srcId, int dstId);
    void writeDotGraph(const char* fname, float height);
	bool voteEdge(int s1, int t1, int s2, int t2, bool nearIncident,
				  subtraj_descriptor subtraj);
	bool addEdge(int s1, int t1, int s2, int t2, subtraj_descriptor subtraj);
    //getGraph
	//	SubTrajClusters getSubTrajClusters();
	
	void getSubTrajClusters(SubTrajClusters &clusters);
	void filterEdgeByVotes();
	float getMinFlowSize(){return minFlowSize_;}
	void setMinFlowSize(float f){minFlowSize_=f;}
protected:
    JGraph graph_;
    // maps (src, dst) to  vId
    std::map<std::pair<int,int> , int > vIndex;

	/**
	 * use Otsu's algorithm to compute optimal value
	 * for minimum cluster size
	 */
	float getOptFlowThresh();
    // maps (src, intersection, dst) to edgeId
    //std::unordered_map<std::tuple<int,int,int> > eIndex;
	float minFlowSize_;
	
	/**
	 * return true if v1 and v2 are nearly incident
	 * i.e. skeleton branch represented by v1 and v2 are nearly 
	 * incident
	 
	bool nearlyAdjacent(vertex_descriptor v1, vertex_descriptor v2);*/
	
	
};

#endif
