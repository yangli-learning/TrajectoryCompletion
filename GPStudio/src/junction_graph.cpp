#include "junction_graph.h"
#include <fstream>
#include <iostream>
#include <QDebug>
#include <boost/graph/graphviz.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <utility>

#include <boost/graph/copy.hpp>

template <class VertexOrEdge>
 void JunctionGraph::VertexWriter::operator()(std::ostream& out, const VertexOrEdge& v) const

{
    float graphScale =  1.f/(1.9f*Parameters::getInstance()->skeleton_sample_cell_size);
    int curveLen = g[v].shape.size();
    int mid = floor(curveLen/2);
    cv::Point2f pos = g[v].shape[mid];

    if ( g[v].id  > g[v].pairId){
        pos.y += 0.75*Parameters::getInstance()->skeleton_sample_cell_size;
    }else{

       pos.y -= 0.75*Parameters::getInstance()->skeleton_sample_cell_size;
    }
    pos *= graphScale;
    out << "[label=\"" << g[v].id //<< "\" ]";
         <<"\" pos=\" " << pos.x <<"," << height *graphScale - pos.y
          << "!\" ]";

}

template <class VertexOrEdge>
void JunctionGraph::EdgeWriter::operator()(std::ostream& out, const VertexOrEdge& e) const {
	int minClusterSize = Parameters::getInstance()->skeleton_min_cluster_size;
	std::string colorStr = (g[e].flow > minClusterSize)? "\"red\"":"\"gray\"" ;
	out << "[label= \" (" << g[e].flow << ") \" color="
		<< colorStr
		<< " ]";
}

 void JunctionGraph::GraphWriter::operator()(std::ostream& out ) const
{
     //out <<
     vertex_iterator vit,vit_end;
     int clusterId=0;

     //for (auto it=vIndex.begin();it!=vIndex.end();++it){
     for (tie(vit,vit_end)=vertices(g);vit!=vit_end;++vit){
         int id = g[*vit].id; int pairId = g[*vit].pairId;
         //int src= it->first.first;
         //int target = it->first.second;
         //int vid = it->second;
         if (id > pairId) continue;
         out <<"subgraph cluster_" << clusterId<<"{ style=filled;color=lightgrey;"
         //vertex_descriptor v = vertex(vid,graph_);
           << id <<"->" <<  pairId <<" [dir=both]}\n";
     }

 }



JunctionGraph::JunctionGraph(){
    qDebug() <<"Initiate empty junction graph. ";
}

JunctionGraph::JunctionGraph(L1SkeletonGraph *g){
    init(g->getGraph());
}
void JunctionGraph::init(L1SkeletonGraph::SGraph  g){
    qDebug() <<" construct junction graph from skeleton graph... ";
	vIndex.clear();
	graph_.clear();
    typedef L1SkeletonGraph::vertex_iterator svertex_iterator;
    typedef L1SkeletonGraph::vertex_descriptor svertex_descriptor;
    typedef L1SkeletonGraph::un_edge_iterator sedge_iterator;
    typedef L1SkeletonGraph::edge_descriptor sedge_descriptor;
    typedef L1SkeletonGraph::adjacency_iterator sadjacency_iterator;
    sedge_iterator it,it_end;

    // add vertices to junction graph
    int i=0;
    for (tie(it,it_end) = edges(g );it!=it_end;++it ){
        svertex_descriptor s = source(*it,g),
                t = target(*it,g);


        vertex_descriptor st = add_vertex(graph_);
        graph_[st].id = i++;
        vertex_descriptor ts = add_vertex(graph_);
        graph_[ts].id = i++;

        vIndex[std::make_pair(g[s].id,g[t].id)] = graph_[st].id;
        vIndex[std::make_pair(g[t].id,g[s].id)] = graph_[ts].id ;
        // set vertex properties
        graph_[st].pairId = graph_[ts].id;
        graph_[ts].pairId = graph_[st].id;

        graph_[st].src = g[s].id;//i++;
        graph_[st].target = g[t].id;//i++;


        graph_[ts].src = g[t].id;//i++;
        graph_[ts].target = g[s].id;//i++;


        graph_[st].shape = g[*it].internalPts;
        graph_[ts].shape = g[*it].internalPts;
        std::reverse( graph_[ts].shape.begin(),graph_[ts].shape.end());
    }

    qDebug()<<"Added " << vIndex.size() << " vertices in junction graph. ";

    // add edges to junction graph
    svertex_iterator vit, vit_end;
    for (tie(vit,vit_end) = vertices(g); vit!=vit_end;++vit){
        // check v's adjacent vertices
        sadjacency_iterator ait, ait_end;
		
        for (tie(ait, ait_end) = adjacent_vertices(*vit,g);
			 ait!=ait_end;++ait){
			
            int i_av =  vIndex[std::make_pair(g[*ait].id, g[*vit].id)];
            vertex_descriptor av = vertex(i_av,graph_);

            sadjacency_iterator bit, bit_end;
            for (tie(bit, bit_end) = adjacent_vertices(*vit,g); bit!=bit_end;++bit){
                if (ait==bit) continue;

                int i_vb =  vIndex[std::make_pair(g[*vit].id, g[*bit].id)];
                vertex_descriptor vb = vertex(i_vb,graph_);
               // edge_descriptor e,bool exists;
                if ( ! edge(av,vb,graph_).second){
                    auto newEdge = add_edge(av,vb,graph_);
					graph_[newEdge.first].flow = 0;
                }
                //add_edge(i_av, i_vb,graph_);
                //graph_[e].Foo

                //tie(eNew, exist) = add_edge(*ait, *bit,);
                //incidentEdges.push_back(eNew);
            }
        }
    }
    qDebug() <<"Imported Junction graph of "<< num_vertices(graph_)
            <<" vertices, " << num_edges(graph_) << " edges";
}

void JunctionGraph::writeDotGraph(const char* fname,float height){
    std::ofstream out(fname, std::ofstream::out );
    write_graphviz(out, graph_,  VertexWriter(graph_ ,height ),EdgeWriter(graph_));//,GraphWriter(graph_) );
    qDebug() <<"write graph to file " << fname;
}

bool JunctionGraph::pairVertices(vertex_descriptor v1, vertex_descriptor v2){
    return graph_[v1].id == graph_[v2].pairId;
}
/*
bool JunctionGraph::nearlyAdjacent(vertex_descriptor v1, vertex_descriptor v2){
	bool nearlyAdj = false;
	if (graph_[v1].pairId = graph_[v2].id){
		return false;
	}
	if (graph_[v1].src 
	//float maxd = 2*Parameters::getInstance()->skeleton_cluster_radius  *min_edge_length_;
	// return false if L2(v1,v2) > maxL2Dist 
	// bfs from v1 to v2. break if bfs distance > maxL2Dist*3
	
	// if no path, return true
	// else compare shortest path distance and L2 distance
	// if L2*3 < shortest path, add edge 
	return nearlyAdj;
}*/
bool JunctionGraph::voteEdge(int s1, int t1, int s2, int t2, 
						     bool nearIncident, subtraj_descriptor sub){
	auto v1= vIndex.find(std::make_pair(s1,t1));
	auto v2 = vIndex.find(std::make_pair(s2,t2));

	vertex_descriptor d1 = vertex(v1->second ,graph_),
		d2 =vertex(v2->second, graph_);
	bool success=false;
	if (v1!= vIndex.end() && v2 != vIndex.end()){
		edge_descriptor e; bool exists;
		tie(e,exists) = edge(d1, d2, graph_);
		if (exists){
			graph_[e].flow+=1;
			graph_[e].sub_trajs.push_back(sub);
			success= true;
		}
			//			qDebug() <<"Note: Add new edge between jgraph vertex " <<
			//	graph_[d1].id << " and " << graph_[d2].id;
			// check wether d1 and d2 are *nearly incident*
			// i.e. check graph distance from d1 to d2, if longer than 2*l2(d1-d2)
			// add edged
		else if (nearIncident){
			tie(e,exists) = add_edge(d1,d2,graph_);
			graph_[e].flow = 1;
			graph_[e].sub_trajs.push_back(sub);
			success= true;			
		}
	}
	//qDebug() <<"vote edge (" << s1<<"," << t1 <<"," << s2 << "," <<t2 <<") success? "<< success;
	return success;
}
bool JunctionGraph::addEdge(int s1, int t1, int s2, int t2,
							subtraj_descriptor subtraj){

	auto v1 = vIndex.find(std::make_pair(s1,t1));
	auto v2 = vIndex.find(std::make_pair(s2,t2));
	if (v1->second == v2->second){
		qDebug()<<"Error: can not add edge "<< v1->second<<","<<v2->second<< endl;
		return false;
	}
	vertex_descriptor d1 = vertex(v1->second ,graph_),
		d2 =vertex(v2->second, graph_);
	edge_descriptor e; bool tmp;
	tie(e,tmp) = add_edge(d1,d2,graph_);
	graph_[e].flow = 1;
	graph_[e].sub_trajs.push_back(subtraj);
	//	addEdge();
	return true;
}

void JunctionGraph::getSubTrajClusters(SubTrajClusters &clusters){
	clusters.clear();
	edge_iterator it,it_end;
	for (tie(it,it_end)= edges(graph_);it!=it_end;++it){
		if (graph_[*it].flow>0){
			//		SubTrajCluster c = graph_[*it].trajs;
			int srcId = graph_[source(*it,graph_)].id;
			int targetId = graph_[target(*it,graph_)].id;
			std::array<int,2> key={srcId, targetId};
			clusters[key] =graph_[*it].sub_trajs;
		}
	}
}
/*
void JunctionGraph::filterEdgeByVotes(){
	float cutoff=getOptFlowThresh();
	qDebug()<<"=======================================";
	qDebug() <<"Original num of edges " << num_edges(graph_);
	qDebug() <<"Filter graph edge with flows lower than "<< cutoff;
	filter_edge_predicate<JGraph> filter(graph_, cutoff);
	typedef filtered_graph<JGraph, filter_edge_predicate<JGraph> > FGraph;

	FGraph fg(graph_, filter);
	qDebug() <<"Filtered graph has " << num_edges(fg) <<" edges.";
	
	JGraph gNew;
*/
	/*
	vertex_descriptor vit,vit_end;
   
	for (tie(vit,vit_end) =vertices(graph_);vit!=vit_end;++vit){
		vertex_descriptor v0 = *vit;//vertex(vi,graph_);
		vertex_descriptor v = add_vertex(gNew);
		gNew[v].id = graph_[v0].id;
		gNew[v].src = graph_[v0].src;
		gNew[v].target = graph_[v0].target;
		gNew[v].pairId = graph_[v0].pairId;
		gNew[v].shape = graph_[v0].shape;	
	}
	FGraph::edge_iterator eit, eit_end;
	for (tie(eit,eit_end) = edges(fg);eit!= eit_end;++eit){
		FGraph::vertex_descriptor s,t;
		s=source(*eit,fg);
		t = source(*eit,fg);
		edge_descriptor e;bool exists;
		tie(e,exists) = add_edge(
		}*/
/*
	copy_graph(fg, gNew);
	graph_ = gNew;
	//print_graph();
	// construct a filter graph based on flows > cutoff
	// rewrite JGraph (keep vIndex...)
	
}*/
void JunctionGraph::filterEdgeByVotes(){
	float cutoff = (minFlowSize_==-1)?getOptFlowThresh():minFlowSize_;
	
	qDebug()<<"=======================================";
	qDebug() <<"Original num of edges " << num_edges(graph_);
	qDebug() <<"Filter graph edge with flows lower than "<< cutoff;
	
	

	

	
	JGraph gNew;
	
	vertex_iterator vit,vit_end;
	//std::map<std::pair<int,int> , int > vIndexNew;
	
	for (tie(vit,vit_end) =vertices(graph_);vit!=vit_end;++vit){
		vertex_descriptor v0 = *vit;//vertex(vi,graph_);
		vertex_descriptor v = add_vertex(gNew);
		gNew[v].id = graph_[v0].id;
		gNew[v].src = graph_[v0].src;
		gNew[v].target = graph_[v0].target;
		gNew[v].pairId = graph_[v0].pairId;
		gNew[v].shape = std::vector<cv::Point2f>(graph_[v0].shape);	
	}
	edge_iterator eit, eit_end;
	
	for (tie(eit,eit_end) = edges(graph_);eit!= eit_end;++eit){
		if ( graph_[*eit].flow >= cutoff){
			vertex_descriptor s,t;
			s=source(*eit,graph_);
			t = target(*eit,graph_);
			edge_descriptor e1;bool exists;
			vertex_descriptor s1 = vertex(graph_[s].id, gNew),
				t1 = vertex(graph_[t].id,gNew);
			
			tie(e1,exists) = add_edge(s1,t1,gNew);
			gNew[e1].flow = graph_[*eit].flow;
			gNew[e1].sub_trajs.insert( gNew[e1].sub_trajs.end(),
									   graph_[*eit].sub_trajs.begin(),
									   graph_[*eit].sub_trajs.end());
			qDebug() <<"edge weight: " << gNew[e1].flow;

		}

	}							 
		//	copy_graph(fg, gNew);
		graph_ = gNew;
	qDebug() <<"Filtered graph has " << num_edges(graph_) <<" edges.";	
	//print_graph();
	// construct a filter graph based on flows > cutoff
	// rewrite JGraph (keep vIndex...)
	
}
float JunctionGraph::getOptFlowThresh(){
	// construct histogram
	edge_iterator eit,eit_end;
	
	cv::Mat flows(1,num_edges(graph_),CV_32S,cv::Scalar(0));
	int e=0;
	for (tie(eit,eit_end)=edges(graph_);eit!=eit_end;++eit){
		flows.at<int>(0,e++) = graph_[*eit].flow;
	}
	double maxFlow;
	double minFlow;
	minMaxLoc( flows,&minFlow, &maxFlow);
	qDebug() <<" max flow: " << maxFlow <<", " << "min flow: " <<minFlow ;
	cv::Mat hist(1,maxFlow+1,CV_32S, cv::Scalar(0));
	int sum=0; // sum of all edge weights
	for (int i=0;i< flows.cols;++i){
		if (flows.at<int>(i)>=1){
			hist.at<int>(flows.at<int>(i)) +=1;
			sum+= flows.at<int>(i);
		}

	}
	int total = flows.cols; //# of pixels (edges)
	//Background=2, Forground=1
	int w2 =0,w1 =0; float mean2=0,mean1=0, max=0.f,thresh1=0.f,thresh2=0.f;
	int between =0,sum2=0,sum1=0;
	for (int i=0;i<hist.cols;++i){
		w2+=hist.at<int>(i);
		if (w2==0) continue;
		w1 = total - w2;
		if (w1==0)
			break;
		sum2 += i *hist.at<int>(i);
		mean2 = sum2/float(w2);
		mean1 = (sum-sum2)/float(w1);
		between = w2*w1*(mean2-mean1)*(mean2-mean1);
		if (between >=max){
			thresh1 = i;
			if (between >max){
				thresh2 = i;
			}
			max = between;
		}
	}
	double thresh = (thresh1+thresh2)/2.0;
	//	double min_flow = maxFlow*Parameters::getInstance()->min_flow_ratio;
	double min_flow = Parameters::getInstance()->skeleton_min_cluster_size;
	qDebug()<<"\nMax variance cluster size threshold: " << thresh
			<<"; actual threshold: "<< min_flow;
	//return std::min(min_flow, thresh);
	return min_flow;
}

/*JunctionGraph::vertex_descriptor JunctionGraph::getVertexById(int i){
    vertex_iterator it,it_end;
    for (tie(it,it_end) = vertices(graph_);it!=it_end;++it){
        if (graph_[*it].id == i){

            return *it;
        }
    }
    return vertex_descriptor();
}*/
