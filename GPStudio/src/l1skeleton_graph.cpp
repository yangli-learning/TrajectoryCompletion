#include <QObject>
#include "l1skeleton_graph.h"
#include <QDebug>
#include <string>
#include <sstream>
#include <fstream>
#include <boost/graph/graphviz.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>

#include <opencv2/opencv.hpp>
L1SkeletonGraph::L1SkeletonGraph(BranchList branches,
                                 BridgeMap bridges,
                                 std::vector<SkeletonPoint> bridgePoints,float height){
    nCompactVertices=0;
    height_ = height;
    int nEdges = addEdgesFromBranches(branches);
	addEdgesFromBridges(branches,bridges,bridgePoints,nEdges);
    assignVertexType();

    removeInternalNodes();    
    writeDotGraph("graph1.dot");

}
void L1SkeletonGraph::mergeSimilarBranches(float minLen, float minAngle){
	// sort vertices by degree, for each vertex v, compute angle between branches
	// remove degree dangling edges if its angle to an incident branch
	// is too small

	qDebug() << "[SGraph] Simplify graph (apply branch angle filter)...";
	std::vector<vertex_descriptor> vs;
	vertex_iterator vit,vit_end;
	float similarDist = minLen;//Parameters::getInstance()->skeleton_min_edge_len;
	for (tie(vit,vit_end)=vertices(graph_);vit!=vit_end;++vit){
		vs.push_back(*vit);
	}
	while (vs.size()>1){
		auto maxit = std::max_element(vs.begin(),vs.end(),
		  [this](const vertex_descriptor &v1, const vertex_descriptor &v2){
						   return out_degree(v1,this->graph_)< out_degree(v2,this->graph_);});
		if (out_degree(*maxit,graph_)<=1)
			break;
		//		qDebug() <<"max degree element: " << graph_[*maxit].id <<" (" << 
		//	out_degree(*maxit,graph_)
		//	 << ") !\n";

		// check adjacent vertices
		adjacency_iterator ait,ait_end;
		std::vector<vertex_descriptor> adjs;
		for (tie(ait,ait_end)=adjacent_vertices(*maxit,graph_);ait!=ait_end;++ait){
			adjs.push_back(*ait);
		}
		// compute matrix of angles
		cv::Mat branchAngles(adjs.size(),adjs.size(),CV_32F, cv::Scalar(1000.f));
		for (size_t i=0;i< adjs.size();++i){
			if (out_degree(adjs[i],graph_)>1)
				continue;
			for (size_t j=0;j<adjs.size();++j){
				if (i!=j){
					float dist = cv::norm(graph_[adjs[i]].p.pos
										  - graph_[adjs[j]].p.pos);
					if (dist < similarDist){
						branchAngles.at<float>(i,j)
							= computeBranchAngle(*maxit,adjs[i],adjs[j]);
					}
				}
			}
		}
		// if lowest angle 
		double minV=0,maxV; cv::Point minL,maxL;
		
		while (float(minV) < minAngle){
			cv::minMaxLoc(branchAngles,&minV,&maxV,&minL,&maxL);
			//		qDebug() << "min val: " << minV <<", max val " << maxV;
		
			if (minV < minAngle){
				vertex_descriptor va = adjs[minL.y];
				remove_edge(*maxit,va,graph_);
				qDebug()<<"delete edge " << graph_[*maxit].id<<"--"<<graph_[va].id
						<< " " << minL.x<<"," <<minL.y << " angle=" << minV
						<<" minAngle=" << minAngle;
				branchAngles.row(minL.y)=cv::Scalar(1000.f);
				branchAngles.col(minL.y) = cv::Scalar(1000.f);
				//				branchAngles.at<float>(minL.y,minL.x)=1000.f;
				//				qDebug()<<" minV=" <<float(branchAngles.at<float>(minL.x,min);
				//break;	<<
			}
		}
		//
		vs.erase(maxit);
	}
	assignVertexType();
	updateEdgeIndex();
	
}
float L1SkeletonGraph::computeBranchAngle(vertex_descriptor v, vertex_descriptor u,
										  vertex_descriptor w){
	cv::Point2f v1 = graph_[u].p.pos - graph_[v].p.pos,v1n,
		v2=graph_[w].p.pos-graph_[v].p.pos,v2n;
	Common::normalize(v1,v1n); Common::normalize(v2,v2n);

	return acos(double(v1n.dot(v2n)))*180.f/(M_PI);
	
}
int L1SkeletonGraph::addEdgesFromBranches(const BranchList &branches){

    int edgeCount=0;
    for (size_t i=0;i< branches.size();++i){
        //iterate over all branch, add edges along each
        const std::vector<Candidate> &branch = branches[i];

        // insert first point on branch to vertex map
        for (size_t j=0; j<branch.size()-1;++j){
            vertex_descriptor v1,v2;
            v1= getVertex(branch[j]);
            v2 = getVertex(branch[j+1]);

            if (!edge(v1,v2,graph_).second){
                auto e = add_edge(v1,v2,graph_);
                graph_[e.first].id = edgeCount++;
                graph_[e.first].internalPts.resize(2);
                graph_[e.first].internalPts[0] = branch[j].pos;
                graph_[e.first].internalPts[1] = branch[j+1].pos;
            }else{
                qDebug()<<"Error: duplicated edge";
            }
            v1=v2;

        }
    }
    qDebug()<<"Adding " << edgeCount << " edges from edges..."  ;
    return edgeCount;
}
bool L1SkeletonGraph::vIndexExists(int skeletonId, std::map<int, int> vertexIndex){
    return (vertexIndex.find(skeletonId)!=vertexIndex.end());
}

bool L1SkeletonGraph::vIndexExists(int skeletonId){
    return (vIndex_.find(skeletonId)!=vIndex_.end());
}


void L1SkeletonGraph::writeDotGraph(const char* fname){
    std::ofstream out(fname, std::ofstream::out );
    write_graphviz(out, graph_,  VertexWriter(graph_,height_),EdgeWriter(graph_) );
    qDebug() <<"write graph to file " << fname;
}
void L1SkeletonGraph::updateEdgeIndex(){
	int id=0;
	un_edge_iterator it,it_end;
    for (tie(it,it_end) = edges(graph_);it!= it_end; ++it )
		graph_[*it].id = id++;
}
void  L1SkeletonGraph::updateEdgeGeometry(int edgeId,
                                          std::vector<cv::Point2f> pts){
										  //,  bool smooth){
     edge_descriptor e;bool exists;
     tie(e,exists)= getEdgeById(edgeId);
     if (exists){
		 updateEdgeGeometry(e,pts);//smooth);
	 }
}
void  L1SkeletonGraph::updateEdgeGeometry(edge_descriptor e,
                                          std::vector<cv::Point2f> pts){
										  //	  bool smooth){
        //qDebug() << "edge id = "<<graph_[e].id << " edgeId " << edgeId ;

		//float l = 0;
	bool smooth = Parameters::getInstance()->skeleton_edge_smoothing;
	 if (pts.size()>=3 && smooth){
		 std::vector<cv::Point2f> smooth_pts;
		 smooth_pts.push_back(pts.front());
		 for (size_t i=0;i<pts.size()-2;++i){
			 smooth_pts.push_back(0.5*pts[i]+0.5*pts[i+2]);
			 //	 l+= cv::norm(smooth_pts[i]-smooth_pts[i+1]);
		 }
		 smooth_pts.push_back(pts.back());

		 size_t n = smooth_pts.size();
		 // l+= cv::norm(smooth_pts[n-2]-smooth_pts[n-1]);
		 graph_[e].internalPts = smooth_pts;
	 }else{
		 graph_[e].internalPts =pts;
	 }
	 graph_[e].length = Common::polylineLength( graph_[e].internalPts);		
}
// get Nth edge (N is different edges[N].id !!!)
std::pair<L1SkeletonGraph::edge_descriptor,bool> L1SkeletonGraph::getEdgeById(size_t edgeId){
    un_edge_iterator it,it_end;
    tie(it,it_end) = edges(graph_);
    edge_descriptor e; int ct=0;
    bool edgeExists = false;
    while( it!=it_end){
        if ( (ct++) == edgeId){
            e = *it;
            edgeExists = true;
            break;
        }
        ++it;
    }
    return std::make_pair(e, edgeExists);
                //edgeExists;
}


bool L1SkeletonGraph::incidentEdges(edge_descriptor e1, edge_descriptor e2){

    return source(e1,graph_) == source(e2,graph_) ||source(e1,graph_) == target(e2,graph_)
            || target(e1,graph_) == source(e2,graph_) || target(e1,graph_)==target(e2,graph_);
}

bool L1SkeletonGraph::incidentEdges(int eid1, int eid2){
    // unsafe
    return incidentEdges(getEdgeById(eid1).first, getEdgeById(eid2).first);
}

L1SkeletonGraph::vertex_descriptor L1SkeletonGraph::getVertex(SkeletonPoint p){
    if ( ! vIndexExists(p.id, vIndex_)){
         int vId = ( int)num_vertices(graph_ );
         vIndex_[p.id] = vId;
         vertex_descriptor v = add_vertex(graph_);
         graph_[v].id = p.id;
         graph_[v].id = p.id; //unique id
         graph_[v].p = p; //skeleton point
         graph_[v].type = SkeletonGraphVertexType::UNKNOWN;

    }
     return vertex(vIndex_[p.id],graph_);
}
L1SkeletonGraph::vertex_descriptor L1SkeletonGraph::getVertex(SkeletonPoint p,
         SGraph &g,   std::map<int, int> &vertexIndex   ){
    if ( ! vIndexExists(p.id, vertexIndex)){
         int vId = ( int)num_vertices(g );
         vertexIndex[p.id] = vId;
         vertex_descriptor v = add_vertex(g);
         g[v].id = p.id;
         g[v].id = p.id; //unique id
         g[v].p = p; //skeleton point
         g[v].type = SkeletonGraphVertexType::UNKNOWN;

    }
     return vertex(vertexIndex[p.id],g );
}

int L1SkeletonGraph::addEdgesFromBridges(BranchList branches,
                                         BridgeMap bridges,
                                         std::vector<SkeletonPoint> bridgePoints,
                                         int edgeCount){
    auto iter= bridges.begin();
    int bridge_i=0;

    while (iter!=bridges.end()) {

        int bridgePtId = iter->first;
        auto bridgePt = find_if(bridgePoints.begin(),bridgePoints.end(),
                                [&bridgePtId](const SkeletonPoint &v)
                                {return v.id==bridgePtId;});
        if (bridgePt == bridgePoints.end()) {
            iter++;
            continue;
        }
        // set bridgePt in vertices_ if it is not yet defined
        vertex_descriptor bridge_ds = getVertex(*bridgePt);
        bridge_i++;

        float branchCosThresh = 0.1;
        //check each connected branch_i of bridgePt
        std::set<int> &branchIds = iter->second;
        for (auto vit= branchIds.begin(); vit != branchIds.end();++vit){

            std::vector<Candidate> &branch_i =  branches[*vit];
            // if bridgePt is not in branch_i, connect it to
            // one of the endpoint of branch_i
            auto bit = L1Skeleton::searchCandidates(branch_i,bridgePt->id);
            if (bit == branch_i.end()){
                // find closes endpoint of branch_i to bridgePt
                SkeletonPoint e1 = branch_i.front(),
                    e2 = branch_i.back();
                SkeletonPoint target;
                cv::Point2f vec1, vec2;
                bool isHead;
                if  (cv::norm(e1.pos-bridgePt->pos)
                     < cv::norm(e2.pos-bridgePt->pos)){
                    // bridge is connect to head
                    target = e1;
                    vec1 = branch_i[1].pos - e1.pos;
                    isHead = false;
                }else{
                    // bridge is connected to tail
                    target = e2;
                    int branchLen = branch_i.size();
                    vec1 =e1.pos - branch_i[branchLen-2].pos ;
                    isHead = true;
                }
                vec2 = bridgePt->pos - target.pos;
                float cosTest =  cosAngleBetween(vec1,vec2) ;
                if ( true || cosTest> branchCosThresh){

                    vertex_descriptor target_ds = getVertex(target);

                    //auto e = add_edge(bridge_ds, target_ds, graph_);
                    if (! edge(bridge_ds, target_ds,graph_).second){
                        auto e = add_edge(bridge_ds, target_ds, graph_);
                        graph_[e.first].id = edgeCount++;
                        graph_[e.first].internalPts.resize(2);
                        graph_[e.first].internalPts[0] = (isHead)?target.pos:bridgePt->pos;
                        graph_[e.first].internalPts[1] = (isHead)?bridgePt->pos:target.pos;
                    }
                }
            }
        }
        iter++;
        bridge_i ++;
    }
    qDebug() << " finished processing bridge points " ;
    return edgeCount;
}
float L1SkeletonGraph::cosAngleBetween(const cv::Point2f &v1, const cv::Point2f &v2){
    float len1 = sqrt(v1.x * v1.x + v1.y * v1.y);
    float len2 = sqrt(v2.x * v2.x + v2.y * v2.y);

    float dot = v1.x * v2.x + v1.y * v2.y;
    if (len1*len2 ==0){
        return 0;
    }else{
        return dot / (len1 * len2);
    }
}

//! should not reference p.pos for vertex position: always use pts
void L1SkeletonGraph::moveVertex( vertex_descriptor v,  cv::Point2f newLoc){
    //update all incident  edges and their internalPts
	cv::Point2f oldLoc = graph_[v].p.pos;
    graph_[v].p.pos = newLoc;
    adjacency_iterator ait,ait_end;
    for (tie(ait,ait_end)=adjacent_vertices(v,graph_); ait!= ait_end; ++ait){
        // replace loc in edge points
         edge_descriptor e; bool exists;
         tie(e,exists)= edge(v,*ait,graph_);
         if (!exists){
             qDebug() << "incident edge doesn't exist for neighbor vertex "
					  << graph_[*ait].id;
             continue;
         }
         std::vector<cv::Point2f> &pts = graph_[e].internalPts;
         //cv::norm(newLoc-  pts[pts.size()-1]  )  < cv::norm(newLoc-pts[0]) ){
         // pts[pts.size()-1] =newLoc;
         //if (  source(e,graph_) ==v){
         if (cv::norm(oldLoc -  pts[pts.size()-1]  ) > cv::norm(oldLoc-pts[0])  ){
                pts[0]   = newLoc;
         }else{
           pts[pts.size()-1] =newLoc;
         }
    }
}
bool L1SkeletonGraph::isValidEdge(edge_descriptor e){
     vertex_descriptor s =source(e, graph_), t=target(e,graph_);
     return edge(s, t,graph_).second;
}
int L1SkeletonGraph::getCurrentEdgeId(edge_descriptor e){
     vertex_descriptor s =source(e, graph_), t=target(e,graph_);
     if (! edge(s, t,graph_).second){
         return -1;

     }
     int id=0;
     un_edge_iterator it,it_end;
     tie(it,it_end) = edges(graph_);

     while( it!=it_end){
         vertex_descriptor i_s = source(*it,graph_),
                 i_t = target(*it ,graph_);
        if( (s==i_s && t==i_t )|| (s==i_t && t==i_s) ){
             break;
         }
         ++id;
     }
    return id;
}

void L1SkeletonGraph::splitEdgeAt(edge_descriptor e, vertex_descriptor v, int segPos ){
    vertex_descriptor s =source(e, graph_), t=target(e,graph_);
    qDebug() << " split edge " << graph_[s].id <<"--" << graph_[t].id <<" by intersection ";
    // check edge exists
    if (! edge(s, t,graph_).second){
        qDebug() <<"!!!!!\nWarning: neighbor edge does not exist, cannot update graph ";
        return  ;
    }

    // check orientation of internal pts for edge e.
    std::vector<cv::Point2f> pts = graph_[e].internalPts;
    if (cv::norm(pts[0]-graph_[s].p.pos )  > cv::norm(pts.back()-graph_[s].p.pos )){
        qDebug() <<"!!!!!!!!!!!!!!!!!!\nError: edge orientation is reversed! ";
        //assert(false);
        std::reverse(pts.begin(),pts.end());
    }

    // add new edges  nEdgeS =(e.src,  v), nEdgeD = ( v, e.target)
    int edgeCt = num_edges(graph_);
    edge_descriptor es,et; bool sucS, sucT;
    tie(es, sucS)=  add_edge( s, v, graph_),
    tie(et, sucT) = add_edge( v, t, graph_);
    // set nEdgeS.internalPts = { nPts[0:nPos], loc}
    graph_[es].id =edgeCt;
    graph_[es].internalPts.insert(graph_[es].internalPts.end(),pts.begin(),pts.begin()+ segPos+1);
    graph_[es].internalPts.push_back(graph_[v].p.pos);

    // set nEdgeD.internalPts = {loc, nPts[nPos+1:N]}
    graph_[et].id = edgeCt+1;
    graph_[et].internalPts.push_back(graph_[v].p.pos);
    graph_[et].internalPts.insert(graph_[et].internalPts.end(), pts.begin()+segPos+1, pts.end());

    // remove nEdge from G. <---This is last!
    qDebug() <<"Removing edge " <<graph_[e].id;
    remove_edge( e, graph_);
    qDebug() <<"Remaining # edges = " << num_edges(graph_);
    return  ;
}
void L1SkeletonGraph::splitEdgeAtIntersection(edge_descriptor e, vertex_descriptor v,
                                              int segPos,
											  cv::Point2f intersectPt,
											  bool keepFirstHalf, bool keepLastHalf ){
    vertex_descriptor s =source(e, graph_), t=target(e,graph_);
    std::vector<cv::Point2f> pts = graph_[e].internalPts;
    if (cv::norm(pts[0]-graph_[s].p.pos )  > cv::norm(pts.back()-graph_[s].p.pos )){
        qDebug() <<"!!!!!!!!!!!!!!!!!!\nError: edge orientation is reversed! ";
             // assert(false);
        std::reverse(pts.begin(),pts.end());
    }
     int edgeCt = num_edges(graph_);
     edge_descriptor newEdge; bool success;
    if (keepFirstHalf){
        qDebug() <<"add edge s-->intersection";
        if (edge(s,v,graph_).second){
			qDebug() <<"edge already exists from "<<graph_[s].id <<" to "
					 << graph_[v].id;
		}else{
			moveVertex(v,intersectPt);
			tie(newEdge, success)=  add_edge( s, v, graph_);
			graph_[newEdge].id = edgeCt++;
			std::vector<cv::Point2f> newPts(pts.begin(),pts.begin()+ segPos+1);
			newPts.push_back(graph_[v].p.pos);
			//graph_[newEdge].internalPts.insert(graph_[newEdge].internalPts.end(),
			//								   pts.begin(),pts.begin()+ segPos+1);
			//graph_[newEdge].internalPts .push_back(graph_[v].p.pos);
			updateEdgeGeometry(newEdge, newPts);
		}


    }
    if (keepLastHalf){
        qDebug() << " add edge intersection --> t";
        if (edge( v,t,graph_).second){
			qDebug() <<"edge already exists from "<<graph_[v].id <<" to "
					 << graph_[t].id;
		}else{
			moveVertex(v,intersectPt);
			tie(newEdge, success) = add_edge( v, t, graph_);
			graph_[newEdge].id = edgeCt++;
			std::vector<cv::Point2f> newPts;
			newPts.push_back(graph_[v].p.pos);
			newPts.insert(newPts.end(),
						  pts.begin()+segPos+1, pts.end());
			//graph_[newEdge].internalPts.push_back(graph_[v].p.pos);
			//graph_[newEdge].internalPts.insert(graph_[newEdge].internalPts.end(),
            //                               pts.begin()+segPos+1, pts.end());
			updateEdgeGeometry(newEdge,newPts);
		}

    }

    // remove nEdge from G. <---This is last!
    qDebug() <<"Removing edge " <<graph_[e].id;
    remove_edge( e, graph_);
    qDebug() <<"Remaining # edges = " << num_edges(graph_);
    return  ;
}

void L1SkeletonGraph::printGraph(){
    qDebug() << "===== Skeleton Graph ========";

    //const VertexIndexMap vIndex = get( propmapIndex, graph_);

    un_edge_iterator it,it_end;
    vertex_descriptor u,v;
    for (tie(it,it_end) = edges(graph_);it!= it_end; ++it ){
        std::ostringstream oss;
        u=source(*it,graph_);
        v = target(*it,graph_);
        oss << "edge " << graph_[*it].id << "(" << graph_[u].id <<", ["
            <<   graph_[u].p.pos.x << "," <<  graph_[u].p.pos.y <<"]: "
              << vertexTypeToString(graph_[u].type) <<" <->"<< graph_[v].id
             <<", ["
           <<   graph_[v].p.pos.x << "," <<  graph_[v].p.pos.y <<"]: "
             << vertexTypeToString(graph_[v].type)  <<" len=" << graph_[*it].internalPts.size()<<" | ";

        for(int i=0;i<graph_[*it].internalPts.size();++i){
            oss  << "(" << graph_[*it].internalPts[i].x << ", " << graph_[*it].internalPts[i].y
                    <<") ";
        }
        oss <<" len=" <<graph_[*it].length;
        qDebug() << oss.str().c_str() <<"\n";
    }
}
void L1SkeletonGraph::writeSkeletonImage(const char* fname,int height, int width){
    qDebug() << "Writing skeleton graph to image " << fname <<"...";
	int scale=5;
	cv::Mat img( (height+1)*scale,(width+1)*scale,CV_8UC3,cv::Scalar(255,255,255));
	

    //const VertexIndexMap vIndex = get( propmapIndex, graph_);

    un_edge_iterator it,it_end;
    vertex_descriptor u,v;
	int i=0;
    for (tie(it,it_end) = edges(graph_);it!= it_end; ++it ){
        std::ostringstream oss;
        u=source(*it,graph_);
        v = target(*it,graph_);
		//draw polyline
		std::vector<cv::Point> polyline;
        for(int i=0;i<graph_[*it].internalPts.size();++i){
			polyline.push_back(cv::Point( int( scale * graph_[*it].internalPts[i].x),
									  int(scale * graph_[*it].internalPts[i].y)));
        }
		const cv::Point *pts = (const cv::Point*) cv::Mat(polyline).data;
		int npts = polyline.size();
		int colorIdx = (i++) % Common::palette_size;
		cv::Scalar color(Common::palette[colorIdx][0],
						 Common::palette[colorIdx][1],
						 Common::palette[colorIdx][2] );
		cv::polylines(img, &pts, &npts, 1,false,color,2);
	  
    }

	
    vertex_iterator vit,vit_end;
    for (boost::tie(vit,vit_end) = vertices(graph_); vit != vit_end; ++vit){

		cv::Point c(int(scale*graph_[*vit].p.pos.x), int(scale*graph_[*vit].p.pos.y));
		
		cv::circle(img, c, scale, cv::Scalar(100,100,100),-1);
	}
	try{
		cv::imwrite(fname, img);
	}catch(std::runtime_error &ex){
        qDebug() << "Exception converting image to PNG format: " << ex.what();
	}
	//write line to file
}
void L1SkeletonGraph::updateVertexPos(){
    un_edge_iterator it,it_end;
    vertex_descriptor u,v;
    for (tie(it,it_end) = edges(graph_);it!= it_end; ++it ){
        std::ostringstream oss;
        u=source(*it,graph_);
        v = target(*it,graph_);
		//        qDebug() <<"update vertex pos for  edge "<< graph_[*it].id << " of len "
		//      <<  graph_[*it].internalPts.size();
        if (  graph_[*it].internalPts.size()>=1 ){
            graph_[u].p.pos = graph_[*it].internalPts.front();
            graph_[v].p.pos = graph_[*it].internalPts.back();
        }
    }
}

void L1SkeletonGraph::assignVertexType(){

	// int unknowns=0,terminals=0,internals=0, intersections=0;

    typename graph_traits<SGraph>::vertex_iterator it,it_end;
    for (boost::tie(it,it_end) = vertices(graph_); it != it_end; ++it){
        SGraph::vertex_descriptor u =  *it ;
        switch ( degree(u, graph_)){
        case 0:
            graph_[u].type = SkeletonGraphVertexType::UNKNOWN;
            //unknowns++;
            break;
        case 1:
            graph_[u].type = SkeletonGraphVertexType::TERMINAL;
            //terminals++;
            break;
        case 2:
            if (isTurningPoint(u)){
                 graph_[u].type  =SkeletonGraphVertexType::INTERSECTION;
				 //intersections++;
            }else{
                 graph_[u].type =SkeletonGraphVertexType::INTERNAL;
				 // internals++;
            }
            break;
        default:
             graph_[u].type = SkeletonGraphVertexType::INTERSECTION;
			 //intersections++;
            break;
        }
    }
    // if node is internal, check the angle between 2 edges,
	// then mark as internal if
    // angle close to 90 degrees .

    //qDebug() << "unknowns: " << unknowns <<" terminals: " << terminals
    //         << " internal: "
    //         << internals << " intersection: " << intersections;
	//    nCompactVertices = terminals+intersections;
}
bool L1SkeletonGraph::isTurningPoint(SGraph::vertex_descriptor v){
    //return false;
    if ( degree(v,graph_)!=2){
        return false;
    }
    cv::Point2f p = graph_[v].p.pos;
    adjacency_iterator adj = adjacent_vertices(v ,graph_).first;
    SGraph::vertex_descriptor n1 = *(adj++), n2 = *(adj);


    cv::Point2f v1 = graph_[ n1].p.pos - p;
    cv::Point2f v2 = graph_[ n2].p.pos - p;
	/*    if ( graph_[n1].id == graph_[n2].id){
        qDebug() <<"error: first and second adjacent vertices are the same! ";
    }else{
        qDebug() << " neighbor 1:" <<  graph_[ n1].id
                 <<" neighbor 2:" << graph_[ n2].id << " p: "<< graph_[v].id << v;
        qDebug() << "angle between v1,v2 " << L1Skeleton::angleBetween(v1,v2)
                 << " cos " <<cos(L1Skeleton::angleBetween(v1,v2)) ;
				 }*/
    if (cosAngleBetween(v1,v2) < -0.5){
        return false;
    }
    return true;

}
void L1SkeletonGraph::traceSimplePath(vertex_descriptor vertex, edge_descriptor edge,
                                     std::vector<int> &edgePoints){
    graph_[edge].visited = true;
    edgePoints.push_back(graph_[vertex].id); //add to edge points in compact graph
    vertex_descriptor  s = source( edge,graph_) , t =target(edge,graph_); // id in the current gaph
    vertex_descriptor next = (graph_[s].id==graph_[vertex].id)? t : s ;

    if (graph_[next].type == SkeletonGraphVertexType::INTERSECTION ||
        graph_[next].type == SkeletonGraphVertexType::TERMINAL){

        edgePoints.push_back(graph_[next].id); // edge points are indexed based on the current graph

    }else if (graph_[next].type == SkeletonGraphVertexType::INTERNAL){
        edge_iterator it, it_end;
        for (tie(it,it_end) = out_edges( next, graph_);it!=it_end;++it){
            if ( graph_[*it].visited){
                continue;
            }
            traceSimplePath(next,*it,edgePoints);
            break;
        }
    }
}
void L1SkeletonGraph::initializeEdgeStatus(){

    un_edge_iterator it,it_end;
    for (tie(it,it_end) = edges(graph_);it!= it_end; ++it ){
        graph_[*it].visited = false;
    }
}
/* special cases:
 * 1. if tracing results in a cycle, break the cycle by removing the last edge
*/
void L1SkeletonGraph::removeInternalNodes( ){
	//update: do not remove internal node if head = tail
    // mark all edges unvisited
	qDebug()<<" \n\n[SGraph] Simplify graph (remove internal nodes)...";
    initializeEdgeStatus();

    // create temporary graph
    SGraph gNew;
    std::map<int,int> vIndexNew;
    int edgeCount=0;
    vertex_iterator it,it_end;
    std::vector<int> nonInternals; // graph ids of noninternal vertices

    // add all non-internal nodes to new graph
    for (tie(it,it_end) = vertices(graph_);it!=it_end;++it){
        if (*it==0)
            continue;
        int graphId =  vIndex_[graph_[*it].id];
        if (graph_[*it].type==SkeletonGraphVertexType::INTERSECTION ||
            graph_[*it].type==SkeletonGraphVertexType::TERMINAL){
            nonInternals.push_back(graphId);
            getVertex(graph_[*it].p, gNew, vIndexNew);
        }

     }
     // trace edges and store in new graph
     for (size_t i=0; i<nonInternals.size();++i){
        vertex_descriptor v= vertex(nonInternals[i], graph_);

        adjacency_iterator it,it_end;
        for (tie(it,it_end) = adjacent_vertices(v,graph_); it!=it_end; ++it ){
            edge_descriptor e = edge(v,*it,graph_).first;
            if (graph_[e].visited){
                continue;
            }
            std::vector<int> edgePoints;
            traceSimplePath(v, e , edgePoints ); 
            vertex_descriptor v1,v2;
            v1=getVertex(graph_[vIndex_[edgePoints.front()] ].p, gNew, vIndexNew);
            v2=getVertex(graph_[vIndex_[edgePoints.back()] ].p, gNew, vIndexNew);
            if (edgePoints.size()>=2){
                std::vector<cv::Point2f> plist = mergeInternalPts( graph_[v].p.pos , edgePoints);
                int id1 = vIndexNew[edgePoints.front()], id2 = vIndexNew[edgePoints.back()];
                //qDebug() <<"Adding edge " << id1 <<" ("<<edgePoints.front() <<") - "
                //           << id2 << " ("<< edgePoints.back() <<")";
                edge_descriptor e0; bool edgeExists;
                tie(e0, edgeExists) =  edge(v1, v2,gNew);
                if ( edgeExists){
                    qDebug() << "Warning: found duplicate edge. new length: "
                             << plist.size() <<" old length: " << gNew[e0].internalPts.size();
                    if ( gNew[e0].internalPts.size()  > plist.size()){
                        qDebug() <<"update to shorter internal points ";
                        auto src= source(e0,gNew);
                        if (gNew[src].id!=gNew[v1].id){
                            std::reverse(plist.begin(),plist.end());
                        }
                        gNew[e0].internalPts =   plist;

                    }
					// qDebug() <<"skip adding new edge";
                    continue;
                }
                // otherwise, add new edge
                 edge_descriptor eNew  =add_edge(id1,id2,gNew).first;
                gNew[eNew].id =edgeCount++;
                //gNew[eNew].internalPts = plist;
				updateEdgeGeometry(eNew,plist);

            }
        }
    }
    graph_= gNew;
    vIndex_= vIndexNew;
    removeDanglingEdges( );
	assignVertexType();
	updateEdgeIndex();
	computeEdgeWeight();

}
void L1SkeletonGraph::computeEdgeWeight(){ 
    un_edge_iterator it,it_end;
    for(tie(it,it_end)=edges(graph_);it!=it_end;++it){

            std::vector<cv::Point2f> &poly = graph_[*it].internalPts;
            graph_[*it].length = Common::polylineLength(poly); 
       }
}

std::vector<cv::Point2f> L1SkeletonGraph::mergeInternalPts(cv::Point2f start,
                                                           std::vector<int> edgePoints){
    std::vector<cv::Point2f> pts;
    if (edgePoints.size()==0)
        return pts;
    std::pair<edge_descriptor,bool> e; vertex_descriptor v1,v2;
    cv::Point2f firstPt = start;
    for(size_t i=0;i<edgePoints.size()-1;++i){

        e = edge(vIndex_[edgePoints[i]], vIndex_[edgePoints[i+1]],graph_ );
        if (e.second){

            std::vector<cv::Point2f> internalPts(graph_[e.first].internalPts);
            if (firstPt !=  graph_[e.first].internalPts.front()){
                std::reverse(internalPts.begin(), internalPts.end());
            }
            firstPt= internalPts.back();

            if (i==0){
                pts.insert(pts.end(),internalPts.begin(),internalPts.end());
            }else if (graph_[e.first].internalPts.size()==2){
                pts.push_back( internalPts.back());
            }else{
                pts.insert(pts.end(),internalPts.begin()+1,internalPts.end());
            }
        }
    }
    return pts;
}

void L1SkeletonGraph::removePath(std::vector<vertex_descriptor> &edgePoints){

    for (size_t i=0; i<edgePoints.size()-1;++i){
        edge_descriptor e;bool edge_exists;
        tie(e,edge_exists)= edge(edgePoints[i],edgePoints[i+1],graph_);
        if (edge_exists){
            remove_edge(e,graph_);
        }else{
            qDebug() <<"Can not remove edge between points"<< i <<" and "
                    <<i+1 << " already removed!";
        }
    }
}
std::vector<cv::Point2f> L1SkeletonGraph::getVertexPositions(){
    std::vector<cv::Point2f> pts;
    vertex_iterator v,v_end;
    for (tie(v,v_end)=  vertices(graph_); v!=v_end;++v){
        if (degree(*v,graph_) > 0){
            pts.push_back(graph_[*v].p.pos);
        }
    }
    return pts;
}

std::vector<std::vector<cv::Point2f> > L1SkeletonGraph::getEdgeSegments(){
    std::vector<std::vector<cv::Point2f> > edgeSet;
    un_edge_iterator it,it_end;
    for (tie(it,it_end)=edges(graph_); it!=it_end; ++it){
        if (graph_[*it].internalPts.size() >=2){
            edgeSet.push_back( graph_[*it].internalPts);
        }else{
            std::vector<cv::Point2f> pts;
            pts.push_back(graph_[source(*it,graph_)].p.pos);
            pts.push_back(graph_[target(*it,graph_)].p.pos);
            edgeSet.push_back(pts );
        }
    }
    return edgeSet;
}


void L1SkeletonGraph::removeDanglingEdges( ){
    //std::vector<edge_iterator> edgeSet;
   //std::set<int> removedEdges;

    edge_predicate<SGraph> filter(graph_);
    typedef filtered_graph<SGraph, edge_predicate<SGraph> > FGraph;
    FGraph fg(graph_, filter);
    if (num_edges(fg) == num_edges(graph_))
        return; // no change to edges
    SGraph gNew(num_vertices(graph_));
    int i=0; vertex_iterator vit,vit_end;
    for (tie(vit,vit_end) =vertices(gNew); vit!=vit_end;++vit){
        vertex_descriptor v =vertex(i++ ,graph_);
        gNew[*vit].id= graph_[v].id;
        gNew[*vit].type = graph_[v].type;
        gNew[*vit].p =graph_[v].p;

    }
    FGraph::edge_iterator it,it_end;
    for (tie(it,it_end) = edges(fg); it!= it_end; ++it){
        FGraph::vertex_descriptor s,t;
        s =source(*it,fg);
        t = target(*it,fg);
        edge_descriptor e;bool exists;
        assert( vIndex_[fg[s].id] !=  vIndex_[fg[t].id] );
        tie(e,exists) = add_edge(vIndex_[fg[s].id], vIndex_[fg[t].id],gNew);
        gNew[e].internalPts = fg[*it].internalPts;
        gNew[e].id = fg[*it].id;
        gNew[e].visited = false;

    }
    //vIndex_ remains the same since vertices did not change
    graph_=gNew;
    qDebug() << "[SGraph] remove dangling edges";
}

std::pair<L1SkeletonGraph::edge_descriptor,bool>
L1SkeletonGraph::addEdge(vertex_descriptor src, vertex_descriptor dst,
						 std::vector<cv::Point2f> pts){
	bool success =false;
	edge_descriptor e;
    if (graph_[src].id!= graph_[dst].id && !edge(src,dst,graph_).second){
		e = add_edge(src,dst,graph_).first;
		graph_[e].id = num_edges(graph_);
		graph_[e].internalPts = pts;
		graph_[e].visited = false;
		graph_[e].length = Common::polylineLength(pts);
		success = true;
	}
	return std::make_pair(e, success);
}
void L1SkeletonGraph::removeEdge( vertex_descriptor src, vertex_descriptor dst){
    if (!edge(src,dst,graph_).second){
        qDebug() <<"\nError: can not remove edge " << graph_[src].id<<","<< graph_[dst].id
                <<": edge not found!";
        return;

    }
    remove_edge(src,dst,graph_);

}

float L1SkeletonGraph::shortestpath(vertex_descriptor src, vertex_descriptor target){
    std::vector<int> vertex_index_map(num_vertices(graph_));
    for (size_t i=0; i<vertex_index_map.size(); ++i) {
        vertex_index_map[i] = i;
    }

    std::vector<vertex_descriptor> parents(num_vertices( graph_),src);
   std::vector<float> distances(num_vertices(graph_), 0.0f);
    auto dmap = make_iterator_property_map(distances.begin(),
                                           get(vertex_index, graph_));
    auto pmap = make_iterator_property_map(parents.begin(),
                                           get(vertex_index, graph_));
    auto wmap = get(&edge_info::length,  graph_);
    dijkstra_shortest_paths(graph_, src,
                             weight_map(wmap).predecessor_map(pmap).distance_map(dmap));
    return distances[target];

}
float L1SkeletonGraph::shortestpath(vertex_descriptor src, vertex_descriptor target,
									std::vector<vertex_descriptor> &path){
    std::vector<int> vertex_index_map(num_vertices(graph_));
    for (size_t i=0; i<vertex_index_map.size(); ++i) {
        vertex_index_map[i] = i;
    }

    std::vector<vertex_descriptor> parents(num_vertices( graph_),src);
   std::vector<float> distances(num_vertices(graph_), 0.0f);
    auto dmap = make_iterator_property_map(distances.begin(),
                                           get(vertex_index, graph_));
    auto pmap = make_iterator_property_map(parents.begin(),
                                           get(vertex_index, graph_));
    auto wmap = get(&edge_info::length,  graph_);
    dijkstra_shortest_paths(graph_, src,
                             weight_map(wmap).predecessor_map(pmap).distance_map(dmap));
    path.clear();
    if (distances[target] < 1e6){
    	vertex_descriptor v = target;
    	while (v!= src){
    		path.push_back(v);
    		v = pmap[v]; 
    	}
    	path.push_back(src);
    	std::reverse(path.begin(),path.end());
    	//int vt=graph_[target].m_index;
    	//int vs=graph_[src].m_index; 	path.insert(path.begin(),parent);
    }                         
                             
    return distances[target];

}


void L1SkeletonGraph::dijkstras(vertex_descriptor vertex_from,
								std::vector<float>  &distances,
                                std::vector<vertex_descriptor> &parents ){
    std::vector<int> vertex_index_map(num_vertices(graph_));
    for (size_t i=0; i<vertex_index_map.size(); ++i) {
        vertex_index_map[i] = i;
    }

    auto dmap = make_iterator_property_map(distances.begin(),
                                           get(vertex_index, graph_));
    auto pmap = make_iterator_property_map(parents.begin(),
                                              get(vertex_index, graph_));
    auto wmap = get(&edge_info::length,  graph_);
    dijkstra_shortest_paths(graph_, vertex_from,
                 weight_map(wmap).predecessor_map(pmap).distance_map(dmap));

}
bool L1SkeletonGraph::bfs(vertex_descriptor src, vertex_descriptor target,
						  int max_depth){
   graph_traits<SGraph>::vertices_size_type d[num_vertices(graph_)];
   std::fill_n(d, num_vertices(graph_), 0);
   bool found = false;
 // custom_bfs_visitor<SGraph::vertices_size_type*> vis(5,d);
  try{
    breadth_first_search(graph_, src,
		visitor(make_bfs_visitor(searchVisitor( d, graph_[target].id,max_depth,
												on_tree_edge()))));
  }catch(BfsException &e){
	  found = e.success();
  }
  return found;
}
bool L1SkeletonGraph::bfs(int src, int target, int max_depth){
	vertex_descriptor s = vertex(src,graph_), //src/target are original vertex id
		t = vertex(target,graph_);
	return bfs(s,t,max_depth);
}
/*float L1SkeletonGraph::graphDistance(vertex_descriptor src, vertex_descriptor dst){
	return MAX_F;
}*/
std::pair<bool,std::vector<L1SkeletonGraph::vertex_descriptor> > 
	L1SkeletonGraph::nearAdjacent(vertex_descriptor v1, vertex_descriptor v2){

  	std::vector<vertex_descriptor> path;
 	float l = shortestpath(v1,v2,path);
 	bool near_adj=  l > 2 * cv::norm(graph_[v1].p.pos-graph_[v2].p.pos) ;
	return std::make_pair(near_adj,path);
}


void L1SkeletonGraph::collapseEdge(edge_descriptor edgeToCollapse){
	vertex_descriptor v1=source(edgeToCollapse,graph_),
		v2 = target(edgeToCollapse,graph_);
	qDebug() <<"Collapsing edge  " << graph_[v1].id <<"--" << graph_[v2].id
			 << " of length " << graph_[edgeToCollapse].length;
	cv::Point2f newpos = 0.5f *(graph_[v1].p.pos + graph_[v2].p.pos);
	graph_[v1].p.pos = newpos;	//also update internal pts of v1's incident edges	
	adjacency_iterator ai, ai_end;
	
	for (tie(ai,ai_end) = adjacent_vertices(v2,graph_);ai != ai_end; ++ai){
		edge_descriptor incidentEdge = edge(v2,*ai, graph_).first;
		if (*ai != v1){
			edge_descriptor e; bool exists;
			tie(e,exists) = edge(v1, *ai, graph_);
			if (!exists){
				e = add_edge(v1, *ai, graph_).first;
				graph_[e].id = graph_[incidentEdge].id;
				graph_[e].visited = false;
				//updateEdgeGeometry()
				graph_[e].internalPts = graph_[incidentEdge].internalPts;
			}
		}		
	}
	
	clear_vertex(v2,graph_);
	remove_vertex(v2,graph_);

	//update internal pts of all edges incident to v1
	for (tie(ai,ai_end) =adjacent_vertices(v1,graph_);ai!=ai_end;++ai){
		edge_descriptor e = edge(v1, *ai, graph_).first;
		std::vector<cv::Point2f> pts = graph_[e].internalPts;
		if (cv::norm(pts[0]-newpos) < cv::norm(pts[pts.size()-1]-newpos)){
			pts[0] = newpos;
		}else{
			pts[pts.size()-1] = newpos;
		}
		//apply simple smoothing
		updateEdgeGeometry(e,pts);		
		updateEdgeGeometry(e,graph_[e].internalPts);
		//	updateEdgeGeometry(e,graph_[e].internalPts);
	}

}

void L1SkeletonGraph::collapseEdgesShorterThan(float minLen){
	qDebug() <<"[SGraph] remove edges shorter than " <<minLen ;
	un_edge_iterator it,it_end;
	edge_descriptor shortest;float shortestlen = Common::inf;
	for (tie(it,it_end) = edges(graph_);it!=it_end;++it){
		if (graph_[*it].length < 1e-6){
			qDebug() <<"[SGraph] error in edge length " << graph_[*it].length<< endl;
		}else if (graph_[*it].length  < shortestlen){
			shortest = *it;
			shortestlen = graph_[*it].length;
		}
	}
	if (shortestlen < minLen){
		//		qDebug() <<"[SGraph] remove shortest edge "<< source(shortest,graph_)
		//		 <<"-"<< target(shortest,graph_) << " of len " << shortestlen<< endl;
		collapseEdge(shortest);
	
		collapseEdgesShorterThan(minLen);
	}
	//reset vIndex
	vertex_iterator vit,vit_end;
	std::map<int,int> newVIndex;int i=0;
	for (tie(vit,vit_end)= vertices(graph_); vit!=vit_end;++vit){
		newVIndex[graph_[*vit].id] = i++;
	}
	vIndex_ = newVIndex;
    assignVertexType();
	updateEdgeIndex();
	computeEdgeWeight();
}
