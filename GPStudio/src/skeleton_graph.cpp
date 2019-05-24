#include "skeleton_graph.h"
#include <QDebug>
#include <string>
#include <sstream>
//-----------definitions for CompactSkeletonGraph ------------
/*
CompactSkeletonGraph::~CompactSkeletonGraph(){
	qDebug() << "calls compact skeleton graph destructor" ;
	//edge_map_.clear();

	qDebug() << "end calling compact skeleton graph destructor";
	///		edge_map_.erase(edge_map_.begin(), edge_map_.end());
	std::map<EdgeMapKey,std::vector<int> >::iterator iter;

	for (iter= edge_map_.begin();iter!=edge_map_.end();++iter){
		qDebug() << "clearing edge (" << iter->first.first << "," 
				 << iter->first.second << ") point set of size %i "
				 <<   iter->second.size();
		(iter->second).clear();
		edge_map_.erase(iter);

	}
}
*/
int CompactSkeletonGraph::getVertexId(int pointId){
	std::vector<int>::iterator it = find(vertex_map_.begin(),
										 vertex_map_.end(),pointId);
	if (it != vertex_map_.end()){
		return it-vertex_map_.begin();
	}else{
		return -1;
	}
}

void CompactSkeletonGraph::printEdgeMap(){
	std::map<EdgeMapKey,std::vector<int> >::iterator iter;
	for (iter= edge_map_.begin();iter!=edge_map_.end();++iter){
		std::ostringstream oss;
		oss <<"("<<(iter->first).first <<"," << (iter->first).second <<") "
			<< "w: " << iter->second.size()<<" { " ;
		if (iter->second.size()==0) qDebug() <<"Error: empty vertex list";
		for (size_t i=0; i < iter->second.size() ;++i){
			oss << iter->second[i] <<" ";
		}
		oss <<"}";

		qDebug() << oss.str().c_str();
	}
	//print extra vertex info
	qDebug()<<"Total # of vertices: " <<  vertex_map_.size() <<" nvertices= " 
			<< nVertices ;
   
}
void CompactSkeletonGraph::printVertexMap(){
	std::ostringstream oss;
	for (size_t i=0; i <vertex_map_.size();++i){
		cv::Point2f p = points_[vertex_map_[i]];
		oss << " " << i <<"/"<<vertex_map_[i] <<"(" << p.x <<","
			<< p.y <<")";
		if (((int)i)%4==0) 
			oss << "\n";

	}
	qDebug() << oss.str().c_str();
}

void CompactSkeletonGraph::printVertex(int i){
	std::ostringstream oss;
	oss<<" vertex " << i<< " [out]";
	DGraphEdge *e = vertices[i].outHead;
	while(e){
		oss << " "<< e->target;
		e = e->nextOut;
	}
	oss << " [in] ";
	e = vertices[i].inHead;
	while(e){
		oss << " " << e->source;
		e = e->nextIn;
	}
	qDebug() << oss.str().c_str();
}

//-----------definitions for SkeletonGraph -------------------
SkeletonGraph::SkeletonGraph(BranchList branches, 
							 BridgeMap bridges,
							 std::vector<SkeletonPoint> bridgePoints,
							 int n):DGraph(n){
	initialize(branches, bridges, bridgePoints);
	//print();
}

int SkeletonGraph::addVertex(int i,SkeletonPoint p){

	qDebug() <<"add vertex #" << i << " with skeleton id " << p.id;
	if (reverseVertexIndex_.find(p.id)!=reverseVertexIndex_.end()){
		//vertex has already been added
		return i;
	}
	vertices_[i] = p;
	reverseVertexIndex_[p.id] = i;
	return i+1; 
}

void SkeletonGraph::initialize(BranchList branches, BridgeMap bridges,
							   std::vector<SkeletonPoint> bridgePoints){
	//initialize vertex map
	vertices_.resize(nVertices);

	qDebug() << "initialize graph from " << branches.size() << " branches "
			 << bridges.size() << "bridges and " << bridgePoints.size() 
			 << " bridge points. ";
	int ct = addEdgesFromBranches(branches);
	ct = addEdgesFromBridges(branches, bridges, bridgePoints,ct);

	// identify node type based on degrees by iterating 
	// through all vertices in the graph
	assignVertexType();
	//printVertexMap();

}
int SkeletonGraph::addEdgesFromBridges(BranchList branches, 
										BridgeMap bridges,
						  std::vector<SkeletonPoint> bridgePoints,
										int nExistingVertices){
	// join bridge points with branch endpoints*/
	int ct =nExistingVertices;
	auto iter= bridges.begin();
	int bridge_i=0;

	while(iter!=bridges.end()){

		int bridgePtId = iter->first;
		auto bridgePt = find_if(bridgePoints.begin(),bridgePoints.end(),
								[&bridgePtId](const SkeletonPoint &v)
								{return v.id==bridgePtId;});
		if (bridgePt == bridgePoints.end()) {
			iter++;
			continue;
		}
		// set bridgePt in vertices_ if it is not yet defined
		if (vIndex(bridgePt->id) == -1){
			ct = addVertex(ct,*bridgePt);
		}
		bridge_i++;	
	
		float branchCosThresh = -0.9;//Parameters::getInstance()->skeleton_branch_cos_thresh;
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
				if  (cv::norm(e1.pos-bridgePt->pos) 
					 < cv::norm(e2.pos-bridgePt->pos)){
					// bridge is connect to head
					target = e1;
					vec1 = branch_i[1].pos-e1.pos;
					qDebug() <<"bridge is head " ;
				}else{
					// bridge is connected to tail
					target = e2;
					int branchLen = branch_i.size();
					vec1 = e1.pos-branch_i[branchLen-2].pos;
					qDebug() <<"bridge is tail " ;

				}
				vec2 = bridgePt->pos - target.pos;
					
				qDebug() <<"angle between " << vec1.x <<"," << vec1.y
						 << " and " << vec2.x <<"," << vec2.y;
				float angle = L1Skeleton::angleBetween(vec1,vec2);
				float cosTest = cos( L1Skeleton::angleBetween(vec1,vec2));
				if (cosTest<= branchCosThresh){
					// add edge between bridgePt and target
                    int targetVid = vIndex(target.id), bridgeVid = vIndex(bridgePt->id);
                    try{
                    addNewEdge(targetVid, bridgeVid,1);
                    addNewEdge(bridgeVid, targetVid,1);
                    }catch (int e){
                        qDebug() << "Error:can not add edge between "<<targetVid<<" and " << bridgeVid;
                    }
                    qDebug() << "angle: " << cosTest <<" (" << angle
                             << " stitch bridge point "
                             << bridgePt->id << " with branch "
                             << *vit <<" endpoint " << target.id;
				}else{
					qDebug() << "angle: " << cosTest << " (" << angle
							 << " can not stitch bridge point " 
							 << bridgePt->id << " with branch "
                             << *vit <<" endpoint " << target.id;
				}
			}
		}
		iter++;
	}
	qDebug() << " finished processing bridge points " ;
	return ct;
}
int SkeletonGraph::addEdgesFromBranches(const BranchList &branches){
	//	addEdgesFromBridges(branches, bridges, bridgePoints);
	//assigmVertexType();

	int ct = 0;
	for (size_t i=0;i< branches.size();++i){
		qDebug() << "process branch " << i ;
		//iterate over all branch, add edges along each
		const std::vector<Candidate> &branch = branches[i];

		// insert first point on branch to vertex map
		ct=addVertex(ct, branch[0]);

		for (size_t j=0; j<branch.size()-1;++j){
			ct = addVertex(ct,branch[j+1]);
			int s = vIndex(branch[j]), t = vIndex(branch[j+1]);
			if (!edgeExists(s,t)){
				addNewEdge(s,t,1);
				addNewEdge(t,s,1);
			}
			
		}
	}
	return ct;
}
bool SkeletonGraph::isTurningPoint(int i){
	DGraphVertex v = vertices[i];
	//	int n1 = v.outHead->target;
	//int n2 = v.outTail->target;
	//int id1 = reverseVertexIndex_.find (n1);
	cv::Point2f p = vertices_[i].pos;
	if (v.outHead->target <0 || v.outHead->target >=vertices_.size()){

		qDebug() << "error: outHead vertex id out of bound " << i ;
		return false;
	}
	if (v.outTail->target <0 || v.outTail->target >=vertices_.size()){

		qDebug() << "error: outTail vertex id out of bound " << i ;
		return false;
	}


	cv::Point2f n1 = vertices_[ v.outHead->target].pos;
	cv::Point2f n2 = vertices_[ v.outTail->target].pos;
	cv::Point2f v1 = n1-p, v2 = n2-p; //
	if (v.outHead->target == v.outTail->target){
		qDebug() <<"error: v.outHead and tail are the same. outsize: " << 
			v.outSize;
	}else{
		qDebug() <<"outsize = " << v.outSize<< " head:" << v.outHead->target
				 <<" tail:" << v.outTail->target << " pi " << i;
		qDebug() << "angle between v1,v2 " << L1Skeleton::angleBetween(v1,v2)
				 << " cos " <<cos(L1Skeleton::angleBetween(v1,v2)) ;
	}
	if (cos(L1Skeleton::angleBetween(v1,v2)) < -0.7){
		return false;
	}
	return true;
	
}

void SkeletonGraph::assignVertexType(){
	
	vertex_types_.resize(nVertices);

	const DGraphEdge *edge;
	int unknowns=0,terminals=0,internals=0, intersections=0;
	for (size_t i=0; i< nVertices ; ++i){

		switch (vertices[i].outSize){
		case 0:			
			vertex_types_[i] = SkeletonVertexType::UNKNOWN;
			unknowns++;
			break;
		case 1:
			vertex_types_[i] = SkeletonVertexType::TERMINAL;
			terminals++;
			break;
		case 2:
			if (isTurningPoint(i)){
				vertex_types_[i] =SkeletonVertexType::INTERSECTION;
				intersections++;
			}else{ 
				vertex_types_[i] =SkeletonVertexType::INTERNAL;
				internals++;
			}
			break;
		default:
			vertex_types_[i] = SkeletonVertexType::INTERSECTION;
			intersections++;
			break;
		}
	}
	// if node is internal, check the angle between 2 edges, then mark as internal if
	// angle close to 90 degrees .

	qDebug() << "unknowns: " << unknowns <<" terminals: " << terminals
			 << " internal: " 
			 << internals << " intersection: " << intersections;
	nCompactVertices = terminals+intersections;
}

void SkeletonGraph::printVertexMap(){
	std::string nodetypes[] = {"unknown","terminal","internal","intersection"};
	
	for (size_t i=0; i < vertices_.size();++i){

		qDebug() << i << " id=" << vertices_[i].id << " pos=("
				 <<vertices_[i].pos.x << "," << vertices_[i].pos.y 
				 << ") " << nodetypes[(int)vertex_types_[i]].c_str();
	}
}


void SkeletonGraph::traceSimplePath( DGraphEdge *edge, 
									 std::vector<int> &edgePoints){

	int s = edge->source, t = edge->target; // id in the current gaph
	edgePoints.push_back(s); //add to edge points in compact graph

	if (vertex_types_[t] == SkeletonVertexType::INTERSECTION || 
		vertex_types_[t] == SkeletonVertexType::TERMINAL){
		
		edgePoints.push_back(t); // edge points are indexed based on the current graph

	}else if (vertex_types_[t] == SkeletonVertexType::INTERNAL){
		DGraphEdge *e = vertices[t].outHead;
		if (e->target == s)
			e = e->nextOut;
		traceSimplePath(e,edgePoints);
	}
		
	
}
bool SkeletonGraph::containsDuplicates(std::vector<int> endPoints){
	
	for (int i=0;i<endPoints.size()-1;++i){
		if (endPoints[i] == endPoints[i+1])
			return true;
			
	}
	return false;
}

void SkeletonGraph::buildCompactGraph(CompactSkeletonGraph *compact){
	// create compact skeleton graph
	//compact = new CompactSkeletonGraph(nCompactVertices);

	// loop over all non-internal vertices, index them in vertex map
	for (int i=0; i < nVertices; ++i){
		compact->points_.push_back(vertices_[i].pos);
		if (vertex_types_[i] == SkeletonVertexType::INTERSECTION 
			|| vertex_types_[i] == SkeletonVertexType::TERMINAL){
			compact->vertex_map_.push_back(i);
			
		} 
	}
	// compute edges in the compact graph by tracing from non-internal 
	// vertices.
	for (int i=0; i< nCompactVertices;++i){
		int index =compact->vertex_map_[i];
		DGraphEdge *edge=vertices[index].outHead;
		
		while (edge){

			std::vector<int> edgePoints;
			traceSimplePath(edge,edgePoints);
			//check duplicates in edgePoints
			/*
			if (containsDuplicates(edgePoints)){
				qDebug() <<"Error: duplicate point in path";
				exit(1);
				}*/
			

			int dst = compact->getVertexId(edgePoints.back());
			
			int src = compact->getVertexId(edgePoints.front());
			if (src >=0 && dst >=0 && edgePoints.size()>0){
				/**	compact->addNewEdge(src, dst,edgePoints.size());**/
				(compact->edge_map_)[std::make_pair(src,dst)] = edgePoints;
			}
			edge = edge->nextOut;
		}
		
	}
	// iterate edge map, remove noisy edges>>>>>>>>>>>>>>>>>

	std::vector<std::pair<int,int> > edgeToRemove;
	findDanglingEdges(compact, edgeToRemove);
	
	for (auto iter=edgeToRemove.begin(); iter!= edgeToRemove.end();
		 ++iter){
		int k = compact->edge_map_.erase(*iter);
		qDebug()<< "successful removed edge "<<iter->first<<","<< iter->second
				<<" ? " << k;
				}
	//	compact->printEdgeMap();
	//	compact->printVertexMap();

}

void SkeletonGraph::findDanglingEdges(CompactSkeletonGraph* compact, 
					   std::vector<std::pair<int,int> > &edgeToRemove){
	std::map<EdgeMapKey, std::vector<int> >::iterator iter;
	for (iter=compact->edge_map_.begin();iter!=compact->edge_map_.end();
		 ++iter){

		int i = iter->first.first, j= iter->first.second;
		// if edge (i,j) in compact graph has length 1 (2 points)
		if (iter->second.size() == 2 ){//!!@@@@ always add edge
			int u = iter->second[0], v = iter->second[1];
			if (vertex_types_[u] == SkeletonVertexType::TERMINAL ||
				vertex_types_[v] == SkeletonVertexType::TERMINAL){
				edgeToRemove.push_back(iter->first);
				continue;
			}
				
		}
		compact->addNewEdge(i,j,iter->second.size());

	}
	qDebug()<<"Mark " << edgeToRemove.size() <<" edges to be removed.";
}

//----------------------------------------------------------------
//Note: make sure has initialized points_ in compact_skeleton_graph
