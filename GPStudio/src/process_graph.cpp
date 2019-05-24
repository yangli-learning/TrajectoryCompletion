#include "process_graph.h"
 
#include <opencv2/flann/flann.hpp>
#include "parameters.h"
#include <QDebug>
#include <sstream>
#include <string>
#include <boost/geometry/geometries/point_xy.hpp>
#include "gps_trajectories.h"
#include "gps_trajectory.pb.h"
#include "l_shape_iterator.h"
#include <utility>
#include <fstream>

ProcessGraph::ProcessGraph(CompactSkeletonGraph *g, float h, float w):
	graph_(g),height_(h),width_(w){
	step_=0.25f*Parameters::getInstance()->skeleton_sample_cell_size;//0.5*
	process();
	clusterEndpoints();
}

/*
~ProcessGraph(){
    delete
}*/

void ProcessGraph::process(){
	bool finished = false; // finished implies all edges are marked 
	                       // as final or has exceeded bounding box
	initJunctionEdges();
	int maxIter=25;
	float step = step_;
	int i=1;
	while (!finished && i<=maxIter){
		qDebug() <<"=========== iter " << i <<" ===============";
		finished = extendJunctionEdges(step*i); 
		checkIntersections();
		i++;
	}
}
void ProcessGraph::initJunctionEdges(){

	for (auto it = graph_->edge_map_.begin();
		 it!=graph_->edge_map_.end(); ++it ){

		if (it->first.first == it->first.second)
			continue;

		DGraphVertex v1 = graph_->vertices[it->first.second];//edge target
		ExtendedEdge e={it->first.first, it->first.second, cv::Point2f(0,0),
						v1.outSize >1,0};
		extended_edges_.push_back(e);
	}
	qDebug() << "initialized " << 
		extended_edges_.size() << 
		" out of " <<graph_->edge_map_.size()<<" edges to be extended ";
}

bool ProcessGraph::extendJunctionEdges(float len){
	//initialize the junction_edges_
	bool finished = true;
	for (int i=0; i< extended_edges_.size();++i){
		ExtendedEdge &ee = extended_edges_[i];		
		if (ee.final) 
			continue;
		finished = false;
		std::pair<int,int> key;
		key.first = ee.head;
		key.second = ee.tail;

		auto edgeIt = graph_->edge_map_.find(key);
		if ( edgeIt==graph_->edge_map_.end()){
			qDebug() <<"can not find " << key.first <<","<< key.second<<
				"in edge map!";
			continue;
		}
		std::vector<int> &polyline = edgeIt->second;	
		size_t n = polyline.size();
		int tailId = polyline.back();
	   

		if (ee.dir == cv::Point2f(0,0)){
			// add one vertex to end of polyline, update points_ 
			// and vertices_ 		
			cv::Point2f dir = (n>2)? graph_->points_[tailId]-graph_->points_[polyline[n-3]]:
				graph_->points_[tailId]- graph_->points_[polyline[n-2]];

			Common::normalize(dir,ee.dir); // = 
			graph_->points_.push_back(graph_->points_[tailId]);
			polyline.push_back((int)(graph_->points_.size())-1);
			graph_->vertex_map_[ee.tail] = polyline.back();

		}
		int tailVertexId = polyline.back();

		cv::Point2f s = graph_->points_[tailVertexId];
		if (s.x < 0 || s.y <0 || s.x >= width_ || s.y >= height_){
			qDebug()<<" outside bounding box";
			ee.final = true;
		}else{
			graph_->points_[tailVertexId] += len*ee.dir;	
		}
		
	}	
	return finished;
}

void ProcessGraph::checkIntersections(){
	// fill edge linestrings and query linestring
	std::vector<value> segments, queries;
	auto it = extended_edges_.begin();
	for ( ;it != extended_edges_.end();++it){

		std::vector<int> pointIds = graph_->edge_map_[std::make_pair(it->head, 
																it->tail)];
		int i = it-extended_edges_.begin();
		segment segQuery;

		for (size_t j=0; j < pointIds.size()-1; ++j){
			cv::Point2f &v1 = graph_->points_[pointIds[j]];
			cv::Point2f &v2 = graph_->points_[pointIds[j+1]];
			point p1(v1.x,v1.y), p2(v2.x,v2.y);
			segment seg( p1, p2 );
			segments.push_back(std::make_pair(seg,i));
			
			if (j== pointIds.size()-2)
				segQuery = seg; 
		}
		if (! it->final){
			queries.push_back(std::make_pair(segQuery,i));
		}
	}
	qDebug() << segments.size() << " polylines " << queries.size()
			 <<" segments";

    float r = step_*2;
	rtree rt(segments.begin(), segments.end());

	for (auto qit = queries.begin(); qit !=queries.end();++qit){
		// find all nearby linestrings to query segment
		auto eQuery =std::make_pair(extended_edges_[qit->second].head, 
									extended_edges_[qit->second].tail);
		rtree::const_query_iterator nit;
		for (nit=rt.qbegin(bgi::nearest(qit->first,10)); nit!= rt.qend();
				 ++nit){
			bool stopSearch =false;
			auto eNeighbor=
				std::make_pair(extended_edges_[nit->second].head, 
							   extended_edges_[nit->second].tail);

			if (different_id(*nit,qit->second) &&
				                          eQuery.first!=eNeighbor.second){
				double d = bg::distance(qit->first,nit->first);
				qDebug() <<"distance between a segment on " <<  qit->second
						 << " and endpoint " << nit->second << " distance is " << d;
				if (d>r){
					extended_edges_[qit->second].final = false;
					break;
				}
				ExtendedEdge &ee = extended_edges_[qit->second];
				extended_edges_[nit->second].nNeighbors++;
				std::vector<int> edgePoints = graph_->edge_map_[std::make_pair(ee.head,ee.tail)];
				std::vector<point> output; 
				bg::intersection(qit->first, nit->first, output);

				if (output.size()>0){
					// set segment endpoint to intersection,update points_
					point p =output[0];
                    if ( bg::distance(qit->first.first,p) ==0){
						qDebug() <<" intersection is on other endpoint d=" 
							   << bg::distance(qit->first.first,p);
						extended_edges_[qit->second].final = false;
						continue;
					}

					graph_->points_[edgePoints.back()] = boostToCVPoint(p);
					extended_edges_[qit->second].final = true;
                    //stopSearch = true; //! @todo check the consequence
				}else {
					//project to segment
					cv::Point2f p (qit->first.second.get<0>(),qit->first.second.get<1>()),
						e1(nit->first.first.get<0>(), nit->first.first.get<1>()),
						e2(nit->first.second.get<0>(),nit->first.second.get<1>());

					projectPointOnSegment(p,e1,e2,graph_->points_[edgePoints.back()]);
					extended_edges_[qit->second].final = true;					
                    stopSearch = true;
				}
				//check if the projection is on a query segment
				cv::Point2f p1 = boostToCVPoint(nit->first.first);
				cv::Point2f p2 = boostToCVPoint(nit->first.second);
				std::vector<int> npts=graph_->edge_map_[eNeighbor];
					
				if (p2 == graph_->points_[npts.back()]){

					qDebug()<<"Note: segment is an extension";

					// test the angle of queries
					float ang=L1Skeleton::angleBetween(extended_edges_[nit->second].dir, 
													   extended_edges_[qit->second].dir);
				   	if (abs(cos(ang))>0.9){
						qDebug()<<"similar: angle, connect";
						graph_->points_[npts.back()] =graph_->points_[edgePoints.back()];
						//boostToCVPoint(p);
						extended_edges_[nit->second].final =true;
					}
				}
				/*-----------------------------------
				// if distance from intersection to either endpoint is further than r,
				// split edge r at the intersection
				cv::Point2f v3 = graph_->points_[edgePoints.back()];
				std::vector<int> neighborEdgePts = graph_->edge_map_[eNeighbor];
				cv::Point2f v1 = graph_->points_[neighborEdgePts.front()];
				cv::Point2f v2 = graph_->points_[neighborEdgePts.back()];
				if (cv::norm(v3,v1)>R && cv::norm(v3,v2)>R){
					splitExtendedEdge(eNeighbor.first,eNeighbor.second,eQuery.second,
									  );					
				}
				//-------------------------------------*/
				if (stopSearch) break;
			}
		}
	}
}
/*
void ProcessGraph::splitExtenededEdge(int v1, int v2, int v3, cv::Point2f insertP ){
	//modify graph_ connectivity
	//edit edge_map
	graph_->removeDirectedEdge(v1,v2);
	graph_->addEdge(v1,v3,1);
	graph_->addEdge(v3,v2,1);
	graph_->edge_map_.erase(std::make_pair(v1,v2));
	std::vector<int> seg1,seg2;
	graph_->edge_map_[std::make_pair(v1,v3)] = seg1;
	graph_->edge_map_[std::make_pair(v3,v2)] = seg2;
}
*/
std::vector<std::array<int,2> >  ProcessGraph::getSortedClusterKeys(){
	std::vector<std::array<int,2> > keys;
	for (auto it=trajectory_clusters_.begin();it!=trajectory_clusters_.end();++it){
		keys.push_back(std::array<int,2>(it->first));
	}
	
	std::sort(keys.begin(),keys.end(),[this](std::array<int,2> a,  std::array<int,2> b){
			return trajectory_clusters_[a].size() > trajectory_clusters_[b].size();
		});
	return keys;
}
bool ProcessGraph::projectPointOnSegment(cv::Point2f in, 
										 cv::Point2f p1,
										 cv::Point2f p2, 
										 cv::Point2f &out){
	if (p1==p2){
		out =  p1;//mergeTopPoints(in,p1);
		return false;
	}
	cv::Point2f e = p2-p1;
	cv::Point2f v = in-p1;
	cv::Point2f en, vn;
	Common::normalize(e,en);
	Common::normalize(v,vn);
	float projLen = v.x*en.x + v.y*en.y;
	float eLen = cv::norm(e);
	if (projLen <=0 || projLen>=eLen){
		out =  (projLen <=0)? p1:p2;
		return false;
	}else{
		out = p1+en*projLen;
		return true;
	}
}

void ProcessGraph::clusterEndpoints(){
    //sort extended_edges_ in descending order of # neighbors
	sort(extended_edges_.begin(),extended_edges_.end(), 
		 [](const ExtendedEdge &e1,const ExtendedEdge &e2){
			 return e1.nNeighbors>e2.nNeighbors;
		 });
    float r = Parameters::getInstance()->skeleton_cluster_radius
                *Parameters::getInstance()->skeleton_sample_cell_size;

	bool visited[extended_edges_.size()] ;//= {false};
	std::vector<pvalue> pointList;
	for (size_t i=0; i < extended_edges_.size() ;++i ){
		std::vector<int> pids =   graph_->getEdgePointIds(extended_edges_[i].head, 
												  extended_edges_[i].tail);

		cv::Point2f ep = graph_->points_[pids.back()];
		pointList.push_back(std::make_pair( point(ep.x, ep.y),i));
	}
	prtree rt(pointList.begin(),pointList.end());
	std::ostringstream oss;
	for (size_t i=0; i < pointList.size() ;++i ){
		if (visited[i]) continue;
		prtree::const_query_iterator nit;
		std::unordered_set<size_t> pset={i};
		for (nit=rt.qbegin(bgi::nearest( pointList[i].first,10));nit!=rt.qend();++nit){
			double d = bg::distance(pointList[i].first,nit->first);
			if (d>r){
				break;
			}
			pset.insert( nit->second);
			int twinId = getTwinExtendedEdge(nit->second)-extended_edges_.begin();
			pset.insert(twinId);
			oss << " " << nit->second<<"|" <<twinId<<"("<<d<<") " ;
			visited[nit->second] = true;

		}
		clusters_.push_back(std::vector<size_t>(pset.begin(),pset.end()));
		qDebug()<< i <<": "<<oss.str().c_str();
		oss.str("");
		cv::Mat m = cv::Mat(pset.size(),pset.size(), CV_32S, cv::Scalar(1));
		cluster_adj_.push_back(m);
		oss << m;
		qDebug() << oss.str().c_str();
		oss.str("");
	}
}

rtree ProcessGraph::makeRTreeFromEdges(){
	std::vector<value> segments;
	auto it = extended_edges_.begin();
	for ( ;it != extended_edges_.end();++it){

		std::vector<int> pointIds = graph_->edge_map_[std::make_pair(it->head, 
																it->tail)];
		int i = it-extended_edges_.begin();

		for (size_t j=0; j < pointIds.size()-1; ++j){
			cv::Point2f &v1 = graph_->points_[pointIds[j]];
			cv::Point2f &v2 = graph_->points_[pointIds[j+1]];
			point p1(v1.x,v1.y), p2(v2.x,v2.y);
			segment seg( p1, p2 );
			segments.push_back(std::make_pair(seg,i));
		}
	}
	rtree rt(segments.begin(), segments.end());
	return rt;
}

float ProcessGraph:: computeDistRatio(point p, 
									  std::vector<value> result_n){
	IndexVal<float> nval;
	nval.index=-1;nval.val = Common::inf;
	std::vector<IndexVal<float> > branchDist(extended_edges_.size(),nval);
	for (size_t i=0; i<result_n.size();++i){
		int id = result_n[i].second;
		float d = bg::distance(p,result_n[i].first);
		if (branchDist[id].val > d){
			
			branchDist[id].index= id;
			branchDist[id].val=d;
		}
	}
	partial_sort(branchDist.begin(), branchDist.begin()+3,
				 branchDist.end());
	/*
		qDebug()<<" d0=" << branchDist[0].index<<"|"<< branchDist[0].val
			<<" d1="<<  branchDist[1].index<<"|"<< branchDist[1].val
			<<" d2=" << branchDist[2].index<<"|"<< branchDist[2].val;
	*/
	if (branchDist[2].index==-1){
		return 0;//no ambiguity
    }else if (twinEdge(branchDist[0].index, branchDist[1].index)){
		return branchDist[0].val/branchDist[2].val;
	}else{
		return branchDist[0].val/branchDist[1].val;
	}
}

void ProcessGraph::validateTopologyByGPS(
					osg::observer_ptr<GPSTrajectories>  trajectories,
					osg::Vec3 min_corner,  float resolution){

	GPSTrajectories::GpsTrajs gps_trajs = trajectories->getGpsTrajs();
	qDebug() <<"validate through " << gps_trajs.size() << "trajectories";

	//step 1 project each point on the extended edge using rtree.
	//each point is labeled J (J>=0), if not certain, J=-1

	rtree rt = makeRTreeFromEdges();
    int shortTraj=0, ambTraj=0 ;
    size_t maxTraj = 5000;
	int ct_all=0, ct_discarded=0;
	for  (size_t i = 0, i_end = gps_trajs.size(); i < i_end; ++i) {
		int seg_prev = -1;
		std::shared_ptr<GpsTraj>& gps_traj = gps_trajs[i];
		if (gps_traj->point_size() < 2){
            shortTraj ++;
			continue;
		}

		std::vector<int> project_idx(gps_traj->point_size(),-1);

		std::vector<std::array<int,3> > breakpt; // [index of breakpt, prevId,and nextId]
		for (size_t j = 0, j_end = gps_traj->point_size(); j < j_end; ++j) {
			ct_all++;
			point p = toBoostPoint(&gps_traj->point(j), min_corner, resolution);

			std::vector<value> result_n;
			rt.query(bgi::nearest(p, 3), std::back_inserter(result_n));

			sort(result_n.begin(),result_n.end(),[p](const value& v1, const value &v2){
					return bg::distance(p,v1.first) < bg::distance(p,v2.first);
				});
            // skip if nearest distance is greater than a threshold
			// compute projection distance ratio------------
			float ratio = computeDistRatio(p,result_n);
			float distToHead=0,distToTail=0;
			value nearest = result_n[0];

			if (j>0){
				cv::Point2f cvp (gps_traj->point(j-1).x(), 
								 gps_traj->point(j-1).y()); 

				distToHead =cv::norm(cvp-boostToCVPoint(nearest.first.first));
				distToTail = cv::norm(cvp-boostToCVPoint(nearest.first.second));
			}
            if (ratio < Parameters::getInstance()->skeleton_proj_ratio_thresh ){//&& distToHead <= distToTail){
				int ang_curr;
				try{
					ang_curr = (int32_t)(gps_traj->point(j).head()); 
				}catch(...){
					qDebug()<<"Error: bad angle!";
					ct_discarded++;
					continue; //skip point with bad angle
				}


				cv::Point2f tDir = LShapeIterator::UnitDirection(ang_curr);
				/*cv::Point2f temp = tDir;
				
				if (j < j_end-1){
					temp.x = gps_traj->point(j+1).x() - gps_traj->point(j).x();
					temp.y = gps_traj->point(j+1).y() - gps_traj->point(j).y();
				}
				Common::normalize(temp,tDir );
				*/
				cv::Point2f sDir;
				Common::normalize(boostToCVPoint(nearest.first.second)- boostToCVPoint(nearest.first.first), sDir);
				int nearestId = seg_prev;
                float angle_thresh = Parameters::getInstance()->skeleton_proj_angle_thresh;
                if (tDir.dot(sDir)> angle_thresh) { //>0.8){
					nearestId = (int)nearest.second;
					project_idx[j] =nearestId;
                }else if (tDir.dot(sDir)< -angle_thresh) {//-0.8){
					nearestId =( getTwinExtendedEdge((int)nearest.second)-extended_edges_.begin());
					project_idx[j] = nearestId;
				}else{
                    ambTraj++;
					ct_discarded++;
					continue;
				}		 
				
				//vote on adjacency matrix according to previous edge
				// if seg id has changed from the previous idprev
				// cast vote (idprev,id), update idprev
				if (seg_prev>=0 && seg_prev != nearestId){
					bool success =  vote(seg_prev,nearestId );

					// add subtrajectory info to trajectory_clusters
					if (success){
						qDebug() << "add to trajectory cluster ";
                        std::array<int,3> breakInfo = {(int)j, seg_prev, nearestId};
                        breakpt.push_back(breakInfo);
                    }else{
                        std::array<int,3> breakInfo = {(int)j,-1,-1};
                        breakpt.push_back(breakInfo);
                    }
				}
				seg_prev = nearestId;		  
            }else{
				ct_discarded++;
            }
        }

		//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        if (breakpt.empty()){
            continue;
        }
		
		clusterSubTrajectories(breakpt, project_idx, i);
		if (i>maxTraj) break;
	}

	std::ostringstream oss;
	for (size_t i=0; i < cluster_adj_.size() ;++i ){
		oss << cluster_adj_[i];
		qDebug() << oss.str().c_str();
		oss.str("");
	}
	qDebug() <<"read "<<ct_all<<" pts, "<< ct_discarded 
			 <<" ambiguous points discarded\n\n";
    qDebug() << shortTraj << " trajectories with <2 vertices, "<< ambTraj <<" trajectories"
                << "with ambiguous angles";

	processTrajClusters();
	for (auto it = trajectory_clusters_.begin(); it!=trajectory_clusters_.end();
		 ++it){
		qDebug() <<it->first[0]<<","<<it->first[1]<<": " <<it->second.size();
	}

}
void ProcessGraph::processTrajClusters(){
	// clear subtrajectories in incorrect clusters
	for (auto it= trajectory_clusters_.begin(); it!=trajectory_clusters_.end();
		 ++it){
		int src = it->first[0];
		int dst = it->first[1];
		int srcTwin =  getTwinExtendedEdge(src)-extended_edges_.begin();
		int dstTwin =  getTwinExtendedEdge(dst)-extended_edges_.begin();
		size_t clusterSize = it->second.size();
		std::array<int,2> twinKey1={srcTwin , dst};
		std::array<int,2> twinKey2={src,dstTwin};
		//	bool erased = false;
		if (trajectory_clusters_.find(twinKey1) != 
			trajectory_clusters_.end()){

			auto &c = trajectory_clusters_[twinKey1];
			if (c.size() <= clusterSize){
				c.clear();
				qDebug() <<"remove: " << twinKey1[0] <<"->"<< twinKey1[1];
				//	trajectory_clusters_.erase(cit);
				//erased = true;
			}
			
		}
		if (trajectory_clusters_.find(twinKey2) != 
				  trajectory_clusters_.end()){
			auto &c = trajectory_clusters_[twinKey2];
			if (c.size() <= clusterSize){
				c.clear();
				qDebug() <<"remove: " << twinKey2[0] <<"->"<< twinKey2[1];
			}
		}
			
	}
	// remove empty clusters
	qDebug() <<"num of clusters (before): " <<  trajectory_clusters_.size();
	for(auto it=trajectory_clusters_.begin(); it!=trajectory_clusters_.end();  ){
		if (it->second.empty()){
			trajectory_clusters_.erase(it++);
		}else{
			++it;
		}

	}
	qDebug() <<"num of clusters (after): " <<  trajectory_clusters_.size();
}
point ProcessGraph::toBoostPoint(const TrajPoint* tp, osg::Vec3  min_corner, 
								 float resolution){
	float x = tp->x();
	float y = tp->y();

	x = (x - min_corner.x()) / resolution;
	y = height_ - (y - min_corner.y()) / resolution;

	return point(x,y);
}
bool ProcessGraph::vote(int i, int j){
	 //	 std::pair<int,int>   e=std::make_pair(i,j);
	 //	 int jhead = extended_edges_[j].head, jtail
	bool success = false;
	 size_t jj  =j;// getTwinExtendedEdge(j)-extended_edges_.begin();//find_if(extended_edges_.begin(), extended_edges_.end(),
	 if (jj < extended_edges_.size()){
		 
		 int ix,jjx;
		 auto cit =  find_if(clusters_.begin(), clusters_.end(), 
							 [i,jj,&ix,&jjx](std::vector<size_t> v){
				 auto eit1 = find(v.begin(),v.end(),i),
				 eit2 = find(v.begin(),v.end(),jj);
				 ix = eit1-v.begin();
				 jjx = eit2-v.begin();
				 return (eit1!=v.end() && eit2!=v.end());
			 });

		 if (cit == clusters_.end()){
			 qDebug()<<"can not find ("<< i<<","<< jj<<")";

		 }else{
			 qDebug()<<"vote:  ("<< i<<","<< jj<<") in  cluster of size " << cit->size() ;
			 int cid = cit-clusters_.begin();
			 if (cid>=0 && cid< cluster_adj_.size()){
			
				 if (ix>=0 && jjx>=0 && ix<  cluster_adj_[cid].rows && jjx<cluster_adj_[cid].cols){
					 cluster_adj_[cid].at<int>(ix,jjx)++;
					 success =true;
				 }else{
					 qDebug()<<"error: invalid i,jj. indices:" << ix <<","<<jjx<<". mat size "<<   cluster_adj_[cid].rows;
				 }
			 }else{ 
				 qDebug()<<"error: invalid cid " <<cid;   

			 }
		 }
	 }else{
		 qDebug()<<"error: cannot find twin";

	 }
	 return success;
 }
int ProcessGraph::getMaxClusterSize(){
    int max=0;
    for (auto it=trajectory_clusters_.begin();it!=trajectory_clusters_.end();++it){
        max=std::max(max,(int)it->second.size());
    }return max;

}

void ProcessGraph::clusterSubTrajectories(std::vector<std::array<int,3> > breakpt,
										  std::vector<int> project_idx, int i){
	// [debug] print out break points--------------
	std::ostringstream oss,ossBP;
	ossBP <<"breakpts ";
	for (int m=0; m<breakpt.size();++m){
		ossBP <<", " << breakpt[m][0] << "|"<<breakpt[m][1]<<"|" << breakpt[m][2];
	}
	qDebug() << ossBP.str().c_str() ;

	//---------------------------------------------
	// create forward and backward map of the ambiguous segments
	// blanks: {startPos : length}
	// blanksR: {endPos: length}
	std::map<int,int> blanks,blanksR;
	ossBP.str("");
	ossBP <<"BLANK ";
	  
	int  idx = 0;
	auto iter = std::find(project_idx.begin(), project_idx.end(),-1);
	while (iter!= project_idx.end()){

		idx = iter-project_idx.begin();
		if (blanks.find(idx)==blanks.end()){
			blanks[idx]=1;
		}
		int i = idx;
		while( (++i) < project_idx.size() && project_idx[i]==-1){
			blanks[idx]++;
		}
		ossBP << ", " << idx << ":"<< blanks[idx];
		iter = std::find(project_idx.begin()+i,project_idx.end(),-1);
		   
	}
	qDebug() << ossBP.str().c_str();

	for (std::map<int,int>::iterator it=blanks.begin();
		 it!=blanks.end();++it){
		   
		blanksR[it->first+it->second-1] = it->second;
	}
	//-------- project subtrajectory ---------------------
	// find the beginning and end of subtrajectories
	//	oss <<"traj "<< i<<": ";
	int firstProjectedIdx = breakpt.front()[0];
    int lastProjectedIdx = breakpt.back()[0];
    for (int c=0; c<project_idx.size() ; ++c){

        if (c!=0){
            oss << ", ";
        }
        oss << project_idx[c] ;

        if (project_idx[c] != -1 && c < firstProjectedIdx){
            firstProjectedIdx = c ;
        }
        if (project_idx[c] !=-1 && c>lastProjectedIdx){
            lastProjectedIdx = c;
        }
    }

	qDebug() <<"first project index " <<      firstProjectedIdx 
			 << "last project index " <<  lastProjectedIdx;
	qDebug() << oss.str().c_str();
	for (size_t k=0; k< breakpt.size(); ++k){
		int j = breakpt[k][0];
		std::array<int,2> key = { breakpt[k][1],breakpt[k][2]};
		//find first trajectory point index != -1
		int j_start = (k==0)? firstProjectedIdx: breakpt[k-1][0],
			j_end = (k==breakpt.size()-1)? lastProjectedIdx+1:breakpt[k+1][0];

		qDebug() <<i<< " <> " << j_start<<"|"<<key[0] <<"<>" <<j_end<<"|"<<key[1];
		// check blank: if either end of the subtrajectory is -1
		// move to  blankStart + n/2
		// check if left tail is -1
		if (project_idx[j_start]==-1){
			// find length of blank
			int len;
			if (blanks.find(j_start)!= blanks.end()){
				len = blanks[j_start];
                 j_start = j_start+ceil(len/2);
            }
		}
		if (project_idx[j_end-1]==-1){
			int len;
			if (blanksR.find(j_end-1)!= blanksR.end()){
				len = blanksR[j_end-1];
                 j_end = j_end - ceil(len/2);
			}
		}
        SubTrajectory subTraj = {int(i),j_start,j_end};
		qDebug() <<i<< " <> " << j_start<<"|"<<key[0] <<"<>" <<j_end<<"|"<<key[1];
		if (key[0]==-1 || key[1]==-1){
			continue;
		}
		if (trajectory_clusters_.find(key)==trajectory_clusters_.end()){
			TrajectoryCluster tc;
			trajectory_clusters_[key] = tc;
		}

		trajectory_clusters_[key].push_back(subTraj);
	}
   // std::sort(trajectory_clusters_.begin(),trajectory_clusters_.end(), )

}

std::vector<ExtendedEdge>::iterator ProcessGraph::getTwinExtendedEdge(int i){
	ExtendedEdge e = extended_edges_[i];
	return find_if(extended_edges_.begin(), extended_edges_.end(),[e](const ExtendedEdge &c){
			return (e.head == c.tail && c.tail == e.head);
		});
}
void ProcessGraph::exportGraph(const char* fname){
	std::ofstream out(fname, std::ofstream::out);
	out << toDotString() ;
	out.close();
	qDebug()<<"wrote dot graph to "<< fname;
}
QStringList ProcessGraph::clusterToQStringList(){
    QStringList slist;
        int i=0;
    for (auto it=trajectory_clusters_.begin();it!=trajectory_clusters_.end();++it){
        slist << QString("Cluster %1: %2--%3 (%4)").arg(i).arg(it->first[0]).arg(it->first[1]).arg(it->second.size());
        i++;
    }
    return slist;
}

std::vector<int> ProcessGraph::getClusterTrajectoryIds(int clusterId){
	int ct = 0;
	std::vector<int> ids;
	for (auto it=trajectory_clusters_.begin();it!=trajectory_clusters_.end();++it){
		//        slist << QString("Cluster %1: %2--%3").arg(i).arg(it->first[0]).arg(it->first[1]);

		if (ct==clusterId){
			qDebug()<<"=== Selected cluster " << clusterId<< " ====";
			for  (int i=0; i<it->second.size();++i){
				SubTrajectory st = it->second[i];
				ids.push_back(st.trajId);
			}
			break;
		}
		ct++;
	}
	return ids;
}
/** alternative method is to highlight the entire trajectory **/
void ProcessGraph::printClusterInfo(int clusterId){
	int i = 0;
	for (auto it=trajectory_clusters_.begin();it!=trajectory_clusters_.end();++it){
		//        slist << QString("Cluster %1: %2--%3").arg(i).arg(it->first[0]).arg(it->first[1]);
		if (i==clusterId){
			qDebug()<<"=== Selected cluster " << clusterId<< " of size " 
					<<  it->second.size() <<" ====";
			std::ostringstream oss;
			
			qDebug() << it->second.size();
			for (int i=0;i< it->second.size();++i){
				SubTrajectory st  = it->second[i];
				//int trajId = st.trajId;

				if (i==0){
					oss << st.trajId;
				}else{
					oss<<", " << st.trajId;
				}
				
			}
			qDebug() << oss.str().c_str() ;
			break;
		}
        i++;
    }
}

std::string ProcessGraph::toDotString(){
	float graphScale = 1.f/( 2.f*Parameters::getInstance()->skeleton_sample_cell_size);
	std::ostringstream oss;
	oss << "digraph G{\n" ;
	oss << "graph [splines=curved]\n";
	oss << "node [shape=record]\n";

	// map a label to its record and port id
	// e.g. portMap[label] = pair<recordId,portId>
	std::vector<std::pair<int,int> > portMap(extended_edges_.size());

	// boolean array to mark visited labels
	bool visited[extended_edges_.size()];// = {false};

	for (size_t i=0; i< clusters_.size();++i){
		for (size_t j=0; j<clusters_[i].size();++j){
			int eid = clusters_[i][j];

			if ( !visited[eid] ){
				visited[eid] = true;
				int twinId = getTwinExtendedEdge(eid)-extended_edges_.begin();
				visited[twinId] = true;
				portMap[eid]  = std::make_pair(eid,0);
				portMap[twinId]  = std::make_pair(eid,1);
				oss << eid << "[\n";
				oss << "\tlabel = \"<0> " <<eid <<" | <1> " <<twinId <<"\";\n";
				auto pids = graph_->getEdgePointIds( extended_edges_[eid].head,
													 extended_edges_[eid].tail);
				
				cv::Point2f midpt =0.5f*( graph_->points_[pids.front()] 
										  + graph_->points_[pids.back()]) * graphScale;
				oss <<"\tpos=\"" << midpt.x <<"," <<height_*graphScale-midpt.y<<"!\";\n";
				oss <<"]\n";
			}
		}
		qDebug()<<"write dot graph edge for matrix of size " << cluster_adj_[i].rows
				<<"x"<<cluster_adj_[i].cols;
		for (size_t r=0; r<cluster_adj_[i].rows;++r){
			for (size_t c=0;c< cluster_adj_[i].cols;++c){
				std::string colorStr = (cluster_adj_[i].at<int>(r,c)>1)?
					"\"black\"":"\"#40e0d0aa\"";
				//	if (cluster_adj_[i].at<int>(r,c)>1){
				int src = clusters_[i][r];
				int dst = clusters_[i][c];
				std::array<int,2> key = {src,dst};
				if (trajectory_clusters_.find(key)==trajectory_clusters_.end())
					continue;

				oss << portMap[src].first<<":" <<portMap[src].second<<"->"
					<< portMap[dst].first<<":" <<portMap[dst].second;
				if ( cluster_adj_[i].at<int>(r,c)>1 ){
					oss <<" [ label=\"" << cluster_adj_[i].at<int>(r,c) 
						<<"\" color="<<colorStr
						<<" fontcolor=" <<colorStr << "]\n"; 
				}else{
					oss <<"[ color="<<colorStr<<" ]\n"; 
				}
					//	}
			}
		}
	}

	oss <<"}\n";
	
	return oss.str();
}
