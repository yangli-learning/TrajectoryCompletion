
#include <ostream>
#include <QDebug>
#include <sstream>
#include <fstream>
#include "gps_trajectory.pb.h"
#include "junction_graph.h"
#include "trajectory_completion.h"

#include "junction_network.h"
#include "common.h"
#include <boost/graph/dijkstra_shortest_paths.hpp>
#include <boost/property_map/property_map.hpp>
#include <boost/graph/graph_utility.hpp>
#include "l_shape_iterator.h"

#define N_TEST_TRAJS 5000
TrajectoryCompletion::TrajectoryCompletion(SubTrajClusters *clusters,
										   GpsTrajs * trajs,
										   cv::Point2f minCoor, float height,
										   float resolution):
	clusters_(clusters),trajs_(trajs),min_corner_(minCoor),height_(height),
	resolution_(resolution)
{
	qDebug()<<"---------------------------------";
    step_ = Parameters::getInstance()->skeleton_extension_step*
              Parameters::getInstance()->skeleton_sample_cell_size;//0.5

	//initialize trajClusterMap_
	int clusterId =0;
	for (auto it = clusters_->begin();it != clusters_->end();++it){
		int s = it->first[0], t =it->first[1];
		jEdges_.push_back(std::make_pair(s,t));
		for (size_t i=0; i< it->second.size(); ++i){
			int trajId = std::get<0>(it->second[i]);
				if (trajClusterMap_.find(trajId)==trajClusterMap_.end()){
				trajClusterMap_[trajId] = std::vector<int>();
			}
			trajClusterMap_[trajId].push_back(clusterId);
		
		}
		clusterId ++;
	}
	//initialize unit_directions_
	int nDirections = 24;
	for (int i=0; i<nDirections;++i){
		cv::Point2f p = LShapeIterator::UnitDirection(i*float(360/nDirections));
		unit_directions_.push_back(p);		
	}
}


cv::Point2f TrajectoryCompletion::toCvPoint2f(const TrajPoint *p){
	return JunctionNetwork::toCvPoint2f(p,min_corner_,height_,resolution_);
}

cv::Point2f TrajectoryCompletion::fromCvPoint2f(const cv::Point2f p){
    cv::Point2f xy;
    xy.x =  p.x*resolution_ + min_corner_.x;
    xy.y = (height_- p.y)*resolution_ + min_corner_.y;
    return xy;
}

std::pair<std::vector<cv::Point2f>,bool> TrajectoryCompletion::completeTrajectory(int trajId){

	int num_success_subtraj =0;
	int trajLen = (*trajs_)[trajId]->point_size();
	int nSubTraj=0;
    std::vector<cv::Point2f> vertices;
    auto clusterIter=trajClusterMap_.find(trajId);
	std::ostringstream oss;
	float rate = 0;
    if (clusterIter!=trajClusterMap_.end() ){
		nSubTraj = clusterIter->second.size();
		//---------------------------------------------------------------
		// sort sub-trajectories by confidence score, e.g. # of unique
		// clusters sharing the segment
            CompareClusterSize comp(clusters_, jEdges_);
         sort(clusterIter->second.begin(), clusterIter->second.end(),
             comp);
		//---------------------------------------------------------------
		
        std::vector<std::vector<cv::Point2f> > denseTraj(trajLen-1);
        std::vector<bool > success(trajLen-1,false);
		
        for (size_t cId =0;cId < clusterIter->second.size();++cId){
            //reverse index into clusters to look up subtrajectory index
			subtraj_descriptor st = getSubTrajFromCluster((clusterIter->second)[cId], trajId);
            if (std::get<1>(st) < std::get<2>(st)){

				//	oss <<  "Fitting cluster "<<(clusterIter->second)[cId]<<":"
				//		<< jEdges_[cId].first<<","<<jEdges_[cId].second <<" T[" << std::get<1>(st)
				//	<<","<< std::get<2>(st) <<") \n";
				
				fitSubTrajToClusterSP(st, (clusterIter->second)[cId],denseTraj,success);
			}//else{
				//	oss <<  "Error: Cluster "<<(clusterIter->second)[cId] <<" T[" << std::get<1>(st)
			//			<<","<< std::get<2>(st) <<") has invalid key \n";
			//}
        }
        //flatten denseTraj
        for (auto it=denseTraj.begin();it!=denseTraj.end();++it){
            if (it->size()==0){
                auto traj = (*trajs_)[trajId];
                cv::Point2f p1 =
                     toCvPoint2f(&(traj->point(it-denseTraj.begin()))
                                                 );
                vertices.push_back(p1);
            }else{
                vertices.insert(vertices.end(),it->begin(),it->end());
            }
			//	oss <<  "|T[" << it-denseTraj.begin()<<"]|="<< std::max(1,(int)it->size())
			//	<<"\n";
        }
		for (size_t k=0;k<success.size();++k){
			if (success[k]){
				num_success_subtraj ++;
			}
		}
		rate = num_success_subtraj/float(success.size());
    }
if (trajId %100 ==0){
    qDebug()<<"\ncompleted trajectory "<<trajId <<" of length "<< vertices.size()
			<<", completion rate: " << rate;
}
	//qDebug() << "| length of sparse trajectory: " <<trajLen;
	//qDebug() << "| # of subtrajectories: " << nSubTraj;

	//qDebug()<<  "| Subtrajectory completion: " ;
	//qDebug() << oss.str().c_str();
    return std::make_pair(vertices, rate==1);
}

void TrajectoryCompletion::exportDenseTrajectories(const std::string & basename){
	std::ostringstream oss;
	oss <<basename <<".txt";
	std::string traj_fname(oss.str());
	oss.str("");
	oss <<basename <<".proj"; //projection and meta data to junction graph
	std::string proj_fname(oss.str());
	
	// write dense trajectories
	std::ofstream traj_file(traj_fname,std::ofstream::out);
	int num_completed = 0;
	traj_file << trajs_->size() <<std::endl;
	int n_total = 0;
	
	for (size_t i=0;i<trajs_->size();++i){
		if (denseTrajs_.find(i) == denseTrajs_.end()){
			traj_file << 0 << std::endl;
		}else{
			auto it = denseTrajs_.find(i);
			traj_file << it->second.size();
			for (auto pi = it->second.begin(); pi!= it->second.end();++pi){
                cv::Point2f xy = fromCvPoint2f( *pi);
                traj_file << " "<< xy.x <<" "<<xy.y;
			}
			traj_file << std::endl;
			num_completed++;
		}
		n_total++;
		if (n_total > N_TEST_TRAJS){
			break;
		}
	}
	traj_file.close();
	qDebug() << "*********";
	qDebug() << "Wrote "<< denseTrajs_.size() << " dense trajectories to "
			 << traj_fname.c_str();
	float rate = num_completed/(float)n_total; //trajs_->size(); 
	qDebug() << "Projection rate: " << rate ;
	qDebug() << "Wrote projection information to "  << proj_fname.c_str();
	qDebug() << "*********";
	qDebug() <<"nTotal = " << n_total << endl;
}


void TrajectoryCompletion::completeAll(){
    denseTrajs_.clear();
	int num_success_completion;
	for (auto it=trajClusterMap_.begin();it!=trajClusterMap_.end();++it){
		// process a single trajectory
        int trajId = it->first;
		bool success;
		
        tie(denseTrajs_[trajId], success)  = completeTrajectory(trajId);
		if (success) num_success_completion++;
	}
	qDebug() << "# of success completion: " << num_success_completion
			 <<"(" << num_success_completion/(float)trajClusterMap_.size()
			 <<")";

	//exportClusterGraphCSV( );
	//debugDenseTrajectories();	
}

void TrajectoryCompletion::completeAll2(const char* fname){
	// filter trajectories
	// #
	/*
	int maxTrajId = 0;
	for (auto it=trajClusterMap_.begin();it!=trajClusterMap_.end();++it){
		// process a single trajectory
        int trajId = it->first;
		maxTrajId = std::max(maxTrajId, trajId);
		//	qDebug() <<"trajId: " << it->first<< endl;
		}
	qDebug() << "max traj id: " << maxTrajId <<", size of trajCLusterMap: " << trajClusterMap_.size() << endl;
	*/
	std::ofstream file(fname, std::ofstream::out);
	if ( !file.is_open()){
		qDebug() << "Error: can not open file "<< fname << endl;
		return;
	}
    denseTrajs_.clear();
	file << "traj_id,input_size,sucess " << std::endl;
	int num_success_completion=0;
	//int trajStart =  0, //21,
	//	trajEnd = min(N_TRAJS, trajClusterMap_.size()-1);// 48;
	int total_trajs=0;
	for (auto it=trajClusterMap_.begin();it!=trajClusterMap_.end();++it){
		
		// process a single trajectory
		int trajId = it->first;
		//		if (trajClusterMap_.find( trajId )==trajClusterMap_.end())
		//	continue;
		bool success;
		
		tie(denseTrajs_[trajId], success)  = completeTrajectory(trajId);
		file <<trajId<<","<<(*trajs_)[trajId]->point_size() <<","<<(int)success  << std::endl;
		total_trajs++;
		if (success) num_success_completion++;
		if (total_trajs>= N_TEST_TRAJS){
			break;
		}
	}
	
	qDebug() << "# of success completion: " << num_success_completion
			 <<"(" << num_success_completion/(float)(total_trajs) //trajClusterMap_.size()
			 <<")";
	file.close();
	qDebug() <<"wrote trajectory completion results to " << fname<< endl;
	//exportClusterGraphCSV( );
	
	//debugDenseTrajectories();	
	
}
 void TrajectoryCompletion::exportClusterGraphCSV( ){
    std::ofstream out("clustergraph.txt",std::ofstream::out);

	//int kMax = Parameters::getInstance()->cluster_graph_deg;
	// assume kMax = 4
	int kMax = 4;
	out <<"cluster_id, vertex_id, px, py, n1,n2,n3,n4" <<std::endl; 
	int nClusters = cluster_graphs_.size();
    qDebug() << "export " << cluster_graphs_.size() << " cluster graphs";
     for (int i=0;i<cluster_graphs_.size();++i){
         CGraph g =cluster_graphs_[i];
         vertex_iterator vit,vit_end;
		 int nVertices = num_vertices(g);
		 
         for ( tie(vit,vit_end) = vertices(g); vit!=vit_end;++vit){
             std::shared_ptr<GpsTraj> traj = (*trajs_)[g[*vit].trajId];
             cv::Point2f p =toCvPoint2f( & (traj->point(g[*vit].pId)));
             out << i <<"," <<  g[*vit].id << ","
				 << p.x <<","<<p.y ;
             adjacency_iterator ait,ait_end;
             for (tie(ait,ait_end)=adjacent_vertices(*vit,g);ait!=ait_end;++ait){
                 out << ","<< g[*ait].id ;
             }
			 //pad -1 if cluster vertex has less than 4 neighbors
			 for (int k=0;k<kMax-out_degree(*vit,g);++k){
				 out <<",-1";
			 }
			 
             out<<"\n";

         }
     }
     out.close();
 }
void TrajectoryCompletion::exportClusterGraph( ){
    std::ofstream out("clustergraph.txt",std::ofstream::out);
    out << cluster_graphs_.size()<<"\n";
    qDebug() << "export " << cluster_graphs_.size() << " cluster graphs";
     for (int i=0;i<cluster_graphs_.size();++i){
         CGraph g =cluster_graphs_[i];
         vertex_iterator vit,vit_end;
         out <<num_vertices(g)<< "\n";
         qDebug() << "graph " << i <<" contains " << num_vertices(g) <<" vertices.";
         for ( tie(vit,vit_end) = vertices(g); vit!=vit_end;++vit){
             std::shared_ptr<GpsTraj> traj = (*trajs_)[g[*vit].trajId];
             cv::Point2f p =toCvPoint2f( & (traj->point(g[*vit].pId)));
             out << g[*vit].id << " " << p.x <<" "<<p.y <<" " << out_degree(*vit,g);
             adjacency_iterator ait,ait_end;
             for (tie(ait,ait_end)=adjacent_vertices(*vit,g);ait!=ait_end;++ait){
                 out << " "<< g[*ait].id ;
             }
             out<<"\n";

         }
     }
     out.close();
 }
void TrajectoryCompletion::debugDenseTrajectories(){
	for (auto it=denseTrajs_.begin();it!=denseTrajs_.end();++it){
		qDebug() <<"Trajectory" << it->first <<" (old len = "
				 << (*trajs_)[it->first]->point_size()  << ", new len ="
				 << it->second.size()
				 <<")" ;
		std::ostringstream oss;
		for (auto pi = it->second.begin(); pi!= it->second.end();++pi){
			oss << " ("<<pi->x <<", "<<pi->y<<")"; 
		}
		qDebug() << oss.str().c_str();
	}
}

//fit subtrajectory to cluster using simple projection
void TrajectoryCompletion::fitSubTrajToClusterSP(subtraj_descriptor st, int cid,
												 std::vector<std::vector<cv::Point2f> > &denseTraj,
                                                 std::vector<bool> &success){
	if (cluster_graphs_.find(cid)==cluster_graphs_.end()){
        makeClusterGraph(cid);
	}
	auto g = cluster_graphs_[cid];
    auto vIndex = cluster_graph_indices_[cid];
    //int trajId = std::get<0>(st);
    //std::shared_ptr<GpsTraj> queryTraj =(*trajs_)[trajId];

	// get subtraj index within cluster cid
	std::array<int,2> cKey{jEdges_[cid].first,jEdges_[cid].second};
	auto cluster = (*clusters_)[cKey];
	auto iter = find( cluster.begin(), cluster.end(),st);
   
	size_t tid = iter- cluster.begin();
	for (size_t i=0;i<vIndex[tid].size()-1;++i){
		
		
		vertex_descriptor src =vertex(vIndex[tid][i],g),
			dst = vertex(vIndex[tid][i+1],g);
		
		//	qlist.push_back(point(p.x,p.y));
		if (success[i]){// && denseTraj[i].size()>2){
			// subtraj [j,j+1] is already completed
			// alternatively, complete subtrajectory regardlessly, then
			// take the longer, better path 
			//qDebug()<<"Skip traj point "<< i <<" already completed len="<<denseTraj[i].size();
			continue;
		}
		
		// find shortest path from vIndex[i] to vIndex[i+1]

		std::vector<float> distances(num_vertices(g));
        std::vector<vertex_descriptor> parents(num_vertices(g));
		auto dmap = make_iterator_property_map(distances.begin(),
											   get(vertex_index, g));
        auto pmap = make_iterator_property_map(parents.begin(),get(vertex_index,g));
                                               //get(vertex_index, g));
		auto wmap = get(&ClusterGraphEdge::length,  g);
		
		dijkstra_shortest_paths(g,src,
					weight_map(wmap).predecessor_map(pmap).distance_map(dmap));


        if (parents[dst] == -1|| parents[dst] == dst|| distances[dst]>1e6){
			//qDebug()  << "Can not find path from vertex " << src <<" to "
			//		 <<" vertex " <<  dst;
			success[i] = false;
			continue;
        }/*else if (parents[dst]==src){
            qDebug() <<"direct trajectory in cluster"<<cid<<" src = "
                       << g[src].id <<" dst = " << g[dst].id ;
            edge_descriptor ep;bool ex=false;

            tie(ep,ex) = edge(src,dst,g);
            if (ex ){//&& g[ep].length > 3){!! quick HACK!
                float L = g[ep].length;
                g[ep].length =1000.f;

            success[i]=false;
            //recursive!
            qDebug() <<"recursively finding shortest path in cluster"<<cid<<" src = "
                       << g[src].id <<" dst = " << g[dst].id <<" len= "<< L;
            fitSubTrajToClusterSP(st, cid, denseTraj, success);
            }
            //edge is too long! we set the weight of this edge to extremely high
			}*/
		success[i] = true;
		//std::vector<cv::Point2f> path;
		vertex_descriptor prev=dst;
		denseTraj[i].clear();
        //std::ostringstream oss;oss<< g[prev].id;
        while(prev!=src){
			auto traj = (*trajs_)[g[prev].trajId];
            cv::Point2f p=  toCvPoint2f(&(traj->point(g[prev].pId)) );
			denseTraj[i].insert(denseTraj[i].begin(), p);
			prev = parents[prev];
			//  oss<<", " << g[prev].id ;
		}	
        //qDebug() <<"path from "<< g[src].id <<" to " << g[dst].id <<" dist= "<< distances[dst];
        //qDebug()<<oss.str().c_str();
	}
}
/*
void TrajectoryCompletion::makeClusterGraph(int cid){// CGraph &g, ){

	CGraph g;
	std::vector<std::vector<int> > vIndex;
	//1. read position and heading in trajectory data, construct rtree
	std::array<int,2> ckey = {jEdges_[cid].first, jEdges_[cid].second};
	std::vector<subtraj_descriptor> &cluster=(*clusters_)[ckey];
	std::vector<pvalue> plist;
	std::vector<cv::Point2f> headings ;// = findStableHeadings();
	float angleThresh = Parameters::getInstance()->cluster_graph_ang_thresh;
	vIndex.resize(cluster.size());
	for (size_t i=0; i< cluster.size();++i){
		auto si = cluster[i];
		std::shared_ptr<GpsTraj>& traj = (*trajs_)[std::get<0>(si)];
		int pStart = std::get<1>(si),pEnd = std::get<2>(si);
		for (int k=pStart; k<pEnd;++k){
			vertex_descriptor v = add_vertex(g);
			g[v].id =num_vertices(g)-1;
			g[v].trajId = std::get<0>(si);
			g[v].pId = k;
			cv::Point2f p1 = JunctionNetwork::toCvPoint2f(&(traj->point(k) ),
														  min_corner_,
														  height_, resolution_
														  );
			plist.push_back(std::make_pair(point(p1.x, p1.y),g[v].id));
			vIndex[i].push_back(g[v].id);//ith subtrajectory

			cv::Point2f v1 = Common::unitDirection(traj->point(k).head());
			headings.push_back(v1);
		}
	}
	qDebug() <<"Creating trajectory cluster graph for cluster " <<cid <<" of "
			 << cluster.size() <<" subtrajectories. # points = " << plist.size();
	
	prtree rt(plist.begin(),plist.end());
	
	vertex_iterator vit,vit_end;
	
	//2. validate heading based on neighborhood
	for (tie(vit,vit_end)=vertices(g);vit!=vit_end;++vit){

		// find nearest k neighbors that are in the heading direction of v
		cv::Point2f v0 = headings[g[*vit].pId];
		pvalue p0= plist[g[*vit].id];		
		bool validHeading =  (v0.x==0 && v0.y==0)?false:true;
		//----------------------------------------------------------------------
		//find immediate neighbor vectors
		
		int k=4; 		std::vector<pvalue> knn;
		cv::Mat nHist(1,unit_directions_.size(),CV_32F,cv::Scalar(0));
		rt.query(bgi::nearest(p0.first, k), std::back_inserter(knn));
		for (int i=0; i<knn.size();++i){
			if (knn[i].second == p0.second) continue;
			nHist += headingHistogram( headings[knn[i].second]);
		}
		nHist*= 1/float(knn.size()-1); 
		if (!validHeading){// || nHist.dot( headingHistogram(v0)) < 0.5){
			qDebug() << "current heading " << v0.x <<"," << v0.y
					 <<"is inconsistent with neighbors ("
					 <<  nHist.dot( headingHistogram(v0)) ;
			double minVal,maxVal; cv::Point minLoc,maxLoc;
			cv::minMaxLoc(nHist,&minVal,&maxVal, &minLoc, &maxLoc);
			v0 = unit_directions_[maxLoc.x];
			qDebug() << "new heading= " << v0.x <<","
					 <<v0.y <<" new dot product=" << nHist.dot(headingHistogram(v0));

			if (nHist.dot(headingHistogram(v0))>0.8){
				headings[g[*vit].pId] = v0;
				qDebug() <<"use new heading ";
			}

					 
		}else{
			qDebug() << "current heading " << v0.x <<"," << v0.y
					 <<"is consistent with neighbors ("
					 << nHist.dot( headingHistogram(v0));
		}
	}
	int iterMax = Parameters::getInstance()->cluster_graph_deg;
	for  (tie(vit,vit_end)=vertices(g);vit!=vit_end;++vit){
		int iter =0;
		cv::Point2f v0 = headings[g[*vit].pId];
		bool validHeading =  (v0.x==0 && v0.y==0)?false:true;
		pvalue p0= plist[g[*vit].id];		
		cv::Point2f x0(p0.first.get<0>(),p0.first.get<1>());//position of edge src
		auto traj = (*trajs_)[g[*vit].trajId];
		
		prtree::const_query_iterator nit;
		int K = std::min(50, (int)plist.size()-1);
		for (nit=rt.qbegin(bgi::nearest(p0.first,K));nit!=rt.qend();++nit){
			int nvid = nit->second;
			if (nvid == g[*vit].id) continue;
			// candidate position of edge target
			cv::Point2f n(nit->first.get<0>(),nit->first.get<1>()); 
			if (  Common::cosAngleBetween(n-x0, v0) > angleThresh){
				// add an edge
				vertex_descriptor nvi = vertex(nit->second,g);
				edge_descriptor e; bool tmp;
				tie(e,tmp) = add_edge(*vit, nvi,g);
				g[e].length =cv::norm(n-x0);
				iter++;
			}else if (!validHeading){
			   	vertex_descriptor nvi = vertex(nit->second,g);
				v0 = Common::unitDirection(traj->point(g[nvi].pId).head());
				if (Common::cosAngleBetween(n-x0,v0)>angleThresh){
					edge_descriptor e; bool tmp;
					tie(e,tmp) = add_edge(*vit, nvi,g);
					g[e].length =cv::norm(n-x0)*2;
					iter++;
				}

										   
			}
			
			if (iter==iterMax)
				break;
		}
		
		//	qDebug() <<"** p " << g[vit].id;
	}
	cluster_graphs_[cid] = g;
	cluster_graph_indices_[cid]= vIndex;
	qDebug()<<"finished creating cluster graph of " << num_edges(g) <<" edges";
	}*/
/*
std::vector<cv::Point2f> TrajectoryCompletion::findStableHeading(int cid){
	// search for k nearest neighbor, construct heading pca, l1/l1+l2>0.8
	// or, discretize angles  into 24 bins, find the maximum bin. (if maximum bin is not salient,
	// estimate direction by vectors i+1 and i-1.
	std::vector<cv::Point2f> heading0;
	return heading0;
}
*/

void TrajectoryCompletion::makeClusterGraph(int cid){// CGraph &g, ){

	CGraph g;
	std::vector<std::vector<int> > vIndex;
	
	std::array<int,2> ckey = {jEdges_[cid].first, jEdges_[cid].second};
	std::vector<subtraj_descriptor> &cluster=(*clusters_)[ckey];
	std::vector<pvalue> plist;

	float angleThresh = Parameters::getInstance()->cluster_graph_ang_thresh;
	vIndex.resize(cluster.size());
	for (size_t i=0; i< cluster.size();++i){
		auto si = cluster[i];
		std::shared_ptr<GpsTraj>& traj = (*trajs_)[std::get<0>(si)];
		int pStart = std::get<1>(si),pEnd = std::get<2>(si);
		for (int k=pStart; k<pEnd;++k){
			vertex_descriptor v = add_vertex(g);
			g[v].id =num_vertices(g)-1;
			g[v].trajId = std::get<0>(si);
			g[v].pId = k;
            cv::Point2f p1 = toCvPoint2f(&(traj->point(k) ) );
			plist.push_back(std::make_pair(point(p1.x, p1.y),g[v].id));
			vIndex[i].push_back(g[v].id);//ith subtrajectory
		}
	}
	qDebug() <<"Creating trajectory cluster graph for cluster " <<cid <<" of "
			 << cluster.size() <<" subtrajectories. # points = " << plist.size();
	
	prtree rt(plist.begin(),plist.end());
	vertex_iterator vit,vit_end;

	float neighborPenalty =
		Parameters::getInstance()->cluster_graph_neighbor_penalty;
	for (tie(vit,vit_end)=vertices(g);vit!=vit_end;++vit){
		int kMax = Parameters::getInstance()->cluster_graph_deg;
		// find nearest k neighbors that are in the heading direction of v
		int k=0;

		pvalue p0= plist[g[*vit].id];		
		cv::Point2f x0(p0.first.get<0>(),p0.first.get<1>());//position of edge src
		auto traj = (*trajs_)[g[*vit].trajId];
		cv::Point2f v0 = Common::unitDirection(traj->point(g[*vit].pId).head());
		bool validHeading =  (v0.x==0 && v0.y==0)?false:true;
         if (!validHeading){
             qDebug() << "invalid heading produced by angle "<< traj->point(g[*vit].pId).head();
         }
		prtree::const_query_iterator nit;
		int K = std::min(50, (int)plist.size()-1);
		for (nit=rt.qbegin(bgi::nearest(p0.first,K));nit!=rt.qend();++nit){
			int nvid = nit->second;
			if (nvid == g[*vit].id) continue;
			// candidate position of edge target
			cv::Point2f n(nit->first.get<0>(),nit->first.get<1>()); 
			if (  Common::cosAngleBetween(n-x0, v0) > angleThresh){
				// add an edge
				vertex_descriptor nvi = vertex(nit->second,g);
				edge_descriptor e; bool tmp;
				tie(e,tmp) = add_edge(*vit, nvi,g);
				g[e].length =cv::norm(n-x0);
				k++;
			}else if (!validHeading){
			   	vertex_descriptor nvi = vertex(nit->second,g);
				if (g[nvi].pId<0 || g[nvi].pId >=traj->point_size()){
					qDebug()<<"assertion failed: invalid cluser graph vertex id "
							<< g[*vit].trajId <<":"<<g[nvi].pId;
				}else{
					v0 = Common::unitDirection(traj->point(g[nvi].pId).head());
					if (Common::cosAngleBetween(n-x0,v0)>angleThresh){
						edge_descriptor e; bool tmp;
						tie(e,tmp) = add_edge(*vit, nvi,g);
						g[e].length =cv::norm(n-x0)*neighborPenalty;
						k++;
					}
				}

										   
			}
			
			if (k==kMax)
				break;
		}
	   
	}
	cluster_graphs_[cid] = g;
	cluster_graph_indices_[cid]= vIndex;
	//	qDebug()<<"finished creating cluster graph of " << num_edges(g) <<" edges";
}
void TrajectoryCompletion::fitSubTrajToCluster(subtraj_descriptor st, int cid,
								std::vector<std::vector<cv::Point2f> > &denseTraj){
	int trajId = std::get<0>(st);
	int maxIter = 80;
	float beta = 0.2f;
	float angleThresh= Parameters::getInstance()->skeleton_proj_angle_thresh;
	//GpsTraj traj = gpsTraj_->get

	std::array<int,2> ckey = {jEdges_[cid].first, jEdges_[cid].second};
	std::vector<subtraj_descriptor> &cluster=(*clusters_)[ckey];
	std::vector<point> plist,vlist, qlist;

	std::shared_ptr<GpsTraj> queryTraj =(*trajs_)[std::get<0>(st)];
	
	for (size_t i=0; i< cluster.size();++i){
		auto si = cluster[i];
		std::shared_ptr<GpsTraj>& traj = (*trajs_)[std::get<0>(si)];
		int pStart = std::get<1>(si),pEnd = std::get<2>(si);
		for (int k=pStart; k<pEnd;++k){
            cv::Point2f p1 = toCvPoint2f(&(traj->point(k) ));//is this Needed?
			plist.push_back(point(p1.x, p1.y));
		}
	}
	Rtree rt(plist.begin(),plist.end());

	// skip if this cluster contains too few points
	if  (plist.size() < Parameters::getInstance()->skeleton_min_cluster_size)
		return;
	
	//greedy estimate dense trajectory
	//qDebug() <<"Completing sparse sub-traj of length " << qlist.size()
	//		 << " in cluster of size " << plist.size();


	int j =std::get<1>(st);
    cv::Point2f p = toCvPoint2f(&(queryTraj->point(j) ));
	for (; j<std::get<2>(st)-1;++j){
		
        cv::Point2f pNext = toCvPoint2f(&(queryTraj->point(j+1) ));

		//	qlist.push_back(point(p.x,p.y));
		if (denseTraj[j].size()>1){
			// subtraj [j,j+1] is already completed
			p = pNext;
			//qDebug()<<"skip traj point "<< j;
			continue;
		}
		denseTraj[j].clear();
		// otherwise, fill path between t[j] and t[j+1]
		// get k nearest neighbors from t[j]
		

		cv::Point2f x0 = p,v0= Common::unitDirection(queryTraj->point(j).head());
		/** (skip origional pt in traj?) **/
		denseTraj[j].push_back(x0); 
		
		//qDebug()<<"step = " << step_ ;
		int iters=0;
		while (cv::norm(x0-pNext) > step_ && iters<maxIter){
			
			Rtree::const_query_iterator nit;
			float minEnergy=Common::inf;
			cv::Point2f minP ;
			point p0( x0.x,x0.y );int kk=0;
			for (nit = rt.qbegin(bgi::nearest(p0,200));nit!=rt.qend();
				 ++nit){ //@todo: filter points in opposite direction with satisfies
				cv::Point2f n(nit->get<0>(),nit->get<1>()); //neighbor of p
				cv::Point2f v = n-x0;
				if (Common::cosAngleBetween(v, v0) <angleThresh)
					continue;
				

				if (n == x0) continue; 
				float energy = beta*cv::norm(n-pNext) + (1-beta)*(1- Common::cosAngleBetween(n-x0,v0));
				/*qDebug() <<"T["<<trajId<<"."<< j <<"] "<< kk<<"th neighbor (" << n.x<<","<<n.y <<")"
						 <<" v0=(" << v0.x <<"," << v0.y <<")"
						 <<" dist to target ="<< cv::norm(n-pNext)  <<","
						 <<" cos angle = " << Common::cosAngleBetween(n-x0,v0);
				*/
				kk++;
				if (energy < minEnergy){
					minEnergy  = energy;
					minP = n;
				}

				if (kk>3)
					break;
			}
			v0 = minP-x0;
			x0 = minP;
			if (minP.x==0 && minP.y==0){
				qDebug()<<"Error completing "<<j<<"th point! energy=" << minEnergy <<" iter=" <<iters;
				break;
			}else{
				denseTraj[j].push_back(minP);
			}
			iters++;
		}
		qDebug() << "Completed subtrajectory in "<< trajId<<" of cluster " << cid
				 <<"-" <<j<<" in "
				 << iters <<" iterations.";
	}
	
}
void TrajectoryCompletion::debugTrajClusterMap(){
	qDebug() <<"=================TrajClusterMap=====================";
	for (auto it= trajClusterMap_.begin();it!=trajClusterMap_.end();++it){
		std::ostringstream oss;
		for (size_t i=0; i< it->second.size();++i){
			std::pair<int,int> e = jEdges_[it->second[i]];
			oss << it->second[i] <<"(" << e.first<<"," << e.second << ") " ;
		
		qDebug() << "Trajectory "<<it->first << ": " << oss.str().c_str() ;
		}
	}
}

subtraj_descriptor TrajectoryCompletion::getSubTrajFromCluster(int clusterId,
															   int trajId){
	std::array<int,2> ckey = {jEdges_[clusterId].first, jEdges_[clusterId].second};
	std::vector<subtraj_descriptor> &cluster = (*clusters_)[ckey];
	subtraj_descriptor subtraj{trajId,0,0};
	
	auto iter =
	find_if(cluster.begin(),cluster.end(),
			[trajId](const subtraj_descriptor &d){return trajId == std::get<0>(d);});
	if (iter!= cluster.end()){
		subtraj = *iter; 
	}else{
		//qDebug() <<"Cannot find subtrajectory of cluster " << clusterId<<" in traj "
		//		 << trajId << " of size " << cluster.size();
	}
	return subtraj;
}

cv::Mat TrajectoryCompletion::headingHistogram(std::vector<cv::Point2f> h){
	cv::Mat hist(1,unit_directions_.size(),CV_32F,cv::Scalar(0));
	for (size_t v=0;v<h.size();++v){
		for (size_t i=0;i<hist.cols;++i){
			
			hist.at<float>(0,i) +=  std::max(0.f,unit_directions_[i].dot(h[v]));
		}
	}
	hist*=1/float(h.size());
	return hist;
}
cv::Mat TrajectoryCompletion::headingHistogram(cv::Point2f h){
	cv::Mat hist(1,unit_directions_.size(),CV_32F,cv::Scalar(0));
	for (size_t i=0;i<hist.cols;++i){
		hist.at<float>(0,i) +=  std::max(0.f,unit_directions_[i].dot(h));
	}
	cv::Mat histN;
	cv::normalize(hist,histN);
	return histN;
}
