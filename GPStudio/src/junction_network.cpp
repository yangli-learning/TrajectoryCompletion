#include "junction_network.h"
#include "l_shape_iterator.h"
#include "gps_trajectories.h"
#include "gps_trajectory.pb.h"
#include "l_shape_iterator.h"
#include "viterbi.h"
#include <utility>
#include "parameters.h"

#define _USE_MATH_DEFINES
#include <cmath>
using namespace Common;
 //----------- Definitions for ExtendedEdge -------------------

ExtendedEdge::ExtendedEdge():id(-1),headDir(cv::Point2f(0,0)),
tailDir(cv::Point2f(0,0)),headFinal(false),tailFinal(false),
  pts(std::vector<cv::Point2f>()),nNeighbors(0),invalid(false) {
    // default constructor
}


ExtendedEdge::ExtendedEdge(int _id, cv::Point2f _headDir,cv::Point2f _tailDir,
                           bool _headFinal, bool _tailFinal,
                           std::vector<cv::Point2f> _pts):
    id(_id),headDir(_headDir),tailDir(_tailDir),headFinal(_headFinal),
    tailFinal(_tailFinal),  pts(_pts) ,nNeighbors(0),invalid(false){
	/*    qDebug() <<"initialize edge of id " << id <<" of length " << pts.size()
               <<" headDir = " << headDir.x<<"," << headDir.y << " tailDir = "
			 << tailDir.x << "," << tailDir.y;
	*/
    // initialization constructor
}

ExtendedEdge::~ExtendedEdge(){
    // destructor does nothing
    pts.clear();
}
 ExtendedEdge::ExtendedEdge(const ExtendedEdge &ee):id(ee.id),
headDir(ee.headDir),tailDir(ee.tailDir),headFinal(ee.headFinal),
tailFinal(ee.tailFinal),pts(ee.pts){
    qDebug() <<"Copy constructor is called on edge " << ee.id <<" of input length "<< ee.pts.size()
               <<" stored length " << pts.size();

}
ExtendedEdge & ExtendedEdge::operator=(const ExtendedEdge &ee){
    qDebug()   <<"Assignment operator = is called on edge " << ee.id;
    id = ee.id;
    headDir =  ee.headDir;
    tailDir = ee.tailDir;
    headFinal = ee.headFinal;
    tailFinal = ee.tailFinal;
    pts = ee.pts;
    return *this;
}
//---------- Definitions for JunctionNetwork -----------------
JunctionNetwork::JunctionNetwork(L1SkeletonGraph *g, L1Skeleton_ *sk,
								 float h, float w):
 graph_(g),height_(h),width_(w),skeleton_(sk){

    step_ = Parameters::getInstance()->skeleton_extension_step*
              Parameters::getInstance()->skeleton_sample_cell_size;//0.5*
	min_edge_angle_  = Parameters::getInstance()->skeleton_branch_min_angle;
	min_edge_length_ = Parameters::getInstance()->skeleton_min_edge_len/\
		float(Parameters::getInstance()->map_image_resolution);
	max_missing_edge_length_ = Parameters::getInstance()->\
		skeleton_max_edge_len/float(Parameters::getInstance()->map_image_resolution);
    cv::Mat density = sk->getDensityMatrix();
    double minV,maxV;cv::Scalar mean,stdev;
    cv::minMaxLoc(density,&minV, &maxV);
    cv::meanStdDev(density, mean,stdev,cv::Mat());
    min_vertex_density_ = minV; //  (mean.val[0]- minV)*0.1+minV;
    qDebug() <<"density threshold = " << min_vertex_density_<<" min = "
            << minV<<" max = "
            << maxV << " mean = "
            << mean.val[0] <<  " stdev = " << stdev.val[0];
    process(); //process l1skeleton graph (extend edges and compute intersections)
}

void JunctionNetwork::initJunctionGraph(){
    jgraph_.init( graph_->getGraph());
	//jgraph_.writeDotGraph("junction_graph_00.dot",height_);
}

void JunctionNetwork::process(){
    bool finished = false;
    int maxIter =Parameters::getInstance()->skeleton_extension_max_iter;
    float step = step_;
    int i=0;
    while(!finished && i < maxIter){
        qDebug() <<"=========== iter " << i <<" ===============";
 
        initExtendedEdges();
        finished = growExtendedEdges(step);
        processIntersections();
 
        std::ostringstream oss;
        oss<<"graph"<<i<<".dot";
        graph_->writeDotGraph(oss.str().c_str());

        i++;
    }

    //remove internal nodes again after computing intersection
    graph_->removeInternalNodes();
	graph_->writeDotGraph("graph_clean0.dot");

    clusterEndpoints(Parameters::getInstance()->skeleton_cluster_radius  *min_edge_length_);
    graph_->writeDotGraph("graph_clean1.dot");
	

	graph_->collapseEdgesShorterThan(min_edge_length_);//found invalid edge
	graph_->writeDotGraph("graph_clean2.dot");
	

	graph_->mergeSimilarBranches(0.9*min_edge_length_, min_edge_angle_);//no invalid edge
	graph_->writeDotGraph("graph_clean3.dot");
	
	graph_->removeInternalNodes();//found more invalid edge
	graph_->writeDotGraph("graph_clean4.dot");
}

void JunctionNetwork::initExtendedEdges(){
    extended_edges_.clear();
    L1SkeletonGraph::SGraph g = graph_->getGraph();
    L1SkeletonGraph::un_edge_iterator e, e_end;
    int i=0;
    for (boost::tie(e,e_end) = boost::edges(g); e != e_end;++e){
        std::vector<cv::Point2f>  pts = g[*e].internalPts;
        //if (pts.size() < 2)
        //    continue;
       // ExtendedEdge ee;
        L1SkeletonGraph::vertex_descriptor head = boost::target(*e,g),
                tail = boost::source(*e,g);
        bool headFinal = g[head].type != SkeletonGraphVertexType::TERMINAL,
                tailFinal = g[tail].type != SkeletonGraphVertexType::TERMINAL;

        /*qDebug() << "== processing edge " << i << " (" << g[head].id<<"-"
                 << g[tail].id <<") of size " << pts.size() <<"===";*/
        // initialize pt
        cv::Point2f headDir = pts[pts.size()-1]- pts[pts.size()-2],
                tailDir =  pts[0]-pts[1];
        normalize(headDir,headDir);
        normalize(tailDir,tailDir);
        shared_ptr<ExtendedEdge> ee(new ExtendedEdge( i,
                                                 headDir,tailDir,
                                                      headFinal,
                                                      tailFinal,
                                                      pts));
        if (pts.size() < 2){
            ee->headFinal= ee->tailFinal=true;
        }
        extended_edges_.push_back(ee);
        i++;
    }
}
bool JunctionNetwork::extendEdge(shared_ptr<ExtendedEdge> ee, EndpointType type,
								 float len){
    bool isHead = type==EndpointType::HEAD;
    bool success =true;
    int n = ee->pts.size();
    cv::Point2f &dir = (isHead)? ee->headDir:ee->tailDir;
    cv::Point2f base = (isHead)? ee->pts[n-1]:ee->pts[0];//use longer offset
    cv::Point2f newDir =   updateEndpointDirection( base , dir);
	cv::Point2f edgeDirection = (isHead)?(ee->pts[n-1]-ee->pts[0]):
		(ee->pts[0]-ee->pts[n-1]);
    if ( cosAngleBetween(edgeDirection,newDir) >=0.80){ //~30 degrees
         dir = newDir;
    }else{
		//	dir = 0.5*dir+0.5*edgeDirection;
		len =len*0.5;
	}
    cv::Point2f newEndpoint = base + dir*len;
   // L1SkeletonGraph::SGraph g = graph_->getGraph();
    //edge_descriptor ed = graph_->getEdgeById(ee->id).first ;

    //!!!! problem here !!!!!! updateEdgeGeometry may not be called or incorrect
    if (  newEndpoint!= base && inBBox(newEndpoint)
		  && isValidExtension(newEndpoint) ){
        int oldLen = ee->pts.size();
        if (isHead){
            ee->pts.push_back(newEndpoint);
            //update head.pos!
            //vertex_descriptor vd= target(ed,g);
           // g[vd].p.pos = newEndpoint;
           // graph_->moveVertex(vd,newEndpoint);
        }else{
            ee->pts.insert(ee->pts.begin(), newEndpoint);
            //update tail.pos!
        }
        std::string headstr = (isHead)?"head":"tail";
        //qDebug()<<"extend edge "<< ee->id << " :  "<< headstr.c_str()
        //         <<" from len " <<oldLen << " to " << ee->pts.size();
        graph_->updateEdgeGeometry( ee->id, ee->pts );
    }else{
		if (newEndpoint==base){
			qDebug()<<"Can not compute extending direction";
		}
		/*if (!isValidExtension(newEndpoint)){
			qDebug()<<"Not valid extension.";
			}*/
        if (isHead){
            ee->headFinal = true;
        }else{
            ee->tailFinal = true;
        }
        success = false;
    }
    return success;
}
bool JunctionNetwork::inBBox(cv::Point2f p){
    float padding = std::min(height_,width_)*0.15;
    float minY = -padding, maxY = height_+padding, minX = -padding, maxX = width_+padding;

    return (p.x>minX && p.x < maxX && p.y > minY && p.y < maxY);
}

bool JunctionNetwork::growExtendedEdges(float len){
    // check for intersection, *safely** extend edges not effected by interesections
    bool finished = true;
    checkIntersections();
    for (auto it = extended_edges_.begin(); it!= extended_edges_.end();++it ){
        shared_ptr<ExtendedEdge> ee = *it;
        size_t n = ee->pts.size();
        if (!ee->headFinal ){
           if (  extendEdge(ee, HEAD,len))
               finished = false;
        }
        if (!ee->tailFinal){
           if ( extendEdge(ee,TAIL,len))
               finished =false;
        }


    }
    return finished;

}

cv::Point2f JunctionNetwork::updateEndpointDirection(cv::Point2f p0, cv::Point2f oldDir){
    cv::Mat dist;
    int nDirections;

    if (Parameters::getInstance()->skeleton_use_heading){
        dist = skeleton_->getHeadingDistribution(p0.x , p0.y);
        nDirections = dist.cols;
    }else{
        cv::Mat distRaw = skeleton_->getHeadingDistribution(p0.x, p0.y);
        nDirections =  distRaw.cols;///2; //assume nDirections is even
        dist = distRaw;
    }
    if (nDirections==0){
        return cv::Point2f(0,0);
    }
    double minV,maxV; cv::Point maxLoc, minLoc;
    cv::minMaxLoc(dist, &minV, &maxV, &minLoc, &maxLoc);
    cv::Point2f p =
        LShapeIterator::UnitDirection((nDirections-maxLoc.x)*float(360/(2*nDirections)));
    p.y = -p.y;
    if (p.dot(oldDir) < 0 ){
        p = -p;
    }
    //qDebug() <<"update endpoint " << p0.x << "," << p0.y << " direction = " <<p.x <<"," <<p.y;
    return p;
}

bool JunctionNetwork::isValidExtension(cv::Point2f newHead){
    // Assume heading is NOT normalized
    float density = skeleton_->getDensity(newHead.x, newHead.y, step_*3);
   // qDebug() <<"Check is valid extension.. density=" << density;
    //return density>0;// >= min_vertex_density_;
    return true;
}

void JunctionNetwork::initializeSegmentRTree(std::vector<value> &queries,
                                             rtree &rt){
    std::vector<value> segments;
    for (std::vector<shared_ptr<ExtendedEdge> >::iterator it= extended_edges_.begin() ;
         it != extended_edges_.end(); ++it){
        shared_ptr<ExtendedEdge> ee = *it;
        std::vector<cv::Point2f> pts =  ee->pts;
        size_t edgeId = it-extended_edges_.begin();


        for (size_t j=0; j < pts.size()-1; ++j){
            cv::Point2f &v1 = pts[j];
            cv::Point2f &v2 = pts[j+1];
            point p1(v1.x,v1.y), p2(v2.x,v2.y);
            segment seg( p1, p2 );
            value v = std::make_tuple(seg,edgeId,j);
            segments.push_back(v); //std::make_tuple(seg,edgeId,j));

            if ( (j==0 && !ee->tailFinal)||
                 (j == pts.size()-2 && ! ee->headFinal)){
                queries.push_back(v);
            }
        }
    }
    qDebug() << segments.size() << " search segments. " << queries.size()
             <<" query segments";
    rt = rtree(segments.begin(), segments.end());
}

void JunctionNetwork::checkIntersections(){
    // check upcoming intersection between edge end segments:

    intersections_.clear();
    std::vector<value> queries;
    rtree rt;
    initializeSegmentRTree( queries, rt);
    float r = step_*1.5;
    qDebug() <<"Checking intersection with r = "<< r;
    for (auto qit = queries.begin(); qit !=queries.end();++qit){
        // project query (qSeg) endpoint p to its neighbors
        // find intersections of type I, II and III.

        size_t qEdgeId = std::get<1>(*qit);
        shared_ptr<ExtendedEdge> qe= extended_edges_[qEdgeId];

        // do not check intersection if query edge has been invalidated
        if (qe->invalid){
            continue;
        }
        EndpointType qType = (get<2>(*qit)==0)?EndpointType::TAIL:
                                                 EndpointType::HEAD;


        //find nearest 10 neighbors to query segment
        rtree::const_query_iterator nit;
        segment query=std::get<0>(*qit);
        int i=0;
        for (nit=rt.qbegin(bgi::nearest( query,10)); nit!= rt.qend();
                 ++nit){
            size_t nEdgeId = std::get<1>(*nit);
            shared_ptr<ExtendedEdge> ne= extended_edges_[nEdgeId];

            if ( (!ne->invalid) && qEdgeId != nEdgeId){
                segment qSeg = std::get<0>(*qit), nSeg = std::get<0>(*nit);

                // check wether qEdge and nEdge are incident
                if (graph_->incidentEdges(qEdgeId, nEdgeId)){
                    continue;
                }

                double d = bg::distance( qSeg, nSeg);
                if (d>r){  //<------Need to check the distance function used in bg
                    break;
                }
                // compute projection
                cv::Point2f p  = (qType == EndpointType::HEAD)?
                             boostToCVPoint(qSeg.second):
                             boostToCVPoint(qSeg.first);
                cv::Point2f e1 = boostToCVPoint(nSeg.first),
                        e2 = boostToCVPoint(nSeg.second);
                cv::Point2f newLoc;
                bool onEdge = projectPointOnSegment(p,e1,e2, newLoc);
                float projDist = cv::norm(newLoc- p); //distance from p to its projection
                //qDebug()<<"Debug: distance = " << d <<" projection len = " << projDist;
			   
                //if ( projDist > r)
                //    newLoc = 0.5*e1+0.5*e2;
                // <-- determine intersection type
                int lastSegId = ne->pts.size()-2;
                if  (std::get<2>(*nit) ==0 || std::get<2>(*nit)== lastSegId){
                    //type I and type II intersection
                    qe->headFinal = true;
                    qe->tailFinal = true;
                    qe->invalid = true;
                    if (std::get<2>(*nit) ==0 ){
                        ne->tailFinal =true;
                    }else{
                        ne->headFinal  = true;
                    }
                }else{
                    ne->headFinal = true;
                    ne->tailFinal = true;
                    ne->invalid = true;
                    if (qType ==EndpointType::TAIL ){
                        ne->tailFinal = true;
                    }else{
                        ne->headFinal = true;
                    }
                    // type III intersection
                }
                //if (addIntersection && ne->nNeighbors ==0){//!>TEMP
				shared_ptr<EdgeIntersection> ei(new EdgeIntersection( *qit,
																	  *nit,
																	  qType,
																	  newLoc));
				intersections_.push_back(ei);
                //}
                /*if (endType == EndpointType::HEAD){
                    qe->headFinal = true;
                }else{
                    qe->tailFinal = true;
                } */
                ne->nNeighbors ++;
                break; //<----- alternatively, store a set of nEdgeIds, and only add
                //intersections of unique edge id
            }
        }
    }

	//qDebug()<<"FINISHed checking intersections....";
}

void JunctionNetwork::computeSubTrajClusters(){
	//loop over all edges, add traj to clusters
	jgraph_.getSubTrajClusters(clusters_);
}

void JunctionNetwork::validateTopologyByGPS(
                    osg::observer_ptr<GPSTrajectories>  trajectories,
                    osg::Vec3 min_corner,  float resolution,float minFlow){

    GPSTrajectories::GpsTrajs gps_trajs = trajectories->getGpsTrajs();
    qDebug() <<"validate through " << gps_trajs.size() << "trajectories";

    rtree rt = makeRTreeFromEdges();
    //!!!!!!!!!!!!!!!!<<<<<<
    int shortTraj=0, ambTraj=0 ;
    size_t maxTraj = Parameters::getInstance()->max_n_trajectories;
    int ct_all=0, ct_discarded=0;
    int N=5; //size of nearest neighborhood (candidate search)
	float maxR =30;
    //store candidate weights
	bool useViterbi = false;
    qDebug() << "Built rtree of size " << rt.size();
	size_t i,i_end;
	n_direct_votes_=0; n_indirect_votes_=0;
    for  ( i = 0, i_end = gps_trajs.size(); i < i_end; ++i) {
        //int seg_prev = -1;
        std::shared_ptr<GpsTraj> gps_traj = gps_trajs[i];
        if (gps_traj->point_size() < 2){
            shortTraj ++;
            continue;
        }
        size_t nPoints =gps_traj->point_size();
        cv::Mat distances(nPoints,N,CV_32F, cv::Scalar(inf));
        cv::Mat candidateSegs(nPoints,N,CV_32S, cv::Scalar(-1));
        cv::Mat candidatePos(nPoints,N,CV_32S, cv::Scalar(-1));
		
		//qDebug() <<"\n======= Project Trajectory " <<i<<" of length "
        //       << gps_traj->point_size() <<" =========\n";

        for (size_t j = 0, j_end = nPoints; j < j_end; ++j) {
            ct_all++;
            point p = toBoostPoint(&gps_traj->point(j), min_corner, resolution);

            std::vector<value> result_n;
            rt.query(bgi::nearest(p,2*N) , std::back_inserter(result_n));

            sort(result_n.begin(),result_n.end(),[p](const value& v1, const value &v2){
                    return bg::distance(p,get<0>(v1)) < bg::distance(p,get<0>(v2));
                });
			//std::ostringstream oss;
			  //oss <<" KNN of t" <<i<<"."<<j<<": ";
			int ci=0; std::set<int> edgeIds;
            for (size_t k=0;k<result_n.size();++k){
				float d = bg::distance(p,get<0>(result_n[k]));
				if (d < maxR || ci >= N){
					if ( edgeIds.find(get<1>(result_n[k]) )== edgeIds.end() ){
						edgeIds.insert(get<1>(result_n[k]));
						distances.at<float>(j,ci) = d;

						candidateSegs.at<int>(j,ci) = get<1>(result_n[k]);
						candidatePos.at<int>(j,ci) = get<2>(result_n[k]);
						ci++;
					}
				}else{
					break;
				}
				//  oss << k<<  "-" <<  get<1>(result_n[k]) <<"|"
				//   << get<2>(result_n[k])<< "("<< distances.at<float>(j,k) <<") ";
            }
			// oss <<"\n";
			// qDebug()<< oss.str().c_str() ;
            // skip if nearest distance is greater than a threshold

        }

		if (useViterbi ){
        // project gps trace to graph
			auto projection = projectGPSToGraph(gps_traj, candidateSegs,
											candidatePos, distances, min_corner);
			if (projection.size()>0){
				processVotes(i, gps_traj, candidateSegs, candidatePos, projection);
			}
		}else{
			voteJunctionGraphNaive(i, gps_traj,candidateSegs,candidatePos,
								   distances,min_corner,maxR);
		}
        if (i>maxTraj) break;
    }
    //qDebug() << "[vote] # of direct votes: "<< n_direct_votes_<<" , # of indirect votes: "
    //		 << n_indirect_votes_;
	qDebug() <<"Read " << i << " trajectories. ";
	// remove low vote edges
	//	if (minFlow !=-1){
	jgraph_.setMinFlowSize(minFlow);
		//}
	jgraph_.filterEdgeByVotes();
	
    //qDebug() << "Read "<< shortTraj<<" singletons, "
	//        << ambTraj <<" ambiguous trajectories. "<< ct_discarded <<" are discarded.";
	jgraph_.writeDotGraph("junction_graph_0.dot",height_);
	computeSubTrajClusters();
}
void JunctionNetwork::voteJunctionGraphNaive(int trajId,
		 std::shared_ptr<GpsTraj> traj,
         const cv::Mat &candidateSegs,  const cv::Mat &candidatePos,
		 const cv::Mat &distances, const osg::Vec3 &min_corner, int maxR){

	L1SkeletonGraph::SGraph g=graph_->getGraph();
    size_t nPoints= candidateSegs.rows;
	size_t maxNumCandidates = candidateSegs.cols;
	
	float angle_thresh = Parameters::getInstance()->skeleton_proj_angle_thresh;
	float ratio_thresh = Parameters::getInstance()->skeleton_proj_ratio_thresh;
    std::vector<TrajProjection > projections;
	
	int ambTraj=0, ct_discarded=0;
    for (size_t i=0; i<nPoints;++i){
		int nCand = 0;
        while (nCand < maxNumCandidates) {
			int edgeId = candidateSegs.at<int>(i,nCand),
				segId  = candidatePos.at<int>(i,nCand);
            if (edgeId == -1 || segId ==-1)
                break;
			nCand++;
		}
		if (nCand == 0)
			continue;
		
        float variance = rowVariance<float>(distances.colRange(0,nCand),i,0.f);
		float confidence = 1-sqrt(variance)/maxR;

		if (confidence >ratio_thresh){
			int edgeId = candidateSegs.at<int>(i,0),
				segId  = candidatePos.at<int>(i,0);			

			edge_descriptor e = graph_->getEdgeById(edgeId).first;
            std::vector<cv::Point2f> shape = g[e].internalPts;
            cv::Point2f p1 = shape[segId], p2 = shape[segId+1];
            cv::Point2f p = toCvPoint2f(&(traj->point(i)), min_corner,true);
			cv::Point2f proj;
            projectPointOnSegment(p,p1,p2 ,proj);

			float distToHead=0,distToTail=0;

			if (i>0){
				cv::Point2f p_prev=toCvPoint2f(&(traj->point(i-1)),min_corner,true);
				distToTail =cv::norm(p_prev- p1);//distance from src-->p_prev
				distToHead = cv::norm(p_prev-p2); //distance from target-->p_prev
			}
			int ang_curr;
			try{
				ang_curr = (int32_t)(traj->point(i).head()); 
			}catch(...){
				qDebug()<<"Error: bad angle!";
				ct_discarded++;
				continue; //skip point with bad angle
			}
			//traj direction
			cv::Point2f tDir = LShapeIterator::UnitDirection(ang_curr);
			cv::Point2f sDir = p2-p1;  //segment direction

			bool isForward ; 
			if (Common::cosAngleBetween(tDir,sDir) > angle_thresh) { 
				isForward =true; 
			}else if (tDir.dot(sDir)< -angle_thresh) {
				isForward = false;
      
			}else{
				ambTraj++;
				ct_discarded++;
				continue;
			}	
			if (confidence>ratio_thresh){
				projections.push_back(makeTrajProjection(i, e, proj,
													 isForward));
			}
			
		}
    
	}
		// compute breakpoint
	if (projections.size()>1){
		std::vector<int> breakpts;
		breakpts.push_back(projections[0].pid);
		edge_descriptor  curr, prev = projections[0].edge;
	
		for (size_t i=1; i< projections.size(); ++i){
			curr = projections[i].edge;
			if  (curr != prev){
				breakpts.push_back(projections[i].pid);
				prev = curr;
			}
		}
		breakpts.push_back(projections.back().pid);

		// apply vote
		vote(trajId,projections, breakpts,traj->point_size());
	}
}

void JunctionNetwork::vote(int trajId, std::vector<TrajProjection> projections,
						   std::vector<int> breakpts, int n){

	L1SkeletonGraph::SGraph g = graph_->getGraph();
	int vote_ct =0,n_filtered=0; 
	edge_descriptor curr, prev = projections[0].edge;
	//std::vector<std::vector<int> > allvotes;
	for (size_t i=1;i<projections.size();++i){
		curr = projections[i].edge;

		if (curr !=prev){
			vertex_descriptor s1,t1,s2,t2;
			if (projections[i-1].isForward){
				s1 = source(prev,g);
				t1 = target(prev,g);
			}else{
				s1 = target(prev,g);
				t1 = source(prev,g);
			}
			if (projections[i].isForward){
				s2 = source(curr,g);
				t2 = target(curr,g);
			}else{
				s2 = target(curr,g);
				t2 = source(curr,g);
			}
			auto subtraj = std::make_tuple(trajId,std::max(0,breakpts[vote_ct]),
										   std::min(n,breakpts[vote_ct+2]+1));
			
			vote_ct++;
			//		subtrajs.push_back(subtraj);
			bool nearAdj; std::vector<vertex_descriptor> path;
			tie(nearAdj,path)  = graph_->nearAdjacent(t1,s2);
			bool success = jgraph_.voteEdge(g[s1].id, g[t1].id, g[s2].id, g[t2].id, 
											nearAdj, subtraj);
			if (success){

				 std::vector<cv::Point2f> newPts;
				 newPts.push_back(g[t1].p.pos);
				 newPts.push_back(g[s2].p.pos);

				 if (polylineLength(newPts) < max_missing_edge_length_){
					 n_direct_votes_ ++;
					 edge_descriptor newEdge;bool addedNewEdge=false;
					 tie(newEdge, addedNewEdge)=graph_-> addEdge(t1, s2,newPts);
					 if (addedNewEdge){
						 qDebug() <<"added new temporary edge "<< g[t1].id<<"-" <<g[s2].id
								  <<" of length "<< g[newEdge].length;
					 }
				 }else{//same edge
					 qDebug() <<"skip adding new edge of length " << polylineLength(newPts)
							  <<" (max missing edge len="<< max_missing_edge_length_
							  <<")"<< endl;
				 }
			}else if (!nearAdj & path.size()>=2){
				//qDebug() <<"[near adjacent] shortest path len: "<< path.size() ;
				// accumulate votes based on the shortest path
				
				if (jgraph_.voteEdge(g[s1].id, g[t1].id,
								 g[t1].id,g[path[1]].id,  false, subtraj)){
					n_indirect_votes_++;
				}
					
				for (size_t vi=0;vi<path.size()-2;++vi){
					//jgraph_.voteEdge(g[path[vi]].id, 
				
					if (jgraph_.voteEdge(g[path[vi]].id, g[path[vi+1]].id,
									 g[path[vi+1]].id, g[path[vi+2]].id, 
									 false, subtraj) ){
						n_indirect_votes_++;
					}
					
					//auto subtraj = std::make_tuple(trajId, std::max(0,breakpts[vote_ct]),
					//					   std::min(n,breakpts[vote_ct+2]+1));
				}
				
				if (jgraph_.voteEdge(g[path[path.size()-2]].id, g[s2].id,
							 		g[s2].id,g[t2].id,  false, subtraj)) {
					n_indirect_votes_++;
				}
			}else{
				n_filtered++;
			} 
			prev = curr;
		}
	}
	//qDebug() <<"[vote] # attempted votes: " << vote_ct <<" # filtered: " << n_filtered; 
}
void JunctionNetwork::voteAlt(int trajId, std::vector<TrajProjection> projections,
						   std::vector<int> breakpts, int n){

	L1SkeletonGraph::SGraph g = graph_->getGraph();
	int vote_ct =0,n_filtered=0; 
	edge_descriptor curr, prev = projections[0].edge;
	//std::vector<std::vector<int> > allvotes;
	for (size_t i=1;i<projections.size();++i){
		curr = projections[i].edge;

		if (curr !=prev){
			vertex_descriptor s1,t1,s2,t2;
			if (projections[i-1].isForward){
				s1 = source(prev,g);
				t1 = target(prev,g);
			}else{
				s1 = target(prev,g);
				t1 = source(prev,g);
			}
			if (projections[i].isForward){
				s2 = source(curr,g);
				t2 = target(curr,g);
			}else{
				s2 = target(curr,g);
				t2 = source(curr,g);
			}
			auto subtraj = std::make_tuple(trajId,std::max(0,breakpts[vote_ct]),
										   std::min(n,breakpts[vote_ct+2]+1));
			
			vote_ct++;
			//		subtrajs.push_back(subtraj);
			bool nearAdj; std::vector<vertex_descriptor> path;
			tie(nearAdj,path)  = graph_->nearAdjacent(t1,s2);
			bool success = jgraph_.voteEdge(g[s1].id, g[t1].id, g[s2].id, g[t2].id, 
											nearAdj, subtraj);
			if (success){
				n_direct_votes_ ++;
			}else if (!nearAdj & path.size()>=2){
				//qDebug() <<"[near adjacent] shortest path len: "<< path.size() ;
				// accumulate votes based on the shortest path
				
				if (jgraph_.voteEdge(g[s1].id, g[t1].id,
								 g[t1].id,g[path[1]].id,  false, subtraj)){
					n_indirect_votes_++;
				}
					
				for (size_t vi=0;vi<path.size()-2;++vi){
					//jgraph_.voteEdge(g[path[vi]].id, 
				
					if (jgraph_.voteEdge(g[path[vi]].id, g[path[vi+1]].id,
									 g[path[vi+1]].id, g[path[vi+2]].id, 
									 false, subtraj) ){
						n_indirect_votes_++;
					}
					
					//auto subtraj = std::make_tuple(trajId, std::max(0,breakpts[vote_ct]),
					//					   std::min(n,breakpts[vote_ct+2]+1));
				}
				
				if (jgraph_.voteEdge(g[path[path.size()-2]].id, g[s2].id,
							 		g[s2].id,g[t2].id,  false, subtraj)) {
					n_indirect_votes_++;
				}
			}else{
				n_filtered++;
			} 
			prev = curr;
		}
	}
	//	qDebug() <<"[vote] # direct edges: " << nDirectEdges <<" # indirect edges: "
	//<< inDirectEdges << " # of total possible votes: " << vote_ct; 
}
TrajProjection JunctionNetwork::makeTrajProjection(int pid, edge_descriptor edge,
									 cv::Point2f pos, bool isForward){
	TrajProjection tp;
	tp.pid = pid;
	tp.edge = edge;
	tp.pos = pos;
	tp.isForward = isForward;
	return tp;
	
}

void JunctionNetwork::processVotes(int trajId,std::shared_ptr<GpsTraj> traj,
					  const cv::Mat &candidateSegs,const cv::Mat &candidatePos,
  				   const std::vector<TrajProjection > &projections){
	//identify places that is missing edge
	//1) run depth first search to find missing edge
	//2) add edge

	L1SkeletonGraph::SGraph g = graph_->getGraph();
	int successVotes=0, totalVotes=0;
	edge_descriptor  curr, prev = projections[0].edge;
	//int breakPt = projections[0].pid;
	std::pair<vertex_descriptor,vertex_descriptor> lastvotedEdge;

	// for debug-----------------------------
	std::vector<std::vector<int> > allvotes;
	bool hasNewEdge = false;
	//--------------------------------
	//first step: find breakpoints
	std::vector<int> breakpts;
	breakpts.push_back(projections[0].pid);
	
	for (size_t i=1; i< projections.size()-1; ++i){
		curr = projections[i].edge;
		if  (curr != prev){
			breakpts.push_back(projections[i].pid);
			prev = curr;
		}
	}
	breakpts.push_back(projections.back().pid);
			

	std::ostringstream oss1;
	for (size_t p=0;p<breakpts.size();++p)
		oss1<< breakpts[p]<<" ";
	//qDebug()<<"Break points: " << oss1.str().c_str();
	std::ostringstream oss2;
	for (size_t p=0;p<projections.size();++p){
		TrajProjection proj = projections[p];
		oss2 << proj.pid<<">" << g[proj.edge].id<<" ";
	}
	//qDebug()<<"Projections: "<< oss2.str().c_str();	
	//second step: compute votes
	int vote_ct = 0;
	for (size_t i=1; i< projections.size(); ++i){
		curr = projections[i].edge;
		if  (curr == prev){
			//vTraj = traj->point()
			//candidatePos.at<int>(i-1].first
			if (projections[i].isForward){
				lastvotedEdge.first = source(curr,g);
				lastvotedEdge.second = target(curr,g);
			}else{
				lastvotedEdge.first = target(curr,g);
				lastvotedEdge.second = source(curr,g);
			}
			continue; //same edge
		}

		vertex_descriptor s1 = source(prev,g),
			t1 = target(prev,g),
			s2 = source(curr,g),
			t2 = target(curr,g);
		//		auto subtraj = std::make_tuple(trajId, projections.front().pid, projections.back().pid);
		auto subtraj = std::make_tuple(trajId,breakpts[vote_ct], breakpts[vote_ct+2]);
		vote_ct++;
		//		subtrajs.push_back(subtraj);
		
		//breakPt = projections[i].pid;
		bool success=false;
		if (t1==s2 ){  //|| t1==s2){
			//std::get<2>(subtraj) = projections.back().pid;//???
			success = jgraph_.voteEdge(g[s1].id, g[t1].id, g[s2].id, g[t2].id,true,
									   subtraj);
			lastvotedEdge =std::make_pair(s2,t2);

			std::vector<int> arr = {g[s1].id, g[t1].id, g[s2].id, g[t2].id};
			allvotes.push_back(arr);
			
		

		}else if (s1==s2){

			success = jgraph_.voteEdge(g[t1].id, g[s1].id, g[s2].id, g[t2].id,true,
									   subtraj);
			lastvotedEdge =std::make_pair(s2,t2);
			std::vector<int> arr = {g[t1].id, g[s1].id, g[s2].id, g[t2].id};
			allvotes.push_back(arr);
	

		}else if (t1==t2){
			success = jgraph_.voteEdge(g[s1].id, g[t1].id, g[t2].id,g[s2].id,true,
									   subtraj);
			lastvotedEdge =std::make_pair(t2,s2);
			std::vector<int> arr = {g[s1].id, g[t1].id, g[t2].id, g[s2].id};
			allvotes.push_back(arr);
			
		}else if (s1 == t2){
			success = jgraph_.voteEdge(g[t1].id, g[s1].id, g[t2].id,g[s2].id,true,
									   subtraj);
			lastvotedEdge =std::make_pair(t2,s2);
			std::vector<int> arr = {g[t1].id, g[s1].id, g[t2].id, g[s2].id};
			allvotes.push_back(arr);
			
		}else{
			int v0 =g[lastvotedEdge.first].id, v1 = g[lastvotedEdge.second].id;
			int v2,v3;
			float gapLen =inf;
			if (cv::norm(g[lastvotedEdge.second].p.pos - g[s2].p.pos)
				< cv::norm(g[lastvotedEdge.second].p.pos - g[t2].p.pos)){
				v2 = g[s2].id;
				v3 = g[t2].id;
				gapLen =cv::norm(g[lastvotedEdge.second].p.pos - g[s2].p.pos);
				lastvotedEdge = std::make_pair(s2,t2);
			}else{
				v2 = g[t2].id;				
				v3 = g[s2].id;
				gapLen =cv::norm(g[lastvotedEdge.second].p.pos - g[t2].p.pos);		
				lastvotedEdge = std::make_pair(t2,s2);
			}
			//			if (cv::norm( g[s]
			if (gapLen < min_edge_length_){
				/*
					qDebug() <<"[Voting] Add new edge"
					 << g[prev].id
					 <<"(" << v0 <<","<< v1 <<") "
					 << "->"<< g[curr].id
					 <<"(" << v2 <<","<< v3 <<") ";*/
				hasNewEdge = true;
				std::vector<int> arr = {v0,v1,v2, v3};
				allvotes.push_back(arr);
				try{
					jgraph_.voteEdge(v0 ,v1 ,v2,v3, true,subtraj);
				}catch(...){
					//qDebug() <<"Error: exception caught on voteEdge !!!!!!!!!";
					//	debugEdgeVotes(trajId,allvotes);
				}

			}/*else{
				qDebug() <<"[Voting] Skip adding new edge -- length exceeds ";
				}*/
		}		

		//debugEdgeVotes(trajId,allvotes);
		
		if (success){
			successVotes++;
		}
		/*
		if (hasNewEdge)
			debugEdgeVotes(trajId,allvotes);
		*/
		totalVotes++;
	    prev = curr;
	}
}
void JunctionNetwork::debugEdgeVotes(int id, std::vector<std::vector<int> > votes){

	for (size_t i=0;i<votes.size();++i){
		qDebug() << "Traj " << id << " vote: (" << votes[i][0] <<"->"
				 << votes[i][1] <<") => ("
				 << votes[i][2] << "->" << votes[i][3]<<")";
	}
}

point JunctionNetwork::toBoostPoint(const TrajPoint* tp,
                   osg::Vec3 min_corner,float resolution){
    float x = tp->x();
    float y = tp->y();

    x = (x - min_corner.x()) / (float)resolution;
    y = height_ - (y - min_corner.y()) / (float)resolution;

    return point(x,y);

}

cv::Point2f JunctionNetwork::toCvPoint2f(const TrajPoint* tp, osg::Vec3 min_corner,
                                         float resolution){
    float x = tp->x();
    float y = tp->y();

    x = (x - min_corner.x()) / (float)resolution;
    y = height_ - (y - min_corner.y()) / (float)resolution;

    return cv::Point2f(x,y);
}

cv::Point2f JunctionNetwork::toCvPoint2f(const TrajPoint* tp,
										 cv::Point2f min_corner,
										 float height,
                                         float resolution){
    float x = tp->x();
    float y = tp->y();

    x = (x - min_corner.x) / resolution;
    y = height - (y - min_corner.y) / resolution;

    return cv::Point2f(x,y);
}



rtree JunctionNetwork::makeRTreeFromEdges(){
    std::vector<value> segments;
    //iterate through skeleton graph edges,
    std::vector<std::vector<cv::Point2f> > edges= graph_->getEdgeSegments();

    int i =0 ;
    for (auto it=edges.begin();it!=edges.end();++it){
        for (size_t j=0; j<it->size()-1; ++j){
            cv::Point2f &v1 = (*it)[j];
            cv::Point2f &v2 =  (*it)[j+1];
			if (v1==v2) continue;//avoid segment of length 0
            point p1(v1.x,v1.y), p2(v2.x,v2.y);
            segment seg( p1, p2 );
            segments.push_back(std::make_tuple(seg,i,j));
        }
        i++;
    }
    rtree rt(segments.begin(), segments.end());
    return rt;
}
bool JunctionNetwork::projectPointOnSegment(cv::Point2f in,
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
    normalize(e,en);
    normalize(v,vn);
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
bool JunctionNetwork::projectPointOnLine(cv::Point2f in,
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
    normalize(e,en);
    normalize(v,vn);
    float projLen = v.x*en.x + v.y*en.y;
    float eLen = cv::norm(e);
    bool onSegment = (projLen <0 || projLen>eLen)? false:true;
    out = p1+en*projLen;


    return onSegment;
}
const std::string  JunctionNetwork::valueToStr(const value &v) &{
    std::ostringstream oss;
    oss << std::get<1>(v) <<"(" << std::get<2>(v)<<")";
    return oss.str();
}
std::vector<cv::Point2f>  JunctionNetwork::getIntersectionPoints(){
    std::vector<cv::Point2f> pts;
    for (auto it = intersections_.begin();it!=intersections_.end();++it){
        pts.push_back( (*it)->loc);
    }
    return pts;
}

void JunctionNetwork::clusterEndpoints(float r){
	//1. create rtree of endpoints
	qDebug() << "====== Cluster endpoints of radius "  << r <<" ======";
	L1SkeletonGraph::SGraph g = graph_->getGraph();
	L1SkeletonGraph::vertex_iterator vit,vit_end;
	std::vector<pvalue> pointList;


	int i=0;
	for (tie(vit,vit_end) = vertices(g);vit!= vit_end;++vit ){
		cv::Point2f pt = g[*vit].p.pos;
		pointList.push_back(std::make_pair(point(pt.x,pt.y),i++));
		
	}
	bool visited[pointList.size()];
	for (size_t pi=0;pi<pointList.size();++pi){
		visited[pi]=false;
	}
	
	
	prtree rt(pointList.begin(), pointList.end());
	for (size_t k=0; k< pointList.size();++ k){
		if (visited[k]) continue;
		prtree::const_query_iterator nit;
		std::unordered_set<size_t> pset={k};

		vertex_descriptor vk = vertex(pointList[k].second,g);

		
		for (nit = rt.qbegin(bgi::nearest(pointList[k].first, 10));nit!=rt.qend();
			 ++nit){
			vertex_descriptor vn = vertex(nit->second,g);
			if (vn == vk) continue;
			double d = bg::distance(pointList[k].first, nit->first);
			if (d>r){
				break;
			}
			visited[nit->second] = true;
			// test if vertex (index=pointList[k].second) is reachable within
			// 2 hops
			std::vector<cv::Point2f> pts(2);
			pts[0] = g[vk].p.pos; pts[1] = g[vn].p.pos;
			if ( ! graph_->bfs(vk,vn,2)){
				//add edge between vertex pointList[k] and nit
				edge_descriptor e ;bool success;
				tie(e,success) = graph_->addEdge(vk,vn,pts );
				/*
				if (success)
					qDebug() <<"Endpoint clustering: add edge between "
							 << g[vk].id <<" and "
							 << g[vn].id << " of length " <<  g[e].length;
				*/
			}
		}
		//std::vector<int> pids = g->
	}
	graph_->updateEdgeIndex();
	graph_->computeEdgeWeight();
	graph_->assignVertexType();
	// for each unvisited endpoint u, find its r-neihbors
	// add an edge from u to every other neighbor that is
	// eighter not connected, or too far away (use breadth first search of radius2)
	//!!!tobeimplemented!!
}
void JunctionNetwork::processIntersections(){
    if (intersections_.size()>0){
    L1SkeletonGraph::SGraph g = graph_->getGraph();

    // store edge descriptors in a map
    std::map<int,edge_descriptor> edgeMap ;
    for (auto it=extended_edges_.begin(); it!=extended_edges_.end();++it){
         edge_descriptor e;bool exists;
         tie(e,exists) = graph_->getEdgeById((*it)->id);
         if (exists && (*it)->id < extended_edges_.size()){
            edgeMap[(*it)->id] =e;
         }
    }
    // process intersection
    for (auto it=intersections_.begin();it!=intersections_.end();++it){
       value qSeg = (*it)->qSeg, nSeg = (*it)->nSeg;
       shared_ptr<ExtendedEdge> qe = extended_edges_[std::get<1>(qSeg)];
       shared_ptr<ExtendedEdge > ne = extended_edges_[std::get<1>(nSeg)];

       /*qDebug() << "[ProcessIntersection] qSeg = " << valueToStr(qSeg).c_str()
                <<" nSeg = " << valueToStr(nSeg).c_str()
				<<" loc: (" << (*it)->loc.x <<"," << (*it)->loc.y <<")";*/
	   
        bool qIsHead = ! (std::get<2>(qSeg)==0);
        edge_descriptor qd =edgeMap[qe->id];
		if (edgeMap.find(ne->id)==edgeMap.end() ){
			/*qDebug()<<"Error: cannot find neighbor edge " << ne->id
			  << "of the intersection!";*/
			continue;
		}
        edge_descriptor nd =edgeMap[ne->id];
        vertex_descriptor  q1 = (qIsHead)? target(qd,g) :source(qd,g),
                           q2 = (qIsHead)? source(qd,g) :target(qd,g);
        vertex_descriptor n1 =source(nd, g), n2 = target(nd,g);
       if ( std::get<2>(nSeg)==0 ){ //neighbor segment is a tail segment
		   //qDebug() <<"  neighbor segment edge_id=" << g[nd].id  << "iis a tail segment";
			
            /* edge labels:
             *  n1---x----x--->n2
             *    q1
             *    |
             *    x
             *    |
             *    q2
             */

            graph_->moveVertex( n1, (*it)->loc);
          //  graph_->moveVertex( q1, (*it)->loc);
             std::vector<cv::Point2f> pts =g[qd].internalPts;
            if  (qIsHead){
				pts.push_back( (*it)->loc );
            }else{
                std::reverse(pts.begin(),pts.end());
                pts.push_back( (*it)->loc );
               // pts.insert(pts.begin(), (*it)->loc);
            }
            graph_->addEdge(q2, n1, pts);
            graph_->removeEdge(q1, q2);

       }else if (std::get<2>(nSeg) == ne->pts.size()-2){
           //neighbor segment is a head segment
		   //qDebug() <<"  neighbor segment edge_id=" << g[nd].id  << "is a head segment";
            /* edge labels:
             *  n1---x----x--->n2
             *               q1
             *                |
             *                x
             *                |
             *               q2
             */

            graph_->moveVertex( n2, (*it)->loc);
            std::vector<cv::Point2f> pts = g[qd].internalPts;
            if  (!qIsHead){
				std::reverse(pts.begin(),pts.end());
            }
            pts.push_back( (*it)->loc);
            //graph_->moveVertex( q1, (*it)->loc);
            graph_->addEdge(q2,n2, pts);
            graph_->removeEdge(q1, q2);
       }else{
           //qDebug() << "  intersection is on internal segment edge_id=" << g[nd].id  ;
           /* edge labels:
            *  n1---x-----x--->n2
            *          q1
            *          |
            *          x
            *          |
            *          q2
            */
		   if (g[nd].id > extended_edges_.size()){
			   //qDebug() <<"Error: vertex id " << g[nd].id <<" exceeded length "
			   //						 << extended_edges_.size();
				continue;
			}
            graph_->moveVertex(q1,(*it)->loc);
            int l = std::get<2>(nSeg);cv::Point2f s1,s2;
            cv::Point2f qv , nv = ne->pts[l+1] - ne->pts[l];
            if (qIsHead){
                qv =boostToCVPoint(std::get<0>(qSeg).second) - boostToCVPoint(std::get<0>(qSeg).first);
            } else{
                qv =boostToCVPoint(std::get<0>(qSeg).first) - boostToCVPoint(std::get<0>(qSeg).second);
            }

            float cosa =L1SkeletonGraph::cosAngleBetween(qv,nv);
            bool sourceIsTerminal = g[n1].type==SkeletonGraphVertexType::TERMINAL;
            bool targetIsTerminal = g[n2].type==SkeletonGraphVertexType::TERMINAL;
            if (cosa < -0.85 && targetIsTerminal){ // keep only s-->intersection
                graph_->splitEdgeAtIntersection(nd,q1, std::get<2>(nSeg),
												(*it)->loc, true, false );
            }else if (cosa >0.85 && sourceIsTerminal){//keep only intersection ---> target
                 graph_->splitEdgeAtIntersection(nd,q1, std::get<2>(nSeg),
												 (*it)->loc, false, true );
            }else{
                graph_->splitEdgeAtIntersection(nd,q1, std::get<2>(nSeg),
												(*it)->loc, true, true ); //need to check direction of split!
            }
       }
    }
    graph_->updateVertexPos();
    graph_->removeDanglingEdges();
    graph_->assignVertexType();

}else{
        graph_->updateVertexPos();
        graph_->assignVertexType();
    }
}

vertex_descriptor JunctionNetwork::qEndpoint(const EdgeIntersection &intersect)
{
    L1SkeletonGraph::SGraph g = graph_->getGraph();
    edge_descriptor qEdge; bool exists;
    tie(qEdge, exists) = graph_->getEdgeById(std::get<1>(intersect.qSeg));
    if (!exists){
        return vertex_descriptor();
    }else{
        return (intersect.qType == HEAD)? target(qEdge, g):source(qEdge,g);
    }
}
float JunctionNetwork::computeEmissionCost(float dist, float mean, float var){
    if (dist <0){
        return 1e6;
    }
    return -(1.f/sqrt(2*M_PI*var) * exp( -pow(mean-dist,2)/(2*var) ));
}
float JunctionNetwork::computeTransitionCost(float pathlen, float d2,float lambda){
    if (pathlen==0) return 0.f;
    return pathlen ; //std::fabs(pathlen-d2)/float(pathlen);
 }



//void JunctionNetwork::processTrajectory()
std::vector<TrajProjection> JunctionNetwork::projectGPSToGraph(
				   std::shared_ptr<GpsTraj> traj,
                   const cv::Mat &candidateSegs,  const cv::Mat &candidatePos,
                   const cv::Mat &distances,  const osg::Vec3 &min_corner){

    float lambda = Parameters::getInstance()->skeleton_junction_smoothness;
    L1SkeletonGraph::SGraph g=graph_->getGraph();
    size_t nPoints= candidateSegs.rows;
    size_t maxNumCandidates =candidateSegs.cols;

    // compute projection of traj on candidates and variance of projection distance
    //sequence of edgeids between candidates (i,j) and (i+1,k)
    std::map<std::tuple<size_t,size_t,size_t> ,std::list<edge_descriptor> > shortestpaths;

    //qDebug()<<"==Variance on candidate distances =======";

    std::vector<std::vector<cv::Point2f> > candidateProj(nPoints);
    std::vector<float> variance(nPoints);
    for (size_t i=0; i<nPoints;++i){
        //candidateProj[i].resize(maxNumCandidates,cv::Point2f(0,0));
		int nCand = 0;
        for (size_t j=0; j< maxNumCandidates;++j){
            int edgeId = candidateSegs.at<int>(i,j);
            int segId = candidatePos.at<int>(i,j);
            if (edgeId == -1 || segId ==-1)
                break;
            edge_descriptor e = graph_->getEdgeById(edgeId).first;
            std::vector<cv::Point2f> shape = g[e].internalPts;
            cv::Point2f p1 = shape[segId], p2 = shape[segId+1];
            cv::Point2f p = toCvPoint2f(&(traj->point(i)), min_corner,true);
			cv::Point2f proj;
            projectPointOnSegment(p,p1,p2 ,proj);
			candidateProj[i].push_back(proj);
			nCand++;
        }
        variance[i] = rowVariance<float>(distances.colRange(0,nCand),i,0.f);
    }

    // compute emissionCost  and transitionCost
    cv::Mat emissionCost(nPoints, maxNumCandidates, CV_32F, cv::Scalar(1e6));
    std::vector<cv::Mat> transitionCost(nPoints-1);
	
    for (size_t i=0; i < nPoints;++i){
        cv::Point2f p = toCvPoint2f(&(traj->point(i)), min_corner,true);
        for (size_t j=0; j <candidateProj[i].size();++j){			
            float d = distances.at<float>(i,j);
            emissionCost.at<float>(i,j) = computeEmissionCost(d,0.f,variance[i]);
		}


        if  ( i == nPoints-1){
			continue;
		}
		transitionCost[i]=cv::Mat(candidateProj[i].size(),candidateProj[i+1].size(),
							   CV_32F,cv::Scalar(1e6));
		
        cv::Point2f pNext  = toCvPoint2f(&(traj->point(i+1)), min_corner,true);
		for (size_t j=0; j<candidateProj[i].size();++j){

			int edgeId = candidateSegs.at<int>(i,j);
            int segId = candidatePos.at<int>(i,j);

            edge_descriptor e = graph_->getEdgeById(edgeId).first;
            std::vector<cv::Point2f> shapej = g[e].internalPts;

            size_t n =num_vertices(g);
            vertex_descriptor from_vertex1 =source(e, g), from_vertex2 = target(e,g);
            std::vector<float> distances1(n), distances2(n);
            std::vector<vertex_descriptor> parents1(n), parents2(n);
            graph_->dijkstras(from_vertex1 , distances1, parents1);
            graph_->dijkstras(from_vertex2 , distances2, parents2);

			for (size_t k=0; k<candidateProj[i+1].size();++k){
				int segId_next = candidatePos.at<int>(i+1,k);
                int edgeId_next = candidateSegs.at<int>(i+1,k);

                float d2 = cv::norm(candidateProj[i][j] - candidateProj[i+1][k]);
                std::list<edge_descriptor> shortestpath;
                float d = computePathLen(candidateProj[i][j],
								   candidateProj[i+1][k],
								   distances1, distances2,
								   parents1, parents2,
                                   edgeId,edgeId_next,
								   segId,segId_next,shortestpath);
				if (shortestpath.size()>0){
					auto pathKey = std::make_tuple(i,j,k);
					shortestpaths[pathKey] = shortestpath;
				}

				/**
				 * compute Transition cost
				 */
				float weight = (edgeId==edgeId_next)?1:2;
                transitionCost[i].at<float>(j,k)
					= pow(computeTransitionCost(d,d2,lambda),weight);
				
            }
        }
    }
	//------start debug------------------------------------------
    std::ostringstream oss;
    oss << emissionCost << endl;
    for (size_t tt = 0; tt< nPoints-1;++tt ){
        oss <<"\n== transitionCost["<<  tt  << "]===============\n";
        oss << transitionCost[tt];
    }
    qDebug() << oss.str().c_str() ;
	//------end debug------------------------------------------
	
    //use vertibi algorithm to compute optimal projection
     Viterbi solver(emissionCost, transitionCost); //<--do backtracing, shortest path
	 bool success = solver.forward();
	 if (! success){
		 qDebug() << "Unable to solve for optimal projection on trajectory "
				  <<"of length " << nPoints;
		 return  std::vector<TrajProjection>();
	 }
	 /* shortestpath
	 qDebug()<<"==== Shortest Paths ======";
	 for (auto p=shortestpaths.begin();p!=shortestpaths.end();++p){
		 //std::ostringstream oss;
		 std::list<edge_descriptor> paths;
		 for (std::list<edge_descriptor>::iterator l=p->second.begin();
			  l!= p->second.end();++l){
			 //oss << g[*l].id << " "; 
		 }
		 qDebug()<<"path "<< std::get<0>(p->first)<<"," << std::get<1>(p->first)
				 << "," << std::get<2>(p->first)
				 <<" : length=" << p->second.size() <<" "
				 << oss.str().c_str();
	 }*/
	 
     // construct full path (sequence of edge descriptors)
	 //     std::vector<std::pair<int,edge_descriptor> > fullpath;
	 std::vector<TrajProjection> fullpath;
     std::vector<int> opt_selection = solver.getOptimalPath();
     for (size_t i=0;i<opt_selection.size()-1;++i){
         int j = opt_selection[i];
         int k = opt_selection[i+1];

         edge_descriptor selectedEdge;bool exists;
		 if  (j<0 || k<0 || j>=candidateProj[i].size()
			  || k>= candidateProj[i+1].size()){
			 qDebug() << "Error: Optimal selection " << j << " out of bound";
			 continue;
		 }
         tie(selectedEdge,exists)= graph_->getEdgeById(candidateSegs.at<int>(i,j));
         if (!exists)
             continue;

        //Need to check if i and i+1 line up
        auto key = std::make_tuple(i,j,k);
        auto iter = shortestpaths.find(key);

		cv::Point2f projPos0 = candidateProj[i][j],
			projPos1 = candidateProj[i+1][k];
		auto pts = g[selectedEdge].internalPts;
		int pp= candidatePos.at<int>(i,j);
		cv::Point2f vec ;
		if (pp==0){
			vec = pts[pp+1]-pts[0];
		}else if (pp==pts.size()-1){
			vec = pts[pp]-pts[pp-1];
		}else{
			vec = pts[pp+1]-pts[pp-1];
		}
		cv::Point2f vt = LShapeIterator::UnitDirection(traj->point(i).head());
		bool forward =  vt.dot(vec)>0;
        fullpath.push_back(makeTrajProjection(i, selectedEdge,
										projPos0,forward));
        if (iter != shortestpaths.end()){			
            std::list<edge_descriptor> subpath =  iter->second;
			for (auto s=subpath.begin();s!=subpath.end();++s){
				fullpath.push_back(makeTrajProjection(i, *s, projPos0, true));  
			}
        }
        if (i==opt_selection.size()-2){
            tie(selectedEdge,exists) = graph_->getEdgeById(candidateSegs.at<int>(i+1,k));
            if (exists){
                fullpath.push_back(makeTrajProjection(i+1,selectedEdge, projPos1,true)); 
            }
        }
     }
    return fullpath;
}

void JunctionNetwork::debugDijkstras(const std::vector<float> &distances,
							const std::vector<vertex_descriptor> &parents){
	std::ostringstream oss;
	for (size_t i=0; i<distances.size();++i){
		oss << i <<" " << distances[i] << " " << parents[i] << "\n";
	}
	qDebug() << oss.str().c_str();
}

float JunctionNetwork::computePathLen(const cv::Point2f  &proj_j,
									  const cv::Point2f &proj_k,
                    const std::vector<float> &distances1,
									  const std::vector<float> &distances2,
                    const std::vector<vertex_descriptor> &parents1,
                    const std::vector<vertex_descriptor> &parents2,
                    int edgeId, int edgeId_next, int segId, int segId_next,
                     std::list<edge_descriptor> & shortestpath ){
	L1SkeletonGraph::SGraph g = graph_->getGraph();
    edge_descriptor e = graph_->getEdgeById(edgeId).first;
    std::vector<cv::Point2f> shapej = g[e].internalPts;

     edge_descriptor  eNext = graph_->getEdgeById(edgeId_next).first;
    std::vector<cv::Point2f> shapek =  g[eNext].internalPts;

    float d =0.0f;
    // case 1: same segment same edge
    bool sameSeg = segId == segId_next;
    bool sameEdge = edgeId == edgeId_next;
    //minimum path length distance between two projections
    float d2 = cv::norm(proj_j -  proj_k);
     if (sameSeg && sameEdge){
        // compute projection distance
        d = d2;

    }else if (sameEdge){
        d = polylineLength(shapej,std::min(segId,segId_next)+1,std::max(segId,segId_next));

    }else {
		 d= 0;
		 int sj = segId;
		 int sk = segId_next;
		 float d11,d21,d22,d12;
		 d11 =  polylineLength(shapej, 0,sj) + polylineLength(shapek,0,sk)
			 + cv::norm(proj_j-shapej[sj]) + cv::norm(proj_k-shapek[sk]);
		 d12 =  polylineLength(shapej,0,sj)+cv::norm(proj_j-shapej[sj])
			 + polylineLength(shapek, sk+1,shapek.size()-1)
			 + cv::norm(proj_k-shapek[sk+1]);
		 d21 = polylineLength(shapej,sj+1,shapej.size()-1)
			 + cv::norm(proj_j - shapej[sj+1])
			 + polylineLength(shapek,0,sk)+  cv::norm(proj_k-shapek[sk]);
		 d22 = polylineLength(shapej,sj+1,shapej.size()-1)
			 + cv::norm(proj_j - shapej[sj+1])
			 + polylineLength(shapek, sk+1,shapek.size()-1)
			 + cv::norm(proj_k- shapek[sk+1]);
			
		 // bool eForward=true, eNextForward=true;
		 if (graph_->incidentEdges(e,eNext)){
			 if (source(e,g)  == source(eNext,g)){
				 //shapej[0] == shapek[0]){
				 //qDebug() << "case 3.b [s->s]";
				 d = d11;
			 }else if (source(e,g)  == target(eNext,g)){
				 //shapej[0] == shapek[shapek.size()-1]){
				 //qDebug() << "case 3.d [s->t]";
				 d = d12;
			 }else if (target(e,g) == source(eNext,g)){
				 //shapej[shapej.size()-1] == shapek[0]){ //
				 //qDebug() << "case 3.a [t->s]";
				 d = d21; 
			 }else if (target(e,g) == target(eNext,g)){
				 //shapej[shapej.size()-1] ==shapek[shapek.size()-1] ){//
				 // qDebug() << "case 3.c [t->t]";
				 d = d22; 
			 }

			 //float dtemp = d;
			 //d = std::max(d2,d);
		 }else{
			 float maxGapLen = 1.5*min_edge_length_;
			 float e11 = cv::norm(shapej[0] - shapek[0])*1.5,
				 e12 = cv::norm(shapej[0] - shapek[shapek.size()-1])*1.5,
				 e21 = cv::norm(shapej[shapej.size()-1] - shapek[0])*1.5,
				 e22  = cv::norm(shapej[shapej.size()-1]
								 - shapek[shapek.size()-1])*1.5;
			 //e11 = (e11 >maxGapLen)? inf :e11;
			 //e12 = (e12 > maxGapLen)?inf:e12;
			 // e21 = (e21 > maxGapLen)?inf:e21;
			 // e22 = (e22 > maxGapLen ) ?inf:e22;
			 // compute four possible orientations, pick the shortest distance
			 vertex_descriptor from_vertex1 = source(e,g),
				 from_vertex2 = target(e,g);
			 vertex_descriptor to_vertex1 =source(eNext, g),
				 to_vertex2 = target(eNext,g);

			 d11 += distances1[to_vertex1];//e11) ;
			 d12 += distances1[to_vertex2];//e12);
			 d21 += distances2[to_vertex1];//,e21);
			 d22 += distances2[to_vertex2];//e22);
			 float elem[] = {d11,d12,d21,d22};
			 auto iter =std::min_element(elem,elem+4);
			 float maxDist = 1000;
			 if (*iter>= maxDist ){// maxGapLen){ //implies no path!
				 // qDebug() <<"Can not find path from edge " << edgeId <<
				 //" to edge " << edgeId_next;
				 d =  maxDist;
			 }else{
				 size_t minIndex =  iter-elem;
				 if (minIndex ==0){
					 //qDebug() <<"\n[case s->s]  shortest path len=" << distances1[to_vertex1] ;
					 shortestpath =  backtrace(parents1,from_vertex1,to_vertex1);
					 /*qDebug() << "shortest path len="
					   <<  distances1[to_vertex1] <<" N=" << shortestpath.size() ;*/

				 }else if (minIndex ==1){
					 //qDebug() <<"\n[case s->t]  shortest path len=" 

					 shortestpath = backtrace(parents1, from_vertex1, to_vertex2);
					 /*			   qDebug() << "shortest path len="
								   << distances1[to_vertex2] <<" N=" << shortestpath.size();*/
			   
				 }else if (minIndex ==2){
					 //qDebug() <<"\n[case t->s]  shortest path len="
					 shortestpath = backtrace(parents2, from_vertex2,to_vertex1);
					 /*			   qDebug() << "shortest path len=" << distances2[to_vertex1]
								   << " N=" << shortestpath.size();*/

				 }else{
					 //qDebug() <<"\n[case t->t]  shortest path len=" << distances2[to_vertex2] ;
					 shortestpath = backtrace(parents2, from_vertex2,to_vertex2);
					 /*qDebug() << "shortest path len=" << distances2[to_vertex2]
					   <<" N=" << shortestpath.size();*/


				 }
				 d = *iter;
			 }
		  /*
		  float upperDistBound = d2*4.f;//cv::norm(proj_j - proj_k)*2;//2* projectionDist

		  //if (d> upperDistBound){
			  
			  vertex_descriptor vf= (cv::norm(proj_k-g[from_vertex1].p.pos)
									 <cv::norm(proj_k-g[from_vertex2].p.pos))?
				  from_vertex1:from_vertex2;
			  vertex_descriptor vt= (cv::norm(proj_j-g[to_vertex1].p.pos)
									 <cv::norm(proj_j-g[to_vertex2].p.pos))?
				  to_vertex1:to_vertex2;
			  
			  //set shortestpath to to_vertex -> from_vertex
			
			  d = d2*1.5;
			  qDebug() << "Set pseudo distance between " << g[vf].id <<" -> "
					   <<g[vt].id <<" to " << d ;			  
					   }*/
			  //1e2;//@@Need to use bounding box size...
		  //default parent = self, default dist =inf
		  //d = std::max(upperDistBound,d);
		
		  
          /* qDebug() <<"Candidate transition from edge" << edgeId<<"."<< segId
                   << " to edge " << edgeId_next <<"."
                   <<    segId_next
                    << "** Non-incident edges ** pathlen="<< d <<" d2=" << d2
					<< " shortestpath length " << shortestpath.size();*/
		  //}
		 }
	 }
	 return  std::max(d,d2);
}


std::list<edge_descriptor> JunctionNetwork::backtrace(
						  std::vector<vertex_descriptor> parents,
						  vertex_descriptor v_from, vertex_descriptor v_to){
	if (parents[v_to] == -1 || parents[v_to]==v_to){
		auto g = graph_->getGraph();
		//qDebug() <<"Unable to backtrace from target " << g[v_to].id
		//		 << " to source " << g[v_from].id;
		return std::list<edge_descriptor>();
	}
    L1SkeletonGraph::SGraph g = graph_->getGraph();
    vertex_descriptor prev = v_to;
    std::list<edge_descriptor> path;

	//edge_descriptor lastEdge = edge(v_to
	//path.push_front(lastEdeg);

    while (prev!=v_from){

        vertex_descriptor temp = prev;
        prev = parents[prev];
        edge_descriptor e; bool exists;
        tie(e, exists) = edge(prev,temp,g);
        if (exists){
            path.push_front(e);
        }else{
			qDebug() <<"cannot find backtrace edge " << g[prev].id <<"--" << g[temp].id;
		}
    }
    return path;
}


int JunctionNetwork::getMaxClusterSize(){
    int max=0;
    for (auto it= clusters_.begin();it!= clusters_.end();++it){
        max=std::max(max,(int)it->second.size());
    }
	return max;	
}

QStringList JunctionNetwork::clusterToQStringList(){
    QStringList slist;
        int i=0;
    for (auto it= clusters_.begin();it!= clusters_.end();++it){

		slist << QString("Cluster %1: %2--%3 (%4)").arg(i).arg(it->first[0]).arg(it->first[1]).arg(it->second.size());
		i++;

    }
    return slist;	
}
