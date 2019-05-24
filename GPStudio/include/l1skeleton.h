#pragma once
#ifndef L1_SKELETON_H
#define L1_SKELETON_H

#include <opencv2/core/core.hpp>
#include <vector>
#include <opencv2/flann/flann.hpp>
#include "common.h"
//===================== auxillary data structure====================
enum SkeletonPointType{
	NONE=0,
	BRANCH=1,
	BRIDGE=2,
	NOISE=3
};
enum BranchEndpointType{
	FIRST,
	LAST
};
/**
 * @brief The SkeletonPoint struct : a skeleton point
 */
struct SkeletonPoint{
    int id;             // unique id
    cv::Point2f pos;    // position
    SkeletonPointType type;	 //type flag

    SkeletonPoint():id(-1),pos(cv::Point2f(-1,-1)),
		type(SkeletonPointType::NONE) {}

    SkeletonPoint(int index, cv::Point2f position, SkeletonPointType type):
	    id(index), pos(position), type(type){}

	virtual bool operator == (const SkeletonPoint &rhs)const{
		return id == (rhs.id);
	}
};
/**
 * @brief The Candidate struct : a skeleton point that is
 * a branch candidate
 */
struct Candidate: public SkeletonPoint{

	cv::Point2f dir; // PCA direction
    float sigma;     // linearity (sigma) score
	float proj_dist; // projection distance 
	bool visited;   // whether it has been visted
	bool operator == (const Candidate &rhs)const{
		return id == (rhs.id);
	}
    Candidate():SkeletonPoint(),dir(cv::Point2f(-1,-1)), sigma(0),
		proj_dist(-1),visited(false){id=-1;	}

    Candidate(SkeletonPoint pt, cv::Point2f dir, float sigma, 
			  float proj_dist, bool visited):SkeletonPoint(pt), 
		dir(dir),sigma(sigma),proj_dist(proj_dist),	visited(visited){}

    // default comparator based on projection distance
	bool operator < (const Candidate &rhs)const{
		return proj_dist < rhs.proj_dist;
	}
	SkeletonPoint toSkeletonPoint(){
		return SkeletonPoint(id,pos,type);
	}
};
/**
 * @brief Candidate struct comparator based on sigma
 */
struct compareSigma {
  bool operator() (const Candidate &lhs, const Candidate &rhs){
	  if (lhs.visited  == rhs.visited)
		  return lhs.sigma < rhs.sigma;
	  else if (lhs.visited)
		  return true;//visited candidate is always less than unvisited ones
	  else
		  return false;
  }
};

/**
 * @brief Candidate struct comparator based on id
 */
struct compareId {
  bool operator() (const Candidate &lhs, const Candidate &rhs){
	return lhs.id < rhs.id;
  }
};
/**
 * @brief The unvisited struct predicate
 */
struct unvisited{
	bool operator()(const Candidate& c) {
		return !(c.visited);
	}
};


//=================== L1Skeleton base class =====================
/**
 * @brief The L1Skeleton_ class: extracting a l1 medial skeleton
 * from voting image
 */
class L1Skeleton_{
public:
  /**
   * @brief BranchGraph datastructure to store connectivity of branches map<bridgePt,branches>
   * bridgePt: bridge point's unique skeleton id
   * branches: ids of branches adjacent to bridge Pt.
   */
  typedef std::map<int, std::set<int> > BranchGraph;
  /**
   * @brief The SearchResult struct
   */
  struct SearchResult{
     cv::Mat indices;
     cv::Mat distances;
  };
  L1Skeleton_(const cv::Mat &voting_image,
			  const std::vector<cv::Mat> &heading_images);
  virtual ~L1Skeleton_();
  //======== accessors and untility =============================
  std::vector<SkeletonPoint> getSamples(){return I_;}
  std::vector<SkeletonPoint> getSkeleton(){return X_;}
  int size(){return X_.size();} 
  float getH0() const {return h0_;}
  float getSigmaThresh() const {return sigma_thresh_;}
  virtual int nBranches(){return branches_.size();}
  virtual std::vector<Candidate> getBranchById(int i){return  branches_[i];}
  virtual const std::vector<std::vector<Candidate> > & getBranches(){
	  return branches_;
  }
  virtual std::map<int, std::set<int> > getBranchGraph(){return  branchGraph_;}
  virtual float getAdaptiveH(float x, float y);
  virtual float theta(float d, float h);
  virtual float thetaMVN(cv::Point2f x,float h, cv::Mat C){return 0.0;}

  virtual float thetaMVNRotated(cv::Point2f x,float h,cv::Mat C);
  size_t getNumHeadingImages(){return heading_images_.size();}

  virtual cv::Mat getSkeletonCovariance(int idx){return cv::Mat(2,2,CV_32F,cv::Scalar(0));}
  //======= compute skeleton ====================================
  virtual void initialize(float vmin=-1);
  virtual std::vector<IndexVal<float> > findNearestKSkeletonPoints(float px,
                                                                   float py,int k);
  virtual float getSkeletonLinearity(int i);
  virtual void reinitialize(bool incrementH0=false);
  virtual void identifyBranchPoints();
  virtual bool regularize(int max_iter=0);


  virtual std::vector<SkeletonPoint> getBridgePoints();
  virtual Candidate getCandidateInBranch(int branchId, int skeletonId);


  virtual void computeBranchCandidate(std::vector<Candidate> &candidates,
									  bool excludeBranchPoints);
  virtual void visualizeSigma();
  virtual void visualizeDensity(){}
  /**
   * @brief search candidate by skeleton id
   **/
  static std::vector<Candidate>::iterator
      searchCandidates(std::vector<Candidate> &candidates,int i);
  static float angleBetween(cv::Point2f v1, cv::Point2f v2);
  /**
   * compute the average sigma of a skeleton point's k nearest 
   * neighbors
   */
  virtual double computeSkeletonAvgSigma(int skeletonId,int k);

  //=================public methods to reimplement ===============
  virtual SearchResult findSkeletonNeighbors(cv::Point2f p){return SearchResult();}
  virtual bool iterate_step(){return false;}
  virtual std::pair< std::vector<Candidate>, std::vector<Candidate> > 
	  getBranchFromSeed(int){
	  return  std::pair< std::vector<Candidate>, std::vector<Candidate> >();
  }
  virtual std::vector<SkeletonPoint> 
	  mergeSkeletonPoints(std::vector<SkeletonPoint> skeletonPts){
	  return std::vector<SkeletonPoint>();
  }
  virtual cv::Point2f getPrincipalDirection(int i){
	  return cv::Point2f(0,0);
  }
  virtual float getDensity(int x, int y);
  virtual float getDensity(int x, int y,float r);
  virtual cv::Mat getDensityMatrix(){return density_;}
  virtual cv::Mat getHeadingDistribution(int x, int y, bool recursive = true){return cv::Mat();}
  //========= debug/export/visualization ===========
  cv::Mat getSigma(){return sigma_;}
  int getIter(){return iter_;}
protected:
  //========= protected members ====================
  cv::Mat voting_image_;
  std::vector<cv::Mat> heading_images_;
  cv::Mat features_;
  cv::Mat skeleton_features_;
  cv::Mat density_;
  float h0_;
  float sigma_thresh_;
  float skeleton_thresh_;
  cv::Mat alpha_;
  cv::Mat beta_;
  cv::Mat sigma_;
  cv::Mat pca_directions_;

  std::vector<SkeletonPoint> I_; /////
  std::vector<SkeletonPoint> X_;
  int inputSize_;
  int iter_; // current iteration
  cv::Mat Q_;
  cv::flann::GenericIndex<cvflann::L2<float> > *inputKNNIndex_;
  std::vector<std::vector<Candidate> > branches_;
  BranchGraph branchGraph_;

  //========= protected member functions ==================
  void sampleUniform();
  void sample(); // sample input to create initial skeleton points

  virtual void computeDensity(float h);
  virtual std::vector<IndexVal<float> > findNearestKInputPoints(float px,float py,int k);
  virtual SearchResult radiusSearch(cv::flann::GenericIndex<cvflann::L2<float> >  *index,
                                    cv::Mat features, cv::Mat queries,
                                    cv::Mat radii);
  virtual bool skeletonConverged(std::vector<SkeletonPoint> newX);

  virtual void markBridgePoint(int branchId,BranchEndpointType endType);

  virtual void removeFromCandidates(std::vector<Candidate> &candidates, 
							std::vector<Candidate> branch);

  virtual std::vector<Candidate> traceBranch(std::vector<Candidate> &candidates,
											 Candidate seed,
											 std::vector<Candidate> &noise);
  virtual std::vector<Candidate> traceBranchReverse(std::vector<Candidate> &candidates,
									 Candidate seed,
									 std::vector<Candidate> &noise);
  virtual float getTraceStep();
  virtual void iterationPostProcess(const std::vector<SkeletonPoint> &newX, bool &converged){}
  virtual void processBranch(std::vector<Candidate> &candidates,
					 std::vector<Candidate> branch,
					 std::vector<Candidate> noise);
  virtual void computeProjectionDistance(std::vector<Candidate> &candidates,
								 const Candidate &seed,
								 bool reverse=false);
  virtual void computeBridgePoints(int n);
  virtual bool hasUnvisited(const std::vector<Candidate> &candidates);

  virtual void markVisited(std::vector<Candidate> &candidates, 
				   const std::vector<Candidate> &branch);
    float gaussian(float mean, float stdev);


 //========= protected functions to reimplement ====================
  virtual void processBridgePoints(){}

    void printBranchGraph();

};

//===================== L1Skeleton =============================
/**
 * @brief The L1Skeleton class: extracting a l1 medial skeleton
 * from voting image without directly using point heading
 */
class L1Skeleton : public L1Skeleton_ {
public:

  L1Skeleton(const cv::Mat &voting_image,
			 const std::vector<cv::Mat> &heading_image );

  //============== accessors & utility =========================
  float getHScale() const {return hScale_;}

  std::vector<Candidate> getBranchById(int i){return branches_[i];}
  const std::vector<std::vector<Candidate> > & getBranches(){ return branches_;}
  cv::Mat getSkeletonCovariance(int id){return covariances_[id];}
  std::map<int, std::set<int> > getBranchGraph(){return branchGraph_;}

  virtual void visualizeDensity();
  /**
   * @brief return the candidate with given skeletonId in a given branch
   * if not found, return candidate with id = -1
   */
  // Candidate getCandidateInBranch(int branchId, int skeletonId);

  float getSkeletonLinearity(int i){
      return computeSkeletonAvgSigma(i,6); // k=6 is arbitrary, also used in L1Skeleton::computeBranchCandidate
  }

  /**
   * get the largest principal componenent vector for a 
   * skeleton point
   */ 
  cv::Point2f getPrincipalDirection(int skeletonId);
  //================= skeleton construction==================
  /**
   * @brief perform one step in the fix point iteration 
   * return whether solution has converged
   */
  bool iterate_step();
  void initialize( float skeleton_thresh=-1);
  //void reinitialize(bool incrementH0=false);
  //bool regularize();
  //  void identifyBranchPoints();
  /**
   * @brief trace branch from a skeleton point seed
   * @param skeleton_id id of the seed
   * @return return traced candidates and noise points
   */
  std::pair< std::vector<Candidate>, std::vector<Candidate> > 
	  getBranchFromSeed(int skeleton_id);
  // unused
  void recenterBranches();
  cv::Mat getHeadingDistribution(int x, int y, bool recursive = true);
  cv::Mat getHeadingDistributionNS(int x, int y,bool recursive=true);
  int findStableK(cv::Mat neighbor_ids, cv::Mat dists,
                  int point_id, int k_min, int k_max, int width);
  cv::Mat computeCovMatrix(cv::Mat features, cv::Mat indices,
                           cv::Mat dist, int i);
  float thetaMVN(cv::Point2f x, float h,cv::Mat C);

  //==================== Range/KNN query ===========================
  /**
   * @brief find nearest skeleton point to (px,py) using L2 distance
   * return value: vector of (nearest skeleton index, distance)

  std::vector<IndexVal<float> > findNearestKSkeletonPoints(float px, 
														   float py,
                                                           int k);*/
  //std::vector<IndexVal<float> > findNearestKInputPoints(float px,
//														float py,
    //													int k);

  //SearchResult findSkeletonNeighbors(cv::Point2f p){findSkeletonNeighbors(p,false);}
  SearchResult findSkeletonNeighbors(cv::Point2f p);//, bool directional=false);


  /**
   * @brief find all skeleton points with a fix radius from the queries
   * and have similar heading distributions to the queries
   * (distribution error threshold is defined in parameters.ini)
   * @param index flann search index object (kd tree of features)
   * @param features feature vectors of inputs
   * @param queries feature vectors of queries
   * @param radii  maximum distance for each query
   * return the search result (indices and distances in descending order)
   */

  SearchResult directionalRadiusSearch(
      cv::flann::GenericIndex<cvflann::L2<float> >  *index,
      cv::Mat features, cv::Mat queries,
      cv::Mat radii, bool directional=true);

protected:
  //========= protected members ============================

  std::vector<cv::Mat> covariances_;

  /** index of k nearest neighbors of in the queries **/
  cv::Mat Ip_;
  float hScale_;
  SearchResult skeleton_result_;

  //============= protected methods ============================
  std::pair<float,cv::Mat> sigmaV(cv::Mat indices, cv::Mat dists, int i);

  virtual void computeAlpha(cv::Mat features, SearchResult result);
  virtual void computeBeta(cv::Mat features, SearchResult result);
  /**
   * @brief computeSigma
   * @param features skeleton point
   * @param result nearest neighbor search result
   * update covariances_ list of 2x2 covariance matrices at each skeleton point
   */
  virtual void computeSigma(cv::Mat features, SearchResult result );
  std::pair<float, cv::Mat> sigma(cv::Mat indices, cv::Mat dists, int i,int k);// compute sigma_i for all i in I with k points in neighborhood

  void normalizeByDensity(void);

  std::vector<SkeletonPoint> mergeSkeletonPoints(std::vector<SkeletonPoint> 
                                                 skeletonPts);

  void iterationPostProcess(const std::vector<SkeletonPoint> &newX, bool &converged);

  void filterSearchResultByDirection(cv::Mat features,cv::Mat headingDistQ, 
									 cv::Mat &indices, cv::Mat &dists);
  double evalEnergyFunction(SearchResult inputResult,
                            SearchResult sampleResult,
                            float h);
  void filterSearchResultByDirection2(cv::Mat features,cv::Mat headingDistQ,
                                      cv::Mat &indices, cv::Mat &dists);
  //========= deprecated and debug methods =============================
  cv::Point2f addNoise(cv::Point2f p, float maxError);
  void printSample();
  void printBranchGraph();
  void processBridgePoints();


};

#endif // !L1_SKELETON
