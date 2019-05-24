#pragma once
#ifndef L1Skeleton_HEADING_H
#define L1Skeleton_HEADING_H

//#if (CV_VERSION_EPOCH  > 2)
//#include <opencv2/core.hpp>
//#else
//#include <opencv2/core/core.hpp>
//#endif
#include "l1skeleton.h"

class L1SkeletonWithHeading: public L1Skeleton_{
public:
    enum SigmaFunctionType{
        PCA,
        ENTROPY
    } sigma_type_;

    L1SkeletonWithHeading(const cv::Mat &map_image,
                          const std::vector<cv::Mat> &heading_images  );
    void initialize();

    L1Skeleton_::SearchResult findSkeletonNeighbors(cv::Point2f p);
    cv::Point2f getPrincipalDirection(int i);

   // void reinitialize(bool incrementH0);
    bool  iterate_step();
    //virtual void identifyBranchPoints();
    void clusterSample();
    cv::Mat getSkeletonCovariance(int id){return covariances_[id];}

    void computeBranchCandidate(std::vector<Candidate> &branchCandidates,
                                            bool excludeBranchPoints);
protected:
  //========= protected members ============================
  cv::Mat sampleHeadingDist_; 
  cv::Mat getInputHeadingDistribution(int x, int y, bool recursive=true); //! @brief static across all iteration
  cv::Mat getSampleHeadingDistribution(int i, bool recursive=true); //! @brief changes in each iteration
  int nDirections_;
  SearchResult skeleton_result_;
  std::vector<cv::Mat> covariances_;
  float thetaMVN(cv::Point2f x, float h,cv::Mat C);


  int findStableK(cv::Mat neighbor_ids, cv::Mat dists,
                  int point_id, int k_min, int k_max, int width);

  //========= protected member functions ====================
  virtual void computeSigma(cv::Mat features_, SearchResult result);
  virtual void computeBeta(cv::Mat features_, SearchResult result);

  virtual void computeAlpha(SearchResult result);
  double phi(int i, int j, int a);
  double psi(int i, int ii, int a);
  std::vector<SkeletonPoint> mergeSkeletonPoints(std::vector<SkeletonPoint> skeletonPts);

  std::pair<float,cv::Mat> sigmaPCA(cv::Mat indices, cv::Mat dists, int i, int k);
  std::pair<float,cv::Mat> sigmaEntropy(cv::Mat indices, cv::Mat dists, int i,int k);

 /**
  * \brief merge bridge points that are close to each other
  */
 void processBridgePoints();

};


#endif
