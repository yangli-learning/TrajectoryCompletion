#include "l1skeleton_with_heading.h"
#include "l_shape_iterator.h"
#include "parameters.h"
#include <QDebug>


#if (CV_VERSION_EPOCH  > 2)
#include <opencv2/highgui.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#else
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/flann/dist.h>
#endif

L1SkeletonWithHeading::L1SkeletonWithHeading(const cv::Mat &voting_image,
                                             const std::vector<cv::Mat> &heading_images):
	L1Skeleton_(voting_image,heading_images){
    voting_image_ = voting_image;
    nDirections_ =heading_images.size()/2; //assume nDirections is even
    heading_images_.resize(nDirections_);
    sigma_type_= static_cast<SigmaFunctionType>(Parameters::getInstance()->skeleton_linearity_function);
            //SigmaFunctionType::ENTROPY;
    double min, max;
    for (size_t i=0;i<nDirections_;++i){
        heading_images_[i] = heading_images[i] +heading_images[nDirections_+i];
        cv::minMaxLoc(heading_images_[i], &min, &max);

        qDebug()<< "intialize heading image " << i<< ": max= "
                   <<max <<" min= " << min;

    }

}

void L1SkeletonWithHeading::initialize(){
    L1Skeleton_::initialize();
    qDebug()<<"initialize heading... "  ;
    sampleHeadingDist_ = cv::Mat((int)X_.size(),nDirections_,CV_32FC1,cv::Scalar(0));
    cv::Mat skeleton_radii(X_.size(), 1, CV_32FC1); //radii for query points in X

    for (int i=0;i<(int)X_.size();++i){
        cv::Point2f &pos = X_[i].pos;
        skeleton_radii.at<float>(i,0) = 1.5*getAdaptiveH(pos.x,pos.y);
        for (int a=0;a<nDirections_;++a){
            sampleHeadingDist_.at<float>(i,a) = heading_images_[a].at<float>(pos.y,pos.x);
        }
    }
    std::ostringstream oss;
    oss<<sampleHeadingDist_;
    qDebug()<<oss.str().c_str();
    cvflann::KDTreeIndexParams indexParams( 4 );
    cv::flann::GenericIndex<cvflann::L2<float> > skeletonIndex(skeleton_features_,
                                                                indexParams);
    skeleton_result_ = radiusSearch(&skeletonIndex,skeleton_features_,
                                    skeleton_features_,skeleton_radii);
    oss.str("");
    oss<<skeleton_result_.indices ;
    qDebug()<<"skeleton neighbor search index: " << oss.str().c_str() ;
    covariances_.resize(X_.size());

     computeSigma(skeleton_features_, skeleton_result_);
     computeBeta(skeleton_features_, skeleton_result_);
}

bool L1SkeletonWithHeading::iterate_step(){
    bool converged = false;
    if (X_.size() <= 2)
        return true;

    float mu = Parameters::getInstance()->skeleton_mu;//35;

    cv::Mat queries(X_.size(), 2, CV_32FC1);
    cv::Mat radii(X_.size(), 1, CV_32FC1); //radii for query points in X
    for (int r=0; r<queries.rows; ++r){
        queries.at<float>(r,0) = X_[r].pos.x;
        queries.at<float>(r,1) = X_[r].pos.y;
        radii.at<float>(r,0) = getAdaptiveH(X_[r].pos.x, X_[r].pos.y);
    }
    SearchResult result1 = radiusSearch(inputKNNIndex_,features_,
                                        queries,radii);
    computeAlpha(result1);
    //	findNearestKIn
    std::vector<SkeletonPoint> newX(X_.size());
    SearchResult &result2  = skeleton_result_;
    qDebug() <<"mu: " << mu;
    for (int i=0; i< X_.size();++i){
        newX[i] = X_[i];
        if (X_[i].type == BRANCH || X_[i].type==NOISE){
            continue;
        }
        cv::Point2f xi1, xi2;
        double sum1 =0, sum2=0;
        for (int j=0; j<result1.indices.cols;++j){
            int inputIdx = result1.indices.at<int>(i,j);
            if (inputIdx >= 0){
                cv::Point2f q(features_.at<float>(inputIdx,0),
                              features_.at<float>(inputIdx,1));
                xi1 += q*alpha_.at<float>(i,j);
                sum1 += alpha_.at<float>(i,j);
            }
        }
        cv::Mat indices = result2.indices.colRange(1,result2.indices.cols);
        for (int ip=0; ip< indices.cols;++ip){
            int idx = indices.at<int>(i,ip);
            if (idx >= 0){
                cv::Point2f v = X_[i].pos - X_[idx].pos;
                xi2 += v *beta_.at<float>(i,ip);
                sum2 += beta_.at<float>(i,ip);
            }
        }

        if (sum1 > 1e-5 && sum2 > 1e-5 ){
            cv::Point2f newLoc( xi1*(1/sum1) + mu*sigma_.at<float>(i,0)*xi2*(1/sum2));
            if (newLoc.x < voting_image_.cols && newLoc.x >= 0 &&
                newLoc.y < voting_image_.rows && newLoc.y >= 0){

                // make sure new index is within the bounding box
                newX[i].pos =newLoc;
                }
        }
    }

    X_=mergeSkeletonPoints(newX);
    cvflann::KDTreeIndexParams indexParams( 4 );
    cv::Mat newQueryFeature(X_.size(), 2, CV_32FC1);//todo: initialize
    cv::Mat newRadii(X_.size(), 1, CV_32FC1); //radii for query points in X
    for (int i=0; i<X_.size();++i){
        newRadii.at<float>(i,0) = getAdaptiveH(X_[i].pos.x, X_[i].pos.y);
        newQueryFeature.at<float>(i,0) = X_[i].pos.x;
        newQueryFeature.at<float>(i,1) = X_[i].pos.y;
    }
    cv::flann::GenericIndex<cvflann::L2<float> > index3(newQueryFeature,
                                                       indexParams);
    skeleton_result_ = radiusSearch(&index3,newQueryFeature,
                                   newQueryFeature,newRadii);

    computeSigma(newQueryFeature,skeleton_result_);
    computeBeta(newQueryFeature, skeleton_result_);
    return converged;
}

void L1SkeletonWithHeading::computeBeta(cv::Mat features, SearchResult result){
    cv::Mat &dist = result.distances;
    cv::Mat &indices = result.indices;
    assert(dist.cols > 1);
    beta_ = cv::Mat(dist.rows, dist.cols-1, CV_32FC1, cv::Scalar(0) );
    for (int i=0; i<dist.rows;++i){
        //	int val1 = voting_image_.at<uchar>(X_[i].pos.y, X_[i].pos.x);
         cv::Point2f xi=X_[i].pos;
        float h = getAdaptiveH(X_[i].pos.x, X_[i].pos.y);

        for (int j=1; j<dist.cols;++j){
            if (dist.at<float>(i,j) >= 0){

                float d =sqrt(dist.at<float>(i,j));
                int idx = indices.at<int>(i,j);
                cv::Point2f xj= X_[idx].pos;


                if (d>1e4){
                    beta_.at<float>(i,j-1) = 0.f;
                }else{
                    float psiVal=0.0;
                    for (int a=0;a<nDirections_;++a){
                        psiVal+= psi(i,j,a);
                    }
                    beta_.at<float>(i,j-1) = theta(d,h)/(d*d*psiVal + 1e-4);

                }
            }

        }

    }
    return;

}


std::pair<float,cv::Mat> L1SkeletonWithHeading::sigmaPCA(cv::Mat indices, cv::Mat dists, int i,int k){
    float sigma = 0;
    cv::Mat principalVec (1,2,CV_32F,cv::Scalar(0));
    cv::Point2f xi= X_[i].pos;
    cv::Mat pcadata(k,nDirections_,CV_32F,cv::Scalar(0));
    int nNeighbors = 0;
    for (int ip=1; ip<std::min(k, indices.cols); ++ip){
        int id = indices.at<int>(0,ip);
        if (id<0 || dists.at<float>(0,ip) <0){
            //qDebug() <<"can not compute sigma for sample "<< id <<" ip="<<ip;
            break;
        }
        nNeighbors++;
        getSampleHeadingDistribution(id).copyTo(pcadata.row(ip));
    }
    if (nNeighbors>0){
        cv::PCA pca(pcadata.rowRange(0,nNeighbors),cv::Mat(),0,2);
        cv::Mat eigs = pca.eigenvalues;
        cv::Mat eigv = pca.eigenvectors;
        if (eigs.rows==2){
            float e1 =  eigs.at<float>(0,0), e2 = eigs.at<float>(1,0);

            /*set covariance matrix:
            cv::Mat L = (cv::Mat_<float>(2,2) << e1,0, 0, e2);
            cv::Mat V;
            cv::transpose(eigv,V);
            cv::Mat cov = L*(V*L.inv());
            if (covariances_.size()<=i)
                qDebug()<<"Error: can not find covariance matrix " << i ;
            cov.copyTo(covariances_[i]);
            */
            // compute sigma
            sigma = std::max(e1,e2)/(e1+e2);
            qDebug()<<"e1=" << e1<<" e2="<< e2 <<"sigma: "<< sigma;
            if (e1>e2){
                pca.project(eigv.row(0), principalVec.row(0));
            }else{
                pca.project(eigv.row(1), principalVec.row(0));
            }
        }else if (eigs.rows==1){
            sigma = 1;
            std::ostringstream oss;
            oss << eigv;

            qDebug()<<"single eigenvalue " << eigs.at<float>(0,0) << "from "
                   << nNeighbors << " points with eigen values"
                      << oss.str().c_str();
             //zero or 1?   !@todo need to check with data
        }
        qDebug()<<"eigenvalue of sample "<< i <<" has size " <<eigs.rows <<"x"<<eigs.cols
                  <<" sigma = " << sigma;
    }
    return std::make_pair(sigma , principalVec); // later should try the sample heading direction
}
float L1SkeletonWithHeading::thetaMVN(cv::Point2f x, float h,cv::Mat C){
    cv::Mat Xt = (cv::Mat_<float>(1,2) << x.x,x.y);
    cv::Mat X;
    cv::transpose(Xt,X);// = (cv::Mat<float>(2,1) << x.x ,x.y );
    //cv::Mat Cinv;
    //C.inv().convertTo(Cinv,CV_32F) ;
    cv::Mat product = -Xt* (C.inv()*X) ;
    return exp(product.at<float>(0)/(2*h));
}

std::pair<float,cv::Mat> L1SkeletonWithHeading::sigmaEntropy(cv::Mat indices,
                                                             cv::Mat dists,
                                                             int i, int k){
    float sigma = 0;
    cv::Mat principalVec(1,2,CV_32F,cv::Scalar(0));
    cv::Point2f xi= X_[i].pos;
    cv::Mat headingDistrib(1,nDirections_,CV_32F,cv::Scalar(0));
    int nNeighbors = 0;
    for (int ip=1; ip<std::min(k, indices.cols); ++ip){
        int id = indices.at<int>(0,ip);
        if (id<0 || dists.at<float>(0,ip) <0){
            //qDebug() <<"can not compute sigma for sample "<< id <<" ip="<<ip;
            break;
        }
        nNeighbors++;
        headingDistrib += getSampleHeadingDistribution(id);

    }
    //std::ostringstream oss; oss << headingDistrib;
   // qDebug()<<"heading distrib" << oss.str().c_str();
    float sum = cv::sum(headingDistrib)[0];
    cv::Mat headingDistribN(1,nDirections_,CV_32F,cv::Scalar(0));
    if (sum >0){
        headingDistribN = headingDistrib/sum;
    }
    //oss.str("");
    //oss << headingDistribN;
    //qDebug() <<" normalized: "<< oss.str().c_str();
    // compute sigma as the entropy of headingDistribN
    for (int j=0; j<nNeighbors;++j){
        sigma += -1*headingDistribN.at<float>(0,j)*log(headingDistribN.at<float>(0,j));

    }
    sigma = (log(nDirections_)-sigma)/log(nDirections_); //sigma is higher when direction entropy is less
    double maxVal,minVal; cv::Point maxLoc,minLoc;
    minMaxLoc(headingDistribN,&minVal,&maxVal,&minLoc,&maxLoc);
    float angleDegree = (nDirections_-maxLoc.x)*float(360/(2*nDirections_));

    cv::Point2f p =
       LShapeIterator::UnitDirection(angleDegree);
    principalVec.at<float>(0,0) = p.x;
    principalVec.at<float>(0,1) = p.y;

    //set covariance matrix:
    /*if (covariances_.size()<=i){
       qDebug()<<"Error: can not find covariance matrix " << i ;
    }*/
    //cv::Mat SR( 2, 2, CV_32F,cv::Scalar(0));
    /*
    cv::Mat SR64 = (cv::getRotationMatrix2D(cv::Point2f(0,0),
                                          angleDegree,sqrt(k) )).colRange(0,2);
    qDebug()<<"\nangle degree="<<angleDegree ;
    cv::Mat SR32;
    SR64.convertTo(SR32,CV_32F);
    cv::Mat SRt;
    cv::transpose(SR32,SRt);
    cv::Mat cov= SR32*SRt;
    */
    /*
    double a =  (p.y==0)? -1*M_PI: -1*atan(p.y/p.x);//M_PI*(180-angleDegree)/180.f;
    cv::Mat S = (cv::Mat_<float>(2,2) <<  sqrt(k), 0, 0,sqrt(k)),
            R = (cv::Mat_<float>(2,2) << cos(a),-sin(a),sin(a),sin(a));

    cv::Mat T = R*S, Tt;
    cv::transpose(T,Tt);
    cv::Mat cov = T*Tt;
    cov.copyTo(covariances_[i]);
    */
    return std::make_pair(sigma , principalVec);
}
int L1SkeletonWithHeading::findStableK(cv::Mat neighbor_ids, cv::Mat dists,
                                int point_id, int k_min, int k_max, int width  ){
    float maxSigma = 0;
    float optKl =k_min;
    int kl, ku, k;

    for ( kl = k_min; kl < k_max-width+1; ++kl){
        bool stable = true; float stableSigma = 0.f;
        for ( k=kl;k< kl+width; ++k){
            if (sigma_type_ == SigmaFunctionType::PCA ){
                stableSigma = sigmaPCA(neighbor_ids, dists, point_id, k).first;
            }else{
                stableSigma = sigmaEntropy(neighbor_ids,dists, point_id,k).first;
            }
            float thresh = (sigma_type_==SigmaFunctionType::ENTROPY)?0.7:0.85;
            if (k <thresh){
                stable =false;
                break;
            }
        }
        if (stable && stableSigma > maxSigma){
            maxSigma = stableSigma;

            optKl = kl;//floor(kl+width/2);
        }
    }
    return floor(optKl+width/2);
}

void L1SkeletonWithHeading::clusterSample(){
    // input sigma and function computeSkeletonAvgSigma
    // reimplement computeAvgSigma to add in weight
    std::vector<float> votes(X_.size() ,0);
    for (size_t s=0;s<votes.size();++s){
        votes[s] = computeSkeletonAvgSigma(s,5 ); //k=5 is not used!!
    }

    // sort all samples by sigma value
    // compute clusters
   /*
    int height= voting_image_.rows;
    int width = voting_image_.cols;
    double minV,maxV;
    cv::minMaxLoc(sigma_,&minV, &maxV);
    QImage sigmaImg(width, height,QImage::Format_RGB32);
    sigmaImg.fill(QColor(Qt::white).rgb());
    for (size_t i=0;i<X_.size();++i){
        float val = 1 - getSkeletonLinearity(i);
        osg::Vec4 color = ColorMap::getInstance().getContinusColor(val,0,1,true);
        sigmaImg.setPixel(X_[i].pos.x, X_[i].pos.y,
                          qRgb( color[0]*255,color[1]*255,color[2]*255));
    }
    QLabel *l = new QLabel;
    l->setScaledContents(true);
    l->setPixmap(QPixmap::fromImage(sigmaImg));
    l->show();
    */
}

void L1SkeletonWithHeading::computeSigma(cv::Mat features, SearchResult result){
    assert(result.indices.cols > 1);
     cv::Mat indices  = result.indices.colRange(1, result.indices.cols);
    cv::Mat dist = result.distances.colRange(1,result.distances.cols);
     sigma_ = cv::Mat(indices.rows, 1, CV_32F, cv::Scalar(0));
    pca_directions_ = cv::Mat(indices.rows,2,CV_32F,cv::Scalar(0));

    for (int i=0;i<indices.rows;++i){
        int k =   findStableK(indices.row(i),dist.row(i),i,4,15,4);
        qDebug()<<"Found stable k: " << k;
        if (k>0){
            if (sigma_type_ == SigmaFunctionType::PCA ){
                auto stableSigma = sigmaPCA(indices.row(i),dist.row(i),i,k);
                sigma_.at<float>(i,0) = stableSigma.first;
                (stableSigma.second).copyTo(pca_directions_.row(i)); //this is bad, not used
            }else if (sigma_type_ == SigmaFunctionType::ENTROPY){
                auto stableSigma = sigmaEntropy(indices.row(i),dist.row(i),i,k);
                sigma_.at<float>(i,0) = stableSigma.first;
                (stableSigma.second).copyTo(pca_directions_.row(i)); //this is bad, not used
            }
        }
        qDebug()<<"sigma "<< i <<" = " << sigma_.at<float>(i,0);
    }
}
void L1SkeletonWithHeading::computeAlpha( SearchResult result){

    // compute alpha(i,j) where xi is a sample point in X, qj is a input point in Q
    // alpha(i,j) = theta(||xi-qj||)/(||xi-qj||*phi(xi,qj)

    cv::Mat &dist = result.distances;
    cv::Mat &indices = result.indices;

    alpha_ = cv::Mat(dist.rows, dist.cols, CV_32FC1,  cv::Scalar(0));
    for (int i=0; i<dist.rows;++i){
        //	int val1 = voting_image_.at<uchar>(X_[i].pos.y, X_[i].pos.x);
        //	float intensity1 = float(256-val1)/256.f;
        float h = getAdaptiveH(X_[i].pos.x, X_[i].pos.y);
        cv::Point2f xi = X_[i].pos;

        for (int j=0; j<dist.cols;++j){
            if (dist.at<float>(i,j) >= 0){
                float d = sqrt(dist.at<float>(i,j));

                int idx = indices.at<int>(i,j);
                //cv::Point2f xj(features_.at<float>(idx,0),features_.at<float>(idx,1));

                if (d>1e7){
                    alpha_.at<float>(i,j) = 0.f;
                }else{
                    float phiVal = 0.0;
                    for (int a=0;a<nDirections_;++a){
                        phiVal+= phi(i,j,a);
                    }
                    alpha_.at<float>(i,j) = theta(d,h)/(d*phiVal+1e-4);
                }


            }

        }
    }
}
cv::Mat L1SkeletonWithHeading::getInputHeadingDistribution(int x, int y,
                                                           bool recursive){
    bool allZero =true;
    cv::Mat headingDist(1,nDirections_,CV_32F,cv::Scalar(0));

    for (size_t a=0; a< nDirections_;++a){
        headingDist.at<float>(0,a) = heading_images_[a].at<float>(y,x);
        if (headingDist.at<float>(0,a)>1e-5)
            allZero = false;
    }
    if (allZero && recursive){
        // use weighted neighbors
        auto nlist = findNearestKInputPoints(x,y,1);
        if (nlist.size()>0){
            IndexVal<float> nearest = nlist.front();
            cv::Point2f p(features_.at<float>(nearest.index,0),
                      features_.at<float>(nearest.index,1));
            headingDist = getInputHeadingDistribution(p.x, p.y, false);
        }

    }
    return headingDist;
}

cv::Mat L1SkeletonWithHeading::getSampleHeadingDistribution(int i, bool recursive){
	if (i>=0 && i<sampleHeadingDist_.rows){
        cv::Mat dist = sampleHeadingDist_.row(i);
        if (countNonZero(dist) < 1){
            dist = getInputHeadingDistribution(X_[i].pos.x,X_[i].pos.y, true);
        }
        return dist;
	}else{
		qDebug()<<"Error: can not get heading distribution of sample "<<i;
		return cv::Mat(1 ,nDirections_ ,CV_32F,cv::Scalar(0));
	}
}

double  L1SkeletonWithHeading::phi(int i, int j, int a){
    //compute l2 error of heading distribution pixel i and pixel j
    cv::Point2f vi = X_[i].pos;
    cv::Mat hdist_i = getSampleHeadingDistribution(i),
            hdist_j = getInputHeadingDistribution(vi.x, vi.y);
    return sqrt((hdist_i-hdist_j).dot(hdist_i-hdist_j));
}

double  L1SkeletonWithHeading::psi(int i, int ii, int a){
    //compute l2 error of heading distribution pixel i and pixel j
    cv::Point2f vi = X_[i].pos;
    cv::Mat hdist_i = getSampleHeadingDistribution(i),
            hdist_ii = getSampleHeadingDistribution(ii);
    return sqrt((hdist_i-hdist_ii).dot(hdist_i-hdist_ii));
}
L1Skeleton_::SearchResult L1SkeletonWithHeading::findSkeletonNeighbors(cv::Point2f p){
    //SearchResult result;
    int n = X_.size();
    cv::Mat features(n, 2, CV_32F);
    for (int i=0; i<n;++i){
        features.at<float>(i,0) = X_[i].pos.x;
        features.at<float>(i,1) = X_[i].pos.y;
    }
    float radius = getAdaptiveH(p.x,p.y);
    qDebug()<<"Search using r="<<  radius;

    cvflann::KDTreeIndexParams indexParams;
    cv::flann::GenericIndex<cv::flann::L2<float> > index(features,
                                                         indexParams);
    cv::Mat query = (cv::Mat_<float>(1, 2) << p.x, p.y );
    cv::Mat radii = (cv::Mat_<float>(1,1) << radius);
    return radiusSearch(&index, features, query,radii);
}
cv::Point2f L1SkeletonWithHeading::getPrincipalDirection(int i){
  /*  if (int i>=0 && i<pca_directions_.rows){
        return cv::Point2f(pca_directions_.at<float>(i,0),
                           pca_directions_.at<float>(i,1));
    }else{*/
    double maxV,minV; cv::Point maxLoc,minLoc;
    minMaxLoc(sampleHeadingDist_.row(i),&minV, &maxV, &minLoc, &maxLoc);
    qDebug()<< "max location: " << maxLoc.x <<","<< maxLoc.y ;
    cv::Point2f p =
        LShapeIterator::UnitDirection((nDirections_-maxLoc.x)*float(360/(2*nDirections_)));

    return p;
    /*}*/
        //return the direction of max vote?
}


std::vector<SkeletonPoint> L1SkeletonWithHeading::mergeSkeletonPoints(
        std::vector<SkeletonPoint> skeletonPts){
    //!@todo
    std::map<std::pair<int,int>, std::vector<int> > gridFilter;
    float r = Parameters::getInstance()->skeleton_sample_cell_size *0.25f;
    for (int i=0; i<skeletonPts.size();++i){
        std::pair<int,int> key;
        key.first = floor(skeletonPts[i].pos.x / r);
        key.second = floor(skeletonPts[i].pos.y/r);
        gridFilter[key].push_back(i);
    }
    std::vector<SkeletonPoint> result;
    cv::Mat newSampleHeadingDist(sampleHeadingDist_.rows, sampleHeadingDist_.cols,CV_32F,cv::Scalar(0));

    for (auto it=gridFilter.begin();it!=gridFilter.end();++it){
        int groupSize = it->second.size();
        int idx = result.size();
        if (groupSize ==1){
            result.push_back(skeletonPts[it->second[0]]);
            getSampleHeadingDistribution(it->second[0]).copyTo(newSampleHeadingDist.row(idx)) ;
        }else if (groupSize >1){
            cv::Mat heading(1,nDirections_,CV_32F,cv::Scalar(0));
            // compute the average position of t
            // push the first into result, modify its position
            cv::Point2f avg(0,0);
            bool hasBranchPts = false, hasBridgePts = false;
            int bridgeId=-1, branchId=-1;
            for (int j=0; j < groupSize;++j){
                SkeletonPoint pt = skeletonPts[it->second[j]];
                heading += getSampleHeadingDistribution(it->second[j]);

                if (pt.type == SkeletonPointType::BRIDGE){
                    hasBridgePts = true;
                    bridgeId = pt.id;
                }else if (pt.type == SkeletonPointType::BRANCH){
                    hasBranchPts = true;
                    branchId = pt.id;
                }

                avg += pt.pos;
            }
            avg *= 1.f/float(it->second.size());
            SkeletonPoint s = skeletonPts[it->second[0]];
            s.pos = avg;
            s.type  = (hasBridgePts)? SkeletonPointType::BRIDGE :
                ( (hasBranchPts)? SkeletonPointType::BRANCH: s.type);
            s.id  = (hasBridgePts)?bridgeId:
                ( (hasBranchPts)? branchId: s.id);
            heading.copyTo(newSampleHeadingDist.row(idx));
            result.push_back(s);
        }


    }
    sampleHeadingDist_ = newSampleHeadingDist.rowRange(0,result.size());
    //update sampleHeadingDist_(!)
    qDebug() << "******\nbefore: " << skeletonPts.size() << " points, after: " <<
        result.size() << " points.";

    return result;
}


void L1SkeletonWithHeading::computeBranchCandidate(std::vector<Candidate> &branchCandidates,
                                        bool excludeBranchPoints){

    float sigmaThresh; //= Common::percentile(sigma_,p );
    if (sigma_type_ == SigmaFunctionType::PCA){
        sigmaThresh = Parameters::getInstance()->skeleton_sigma_thresh ;
    }else{
        float p = 0.85;//Parameters::getInstance()->skeleton_sigma_thresh ;
        sigmaThresh = Common::percentile(sigma_,p );

    }
    qDebug() << "branch candidate threshold: " << sigmaThresh;
    int k = 5;
    //qDebug() << "\n<<< start candidate <<<" << endl;
    int nB = 0;
    for (int i=0;i< X_.size();++i){
    // existing branch points or noise points can not be candidates
    if ( (X_[i].type == SkeletonPointType::BRANCH
          || X_[i].type == SkeletonPointType::NOISE) && excludeBranchPoints){
        continue;

    }
    nB++;
    // get 5 nearest neighbor of xi in X_
    //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
    float avgSigma = computeSkeletonAvgSigma(i,k);//!!! this needs to only consider non-branch pts ???
    //!@@@ todo sigma_.at<float>(i,0);//
    //>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
    if (avgSigma > sigmaThresh){
        branchCandidates.push_back(Candidate(X_[i], getPrincipalDirection(i),
                                             avgSigma, Common::inf, false));

    }
}
qDebug() << "total skeleton points: " << X_.size() <<" # non-branch pts: "
         << nB << " # candidates: " << branchCandidates.size() << endl;
}
//>>>>>>>>>>>>>>>>>>>todo>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
/*

bool L1SkeletonWithHeading::hasUnvisited(const std::vector<Candidate> &candidates){
    // count # of unvisited candidates
    int totalUnvisited = candidates.size();
    for (int k=0;k<candidates.size();++k){
        if (candidates[k].visited){
            totalUnvisited--;
        }
    }
    if (totalUnvisited==0)
        return false;

    return (std::find_if(candidates.begin(),candidates.end(),unvisited())
            !=candidates.end());
}


void L1SkeletonWithHeading::identifyBranchPoints(){
    // intialize branch point vector
    std::vector<Candidate> branchCandidates;
    computeBranchCandidate(branchCandidates,true);

    // identify branch point from candidates
    int existingNumBranches = branches_.size();
    while (hasUnvisited(branchCandidates)){
        std::vector<Candidate>::iterator sit = max_element(branchCandidates.begin(),
                                       branchCandidates.end(),
                                       compareSigma()); // check unvisited!
        Candidate seed = *sit;
        sit->visited = true;
        std::vector<Candidate> noise;
        std::vector<Candidate> branchForward = traceBranch(branchCandidates, seed,noise),
            branch = traceBranchReverse(branchCandidates,seed,noise);
        branch.pop_back();
        branch.insert(branch.end(),branchForward.begin(),branchForward.end());
        processBranch(branchCandidates,branch,noise);
    }
    computeBridgePoints(existingNumBranches);
    qDebug() <<"finsihed processing all branches";
    /// exclude branch points from X_
    // std::vector<SkeletonPoint> XNew;
    // for (size_t i=0;i<X_.size();++i){
    //     //		if (branchCandidates
    //         if (X_[i].type!= SkeletonPointType::BRANCH){
    //         XNew.push_back(X_[i]);
    //     }
    //     }
    // // note: all type of X_ is still NONE since TraceBranch
    // // only changes branchCandidates
    // // identify bridge points in the new X_: [To be implemented]
    // X_ = XNew;
    ///


}

void L1SkeletonWithHeading::computeBranchCandidate(std::vector<Candidate> &branchCandidates,
                                                    bool excludeBranchPoints){
    float sigmaThresh = Parameters::getInstance()->skeleton_sigma_thresh;
    int k = 6;
    //qDebug() << "\n<<< start candidate <<<" << endl;
    int nB = 0;
    for (int i=0;i< X_.size();++i){
        // existing branch points or noise points can not be candidates
        if ( (X_[i].type == SkeletonPointType::BRANCH
              || X_[i].type == SkeletonPointType::NOISE) && excludeBranchPoints){
            continue;
        }
        nB++;
        // get 5 nearest neighbor of xi in X_
        //<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
        float avgSigma = sigma_.at<float>(i,0);//computeSkeletonAvgSigma(i,k);

        if (avgSigma > sigmaThresh){
            branchCandidates.push_back(Candidate(X_[i], getPrincipalDirection(i),
                                                 avgSigma, Common::inf, false));
            //	qDebug() << X_[i].pos.x<< ", "<< X_[i].pos.y ;
        } //else{
          //  qDebug() << "candidate " << i<<": "
          //           << X_[i].id<<" sigma below threshold " << avgSigma ;
          //           }
    }
    //qDebug ()<< "<<<<<<<<<<<<<<<" << endl;
    qDebug() << "total skeleton points: " << X_.size() <<" # non-branch pts: "
             << nB << " # candidates: " << branchCandidates.size() << endl;
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

void L1SkeletonWithHeading::computeBridgePoints(int existingNumBranches){
    std::vector<int> headBridgePts, tailBridgePts;

    if (X_.size() >= 5){
        for (int b= existingNumBranches;b<branches_.size();++b){
            int bid = markBridgePoint(b,BranchEndpointType::FIRST); // find bridge point of first endpoint
            if (bid >=0) headBridgePts.push_back(bid);
            bid = markBridgePoint(b,BranchEndpointType::LAST); // find bridge point of second endpoint
            if (bid >=0) tailBridgePts.push_back(bid);
        }
        processBridgePoints();
    }
    // merge bridge points
    //debug: output current branchgraph
    //printBranchGraph();
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void L1SkeletonWithHeading::processBranch(std::vector<Candidate> &candidates,
                               std::vector<Candidate> branch,
                               std::vector<Candidate> noise){

    //retain branch if its length is at least 5
    if (branch.size() >= 4){

        int branchId = branches_.size();

        for (size_t i=0;i<branch.size();i++){

            // mark each candidate in the branch as BRANCH
            std::vector<Candidate>::iterator cit = searchCandidates(candidates,branch[i].id);
            cit->type =  SkeletonPointType::BRANCH;
            int id = branch[i].id;
            auto vit =find_if( X_.begin(), X_.end(), [&id](const SkeletonPoint &p)
                               {return p.id == id;});
                vit->type = SkeletonPointType::BRANCH;

        }
        branches_.push_back(branch);
        markVisited(candidates,branch);
        removeFromCandidates(candidates,noise);

        removeFromCandidates(candidates,branch);
    }
}
void L1SkeletonWithHeading::markVisited(std::vector<Candidate> &candidates,
                             const std::vector<Candidate> &branch){
    std::vector<Candidate>::const_iterator it;
    std::vector<Candidate>::iterator found;
    for (it=branch.begin();it!=branch.end();++it){

        found = find(candidates.begin(), candidates.end(),*it);
        if (found!=candidates.end()){
            found->visited = true;
        }
    }
}

/>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
void L1SkeletonWithHeading::computeProjectionDistance(std::vector<Candidate> &candidates,
                                           const Candidate &seed,
                                           bool reverse){
    for (size_t i=0;i<candidates.size();++i){

        if (candidates[i]== seed){
            candidates[i].proj_dist =0;
        }else{
            cv::Point2f vec = candidates[i].pos - seed.pos;

            if (!reverse){
                float dist =vec.dot(seed.dir);
                candidates[i].proj_dist = (dist>=0)? dist:Common::inf;
            }else{
                float dist = vec.dot(-seed.dir);
                // reverse projection distance
                //					qDebug()<<"projdist: " << dist;
                candidates[i].proj_dist = (dist>=0)? dist :Common::inf;
            }
        }
    }
}
//>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/
/*
std::vector<Candidate> L1SkeletonWithHeading::traceBranchReverse(
      std::vector<Candidate> &candidates,
      Candidate seed,
      std::vector<Candidate> &noise){
    float res =Parameters::getInstance()->skeleton_max_trace_step;
    float maxTraceStep = (h0_/Parameters::getInstance()->skeleton_h0)*res;
    float branchCosThresh = Parameters::getInstance()->skeleton_branch_cos_thresh;

    Candidate x0=seed;
    cv::Point2f v0 = -(seed.dir);
    std::vector<Candidate> branch;
    int max_branch_size = 40;
    bool branchEnded;

    do{
        branch.insert(branch.begin(),x0);//		branch.push_back(x0);
        computeProjectionDistance(candidates,x0,true);//reverse trace
        int neighborSize = std::min(30, (int)candidates.size());
        partial_sort(candidates.begin(),candidates.begin()+neighborSize,
                     candidates.end()); //sort by projection distance
        std::vector<Candidate> neighbors(candidates.begin(),
                                         candidates.begin()+neighborSize);
        std::vector<Candidate>::iterator it;
        branchEnded = true;
        for (it=neighbors.begin();it!=neighbors.end();it++){
            if (x0 == *it  )
                continue;
            cv::Point2f vec = it->pos - x0.pos;

            if (it->pos == x0.pos){
                qDebug() << "duplicates: point " << it->id <<" and point "
                         << x0.id ;
            }
            if (cos(angleBetween(vec,-v0)) <= branchCosThresh
                && cv::norm(vec) < maxTraceStep){
                v0 = vec;
                x0 = *it;
                branchEnded = false;
                break;
            }else {
                //!!@@@ add it as noise if it is within
                if ( cv::norm(vec) < maxTraceStep*1.f  )
                  noise.push_back(*it);
            }
        }

    }while(!branchEnded || branch.size() > max_branch_size);

    return branch;
}
std::vector<Candidate> L1SkeletonWithHeading::traceBranch(std::vector<Candidate> &candidates,
                                               Candidate seed,
                                               std::vector<Candidate> &noise){
    	// float maxTraceStep =(Parameters::getInstance()->skeleton_max_trace_step
        //                  + hScale_)*(Parameters::getInstance()->skeleton_sample_cell_size);// h0_*0.5;
    
    float res =Parameters::getInstance()->skeleton_max_trace_step;
    float maxTraceStep = (h0_/Parameters::getInstance()->skeleton_h0)*res;

  
    bool branchEnded;
    float branchCosThresh = Parameters::getInstance()->skeleton_branch_cos_thresh;
    Candidate x0=seed;
    cv::Point2f v0 = seed.dir;
    std::vector<Candidate> branch;
    int max_branch_size = 40;
    do{
        if (branch.size()>0 && branch.back().pos ==x0.pos){
            //do nothing
        }else{
            branch.push_back(x0);
        }

        computeProjectionDistance(candidates,x0);
        int neighborSize = std::min(30, (int)candidates.size());
        partial_sort(candidates.begin(),candidates.begin()+neighborSize,
					 candidates.end()); //sort by projection distance
        std::vector<Candidate> neighbors(candidates.begin(),
										 candidates.begin()+neighborSize);
        std::vector<Candidate>::iterator it;
        branchEnded = true;
        for (it=neighbors.begin();it!=neighbors.end();it++){
            if (x0 == *it  )
                continue;
            cv::Point2f vec = it->pos - x0.pos;

            if (it->pos == x0.pos){
                qDebug() << "duplicates: point " << it->id <<" and point "
                         << x0.id ;
            }
            if (cos(angleBetween(vec,-v0)) <= branchCosThresh
                && cv::norm(vec) < maxTraceStep){
                v0 = vec;
                x0 = *it;
                branchEnded = false;
                break;
            }else {
                // add it as noise if it is within
                if ( cv::norm(vec) < maxTraceStep*0.1f  )
                    noise.push_back(*it);
            }
        }

    }while(!branchEnded || branch.size() > max_branch_size);
    
      // if (noise.size()>=1){
      //   noise.pop_back();
      //   }

    return branch;
	}*/
//==================



void L1SkeletonWithHeading::processBridgePoints(){
    // merge bridge points that are close to each other
    std::map<int, std::set<int> >::iterator iter;

    float res =Parameters::getInstance()->skeleton_max_trace_step;
    float r =  2*res;

    qDebug()<<"merge nearby bridge points (total: " << branchGraph_.size() << ").";
    // iterate over all bridge point and merge nearby points
    for (iter=branchGraph_.begin();iter!=branchGraph_.end();++iter){
        // find bridge point in X_. if not found,continue;
        int bridgeId = iter->first;
        auto bridgePt = find_if(X_.begin(),X_.end(),[&bridgeId]
                                (const SkeletonPoint &v){
                                    return v.id == bridgeId;});
        if (bridgePt== X_.end())
            continue;
        // find neighbors of bridge point
        SearchResult result = findSkeletonNeighbors( bridgePt->pos);

        std::vector<SkeletonPoint> bridgePts;
        cv::Point2f avg;
        // iterate over each neighbor
        for (int i=0; i<result.indices.cols;++i){
            int neighborId = result.indices.at<int>(0,i);//correct index
            if (neighborId <0) continue;

            SkeletonPoint &neighbor = X_[neighborId];//ith neighbor
            if (neighbor.id == bridgeId) continue;

            if( neighbor.type==SkeletonPointType::BRIDGE){
                // merge bridge point together into one bridge point
                bridgePts.push_back(neighbor);
                avg += neighbor.pos;
                break;
            }
        }
		// merge neighbor points
        if (bridgePts.size()>0){
            bridgePt->pos = avg*float(1.f/bridgePts.size());
            // merge all nearby bridge points
            for (unsigned int b=0;b<bridgePts.size();++b){
                // remove merged bridgePoint
                int neighborId = bridgePts[b].id;
                auto bridgePtb = find_if(X_.begin(),X_.end(),[&neighborId]
										 (const SkeletonPoint &v){
											 return v.id == neighborId;
										 });
                if (bridgePtb != X_.end()){
                    bridgePtb->type = SkeletonPointType::NOISE;
                }
                std::set<int> branchIds1 = iter->second;
                //sort(branchIds1.begin(),branchIds1.end());

                // remove bridge point, merge its adacent list
                std::map<int, std::set<int> >::iterator iter2
                    =branchGraph_.find(bridgePts[b].id);

                if (iter2==branchGraph_.end()){
                    qDebug() << "skip unfound branch point " << iter2->first;
                    continue; //
                }
                std::set<int> branchIds2 =iter2->second;
                //sort(branchIds2.begin(),branchIds2.end());

                iter->second.insert(branchIds2.begin(),branchIds2.end());


                branchGraph_.erase(iter2);

                qDebug() <<"merging branches that contain neighbor bridge point "<< iter2->first
                         << "["<< Common::printIntSet(branchIds1).c_str()<<"]";
                qDebug() <<	" with branches that contain "<< iter->first <<
                "["<<	Common::printIntSet(branchIds2).c_str()<<"]";
                qDebug()<<"result: "<< "["<<Common::printIntSet(iter->second).c_str()<<"]";
            }
        }
    }
}

