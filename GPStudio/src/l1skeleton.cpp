#include <cassert>
#include <QDebug>
#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <sstream>
#include <string>
#include <opencv2/core/core_c.h>
#include <opencv2/core/version.hpp>
#if (CV_VERSION_EPOCH  > 2)
#include <opencv2/imgcodecs.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#else
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#include <opencv2/flann/dist.h>
#endif
#include "common.h"
#include "parameters.h"
#include "l1skeleton.h"
#include <QImage>
#include <QColor>
#include <QLabel>
#include <QPainter>
#include <QVBoxLayout>
#define _USE_MATH_DEFINES
#include <cmath>
#include "color_map.h"
//==================================================================
// definition for L1Skeleton_
//==================================================================
L1Skeleton_::L1Skeleton_(const cv::Mat &voting_image,
						 const std::vector<cv::Mat> &heading_images):
	voting_image_(voting_image){
    inputKNNIndex_=0;
	iter_=0;

}
L1Skeleton_::~L1Skeleton_(){
    delete inputKNNIndex_;
}
void L1Skeleton_::sampleUniform(){
    int resolution = Parameters::getInstance()->skeleton_sample_cell_size;
    int maxN = Parameters::getInstance()->skeleton_n_samples_per_cell;

    int height = ceil(voting_image_.rows/resolution);
    int width = ceil(voting_image_.cols/resolution);
    int index = 0;

    for (int r=0;r<height;++r){
        for (int c=0; c<width;++c){

            for (int rr=r*resolution;rr<(r+1)*resolution;++rr){
                for (int cc=c*resolution;cc< (c+1)*resolution;++cc){
                    if (Q_.at<uchar>(rr,cc)==1){
                        I_.push_back(SkeletonPoint{index++,
                                    cv::Point((c+0.5)*resolution, (r+0.5)*resolution),
                                    SkeletonPointType::NONE});
                        goto end_loop;
                    }
                }
            }
        end_loop:
            static_cast<void> (0) ;
        }
    }
    qDebug() <<"Total # of cells: " << height*width << endl;
    qDebug() <<"Sampled " << I_.size() <<" points out of "
             << inputSize_+I_.size()
             << " (%" <<float(I_.size())/ float(inputSize_+I_.size()) *100
             <<")" << " max index = " << index <<  endl;
    X_ = I_;
    return;
}

void L1Skeleton_::initialize(float vmin){
	qDebug() <<"Call L1Skeleton_::initialize for vmin="<< vmin <<endl;
	
	float skeleton_thresh = vmin;
	double minV, maxV;
    cv::minMaxLoc(voting_image_, &minV, &maxV);
	float thresh = std::max(5.f , float(maxV)*skeleton_thresh);//24
	qDebug() <<"Skeleton max thresh " << maxV;
	qDebug() << "Initialize skeleton points with minimum threshold value " << thresh;
    Q_=cv::Mat(voting_image_.rows,voting_image_.cols,CV_8UC1, cv::Scalar(0));
    inputSize_ = 0;
    for (size_t y=0;y<voting_image_.rows;++y){
        //		uchar* p_color_row = voting_image_.ptr<float>(y);
        for (size_t x=0; x<voting_image_.cols;++x){
            if ( voting_image_.at<float>(y,x) > thresh){
                Q_.at<uchar>(y,x) = 1;
                inputSize_++;
            }
        }
    }
	qDebug() << "Input size: " << inputSize_ << endl;
	
    sampleUniform();
   // sample();
    // initialize feature matrix for all input points Q
    features_ = cv::Mat(inputSize_, 2, CV_32FC1);
    int i=0;
    for (int r=0;r<Q_.rows;++r){
        for (int c=0; c< Q_.cols;++c){
            if (Q_.at<uchar>(r,c)==1){
                features_.at<float>(i,0) = c; // x
                features_.at<float>(i,1) = r; // y
                i++;
            }
        }
    }

    // build KNN feature search tree for all input points Q
    cvflann::KDTreeIndexParams indexParams( 4 );
    inputKNNIndex_ =
        new cv::flann::GenericIndex<cvflann::L2<float> >(features_,
                                                         indexParams);
    //float diag = sqrt(pow(voting_image_.rows,2)+pow(voting_image_.cols,2));
    h0_ = Parameters::getInstance()->skeleton_h0;
	sigma_thresh_=	Parameters::getInstance()->skeleton_sigma_thresh;
	std::cout <<"sigma thresh: " << sigma_thresh_<< std::endl;
    computeDensity(h0_);
    // initialize skeleton feature matrix
    skeleton_features_=cv::Mat(X_.size(), 2, CV_32FC1);
    for (size_t i=0; i<X_.size();++i){
        skeleton_features_.at<float>(i,0) = X_[i].pos.x;
        skeleton_features_.at<float>(i,1) = X_[i].pos.y;
       // skeleton_radii.at<float>(i,0) = getAdaptiveH(X_[i].pos.x, X_[i].pos.y);
    }
	qDebug()<<"skeleton feature matrix has size " << skeleton_features_.rows<<"x"
			<< skeleton_features_.cols ;
}
void L1Skeleton_::computeDensity(float h){
    //compute density as weighted average intensity
    cv::Mat tmp;
    double sigma = h;// threshold = 5, amount = 1;
    cv::GaussianBlur(voting_image_, tmp, cv::Size(),sigma,sigma);

    density_ = cv::Mat(inputSize_,1,CV_32FC1,cv::Scalar(1));

    double min, max;
    cv::minMaxLoc(tmp, &min, &max);//voting_image_, &min, &max);
    if (max <= 0)
        return;

    for (int i=0;i < inputSize_;++i){
        int x = features_.at<float>(i,0);
        int y = features_.at<float>(i,1);

        float val = tmp.at<float>(y,x);//voting_image_.at<float>(y,x);//tmp.at<cv::Vec3b>(y,x)[0];
        float intensity =10* float(float(val-min)/max);
        density_.at<float>(i,0) = intensity;
        //	qDebug() << "density of Q["<< i << "] = " << intensity <<" val  = " << val  ;
    }
}
float L1Skeleton_::getDensity(int x,int y){
    std::vector<IndexVal<float> > p = findNearestKInputPoints(x,y,1);
    if (p.size()==0){
        return -1;
    }

    return density_.at<float>(p[0].index,0);
}

float L1Skeleton_::getDensity(int x,int y, float r){
    std::vector<IndexVal<float> > inputPts = findNearestKInputPoints(x,y,40);
    std::sort(inputPts.begin(),inputPts.end());
    //if (p.size()==0){
    //    return -1;
    // }
    cv::Point2f c(x,y);
    float density=-1;
    for (size_t i=0;i< inputPts.size();++i){
        if (sqrt(inputPts[i].val) <r){
            float d = density_.at<float>(inputPts[i].index,0);
            if (d>0){
                density = std::max( density, d );
               // break;
            }
        }
    }
    return density ;

}
float L1Skeleton_::getAdaptiveH(float x,float y){
    float d = getDensity((int)x, (int)y);
    if (d <0) qDebug() <<"Error: density is undefined for p "
                       << x<<","<<y;
    //	qDebug()<<"Density at " << x <<"," << y  <<" is "<< d<< " h= " <<  h0_* pow(d,3);
    //return h0_* pow(1+d,3);
    return h0_;//*pow(1+d,2);

}

float L1Skeleton_::theta(float d, float h){
    return exp(-d*d/ pow(h*0.5f,2));
}
float L1Skeleton_::getSkeletonLinearity(int i){
	if ( i <0 || i>= sigma_.rows){
		return 0;
	}     
    return computeSkeletonAvgSigma(i,5);
}

std::vector<IndexVal<float> > L1Skeleton_::findNearestKInputPoints(float px,
                                                                  float py,
                                                                  int k)
{
    if (k > inputSize_){
        qDebug() <<  "Can not find " << k << " points from "
                 << inputSize_ << "  input points" ;
        return std::vector<IndexVal<float> > ();
    }

    std::vector<IndexVal<float> > result(k);
    cv::Mat query = (cv::Mat_<float>(1,2) << px,py);
    cv::Mat indices(1,k,CV_32S),dist(1,k,CV_32F);
    cvflann::SearchParams searchParams(64);
    inputKNNIndex_->knnSearch(query,indices,dist,k,searchParams);
    for (int i=0;i<k;++i){
        result[i] = IndexVal<float>{indices.at<int>(0,i), dist.at<float>(0,i)};
    }
    return result;

}


L1Skeleton_::SearchResult L1Skeleton_::radiusSearch(
         cv::flann::GenericIndex<cvflann::L2<float> >  *index,
         cv::Mat features, cv::Mat queries,
         cv::Mat radii){
    // initialize search result
     SearchResult result;
    int maxNeighbors = std::min(Parameters::getInstance()->skeleton_max_neighbors, index->size());
    result.indices = cv::Mat(queries.rows, maxNeighbors,CV_32S,cv::Scalar(-1));
    result.distances = cv::Mat(queries.rows, maxNeighbors,CV_32FC1,cv::Scalar(-1));

    cvflann::SearchParams searchParams(64);
    for (int r =0;r<queries.rows;++r){
        cv::Mat p(1, queries.cols, CV_32FC1, queries.ptr<float>(r)),
            indices(1, result.indices.cols, CV_32SC1, result.indices.ptr<int>(r)),
            dists(1, result.distances.cols, CV_32FC1, result.distances.ptr<int>(r));
        // radius search
        index->radiusSearch(p,
                            indices,
                            dists,
                            pow(radii.at<float>(r,0),2)/* * radii.at<float>(r,0)*/,
                            searchParams);

    }
    //qDebug() << "result size: index("<< result.indices.rows <<"x" <<result.indices.cols <<") dist("
    //            << result.distances.rows <<"x" << result.distances.cols <<")";
    return result;

}
std::vector<IndexVal<float> > L1Skeleton_::findNearestKSkeletonPoints(float px,
                                                                     float py,
                                                                     int k)
{
    if (k > X_.size()){
        qDebug() <<  "Can not find " << k << " points from "
                 << X_.size() << "  skeleton points" ;
        return std::vector<IndexVal<float> > ();
    }
    std::vector<IndexVal<float> > result(k);

    cv::Mat features(X_.size(), 2, CV_32F);

    for (int i=0; i<X_.size();++i){
        features.at<float>(i,0) = X_[i].pos.x;
        features.at<float>(i,1) = X_[i].pos.y;
    }
    cvflann::KDTreeIndexParams indexParams;
    cv::flann::GenericIndex<cv::flann::L2<float> > index(features,
                                                         indexParams);
    cv::Mat query = (cv::Mat_<float>(1, 2) << px, py );
    cv::Mat indices(1,k,CV_32S), dist(1,k,CV_32F);
    cvflann::SearchParams searchParams(64);

    index.knnSearch(query,indices,dist,k,searchParams);

    for (int i=0;i<k;++i){
        result[i] = IndexVal<float>{indices.at<int>(0,i), dist.at<float>(0,i)};
    }

    return result;
}

void L1Skeleton_::reinitialize(bool incrementH0){
    if (incrementH0){
		qDebug() <<"increment h0 by "<< Parameters::getInstance()->skeleton_h0/float(5);
        h0_+=Parameters::getInstance()->skeleton_h0/float(4); // 20. 30. 45. 65
		if (sigma_thresh_>=0.5){
			sigma_thresh_ -= 0.04;
		}
    }
    qDebug()<<"reinitialize adaptive L1Skeleton with h0 = " << h0_
			<<", sigma_thresh = " << sigma_thresh_ << endl;
}
bool L1Skeleton_::skeletonConverged(std::vector<SkeletonPoint> newX){

    //cv::Mat l2Error(X_.size(),1,CV_32F);
    float epsilon = 1e-4f;//0.1f;
    bool converged = true;
    for (size_t i=0;i<X_.size();++i){
        float err  = (newX[i].pos-X_[i].pos).dot(newX[i].pos-X_[i].pos);
        //l2Error.at<float>(i,0) = err;
        if (  err > epsilon ){
            converged = false;
            break;
        }
    }
    return converged;
}
std::vector<Candidate>::iterator
L1Skeleton_::searchCandidates(std::vector<Candidate> &candidates,int i){

    return find_if(candidates.begin(),candidates.end(),[&i](const Candidate& c){
            return c.id == i;});
}




bool L1Skeleton_::regularize(int max_iter){
    qDebug() << "Regularize with h0 = " << h0_ <<") and adaptive scale ...";
	int maxIter= max_iter;
	if (maxIter <= 0){
		maxIter =Parameters::getInstance()->skeleton_max_iter;
	}
    bool converged = false; int iter =0;

    while ( iter < maxIter && !converged ){
        converged = iterate_step();
		//	iterationPostProcess(X_,converged);
        iter++;
    }
    qDebug() <<"# of iterations: " << iter << endl;
    return converged;

}
void L1Skeleton_::removeFromCandidates(std::vector<Candidate> &candidates,
                                      std::vector<Candidate> branch){
	int ct  = candidates.size();
	std::sort(candidates.begin(),candidates.end(), compareId());
	std::sort(branch.begin(), branch.end(), compareId());
	std::vector<Candidate>::iterator itA,itB;
	for (itA = candidates.begin(), itB = branch.begin(); 
		 itB < branch.end() && itA <  candidates.end(); )
	{
		if(itA->id < itB->id) 
			++itA; 
		else if(itA->id == itB->id) {
			//increment iterator, erase the pointer before increment
			candidates.erase(itA++);			
		}else
			++itB;
	}

	
}
void L1Skeleton_::visualizeSigma(){
    //int resolution = Parameters::getInstance()->skeleton_sample_cell_size;
    //int height= ceil(voting_image_.rows/resolution);
    //int width = ceil(voting_image_.cols/resolution);
	int height= voting_image_.rows;
	int width = voting_image_.cols;
    double minV,maxV;
    cv::minMaxLoc(sigma_,&minV, &maxV);
    QImage sigmaImg1(width, height,QImage::Format_RGB32);
    QImage sigmaImg2(width, height,QImage::Format_RGB32);
    sigmaImg1.fill(QColor(Qt::white).rgb());
    sigmaImg2.fill(QColor(Qt::white).rgb());
    for (size_t i=0;i<X_.size();++i){
        float val1 = 1 - sigma_.at<float>(i,0);
        float val2 = 1 - getSkeletonLinearity(i);
        osg::Vec4 color1 = ColorMap::getInstance().getContinusColor(val1,0,1,false);
        sigmaImg1.setPixel(X_[i].pos.x, X_[i].pos.y,
                          qRgb( color1[0]*255,color1[1]*255,color1[2]*255));
        osg::Vec4 color2 = ColorMap::getInstance().getContinusColor(val2,0,1,false);
        sigmaImg2.setPixel(X_[i].pos.x, X_[i].pos.y,
                          qRgb( color2[0]*255,color2[1]*255,color2[2]*255));
    }
    double meanV = cv::mean(sigma_)[0];

    QWidget *window = new QWidget;
    QVBoxLayout *layout = new QVBoxLayout;

    QLabel *header = new QLabel(window);
    header->setText(QString("sigma min=%1 max=%2 mean=%3").arg(minV).arg(maxV).arg(meanV));
    //header->maximumHeight(20);

    QLabel *l1 = new QLabel(window);
    QLabel *l2 = new QLabel(window);
    l1->setScaledContents(true);
    l2->setScaledContents(true);

    l1->setPixmap(QPixmap::fromImage(sigmaImg1));
    l2->setPixmap(QPixmap::fromImage(sigmaImg2));

    layout->addWidget(header);
    layout->addWidget(l1);
    layout->addWidget(l2);
    window->setLayout(layout);
    window->show();

}


void L1Skeleton_::computeBranchCandidate(std::vector<Candidate> &branchCandidates,
                                        bool excludeBranchPoints){
	//    float sigmaThresh = sigma_thresh_; 
    int k = 5;
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
		if (avgSigma > sigma_thresh_){
			branchCandidates.push_back(Candidate(X_[i], getPrincipalDirection(i), 
												 avgSigma, Common::inf, false));
		
		}
	} 
	qDebug() << "total skeleton points: " << X_.size() <<" # non-branch pts: "
			 << nB << " # candidates: " << branchCandidates.size() << endl;
}

void L1Skeleton_::identifyBranchPoints(){
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
    bool converged = false;
    iterationPostProcess(X_,converged);
    qDebug() <<"finsihed processing all branches";
}
void L1Skeleton_::computeBridgePoints(int existingNumBranches){
    //std::vector<int> headBridgePts, tailBridgePts;

	if (X_.size() >= 5){
		for (int b= existingNumBranches;b<branches_.size();++b){
            markBridgePoint(b,BranchEndpointType::FIRST); // find bridge point of first endpoint
            //if (bid >=0) headBridgePts.push_back(bid);
            markBridgePoint(b,BranchEndpointType::LAST); // find bridge point of second endpoint
            //if (bid >=0) tailBridgePts.push_back(bid);
		}
	    processBridgePoints();
	}

}

//==================================================================
// definition for L1Skeleton
//==================================================================

L1Skeleton::L1Skeleton(const cv::Mat &voting_image,
					   const std::vector<cv::Mat> &heading_images):
	L1Skeleton_(voting_image,heading_images){
    heading_images_ = std::vector<cv::Mat>(heading_images.begin(),heading_images.end());
	srand (time(NULL));
	double min, max;
	cv::minMaxLoc(voting_image_, &min, &max);
	qDebug() << "image cols = " << voting_image_.cols << 
        ", ros = " << voting_image_.rows << endl;
	if (min == max) {
		qDebug()<<"voting image empty.. cannot extract skeleton\n";
		return;
	}else{
		qDebug()<<"read voting image. max val:"<<max <<" min val:"<<min<< endl;
    }
}
/*
L1Skeleton::~L1Skeleton(void){
    L1Skeleton_::~L1Skeleton_();
    //delete inputKNNIndex_;
}*/

void L1Skeleton::initialize(float vmin){
	qDebug() <<"Call L1Skeleton::initialize for vmin="<< vmin <<endl;
	float skeleton_thresh = vmin;
	if (vmin <0){
		skeleton_thresh =Parameters::getInstance()->skeleton_thresh ;
		qDebug() << "skeleton_thresh = " << skeleton_thresh;

	}
    L1Skeleton_::initialize(skeleton_thresh);
	qDebug() <<"Note: Call initialize in derived class !" << endl;
    cv::Mat skeleton_radii(X_.size(), 1, CV_32FC1); //radii for query points in X
    for (int i=0; i<X_.size();++i){
        skeleton_radii.at<float>(i,0) = getAdaptiveH(X_[i].pos.x, X_[i].pos.y);
    }
    cvflann::KDTreeIndexParams indexParams( 4 );
    cv::flann::GenericIndex<cvflann::L2<float> > skeletonIndex(skeleton_features_,
                                                               indexParams);
    skeleton_result_ = directionalRadiusSearch(&skeletonIndex,
                                                   skeleton_features_,
                                                   skeleton_features_,
                                                   skeleton_radii);
    covariances_.resize(X_.size());
    qDebug()<<"Initialize "<< covariances_.size()
			<<" covariance matrices from ("
			<< skeleton_result_.indices.rows<<","
			<< skeleton_result_.indices.cols << ") skeleton result \n" ;

    computeSigma(skeleton_features_, skeleton_result_);
    computeBeta(skeleton_features_, skeleton_result_);
	iter_++;
}

cv::Mat L1Skeleton::computeCovMatrix(cv::Mat features, cv::Mat indices,
                                             cv::Mat dist, int i){

    cv::Mat C(2,2,CV_32F,cv::Scalar(0));

    cv::Point2f xi= X_[i].pos;
    float h = getAdaptiveH(xi.x,xi.y);
    for (int ip=0; ip< indices.cols; ++ip){

        cv::Point2f xip(features.at<float>(indices.at<int>(0,ip),0),
                        features.at<float>(indices.at<int>(0,ip),1));
        if ( dist.at<float>(0,ip)>=0 ){

            float d = sqrt(dist.at<float>(0,ip));
            if (d >3*h){
                break;
            }

            cv::Mat v = (cv::Mat_<float>(1,2) << (xi-xip).x,
                         (xi-xip).y);

            cv::Mat cov(v.t()*v);
            C+= theta(d,h)*cov;
        }
    }
    return C;
}
void L1Skeleton::computeAlpha(cv::Mat features,
                                      L1Skeleton::SearchResult result){
    // alpha(i,j) where xi is a sample point in X, qj is a input point in Q
    cv::Mat &dist = result.distances;
    cv::Mat &indices = result.indices;

    alpha_ = cv::Mat(dist.rows, dist.cols, CV_32FC1,  cv::Scalar(0));
    for (int i=0; i<dist.rows;++i){
        //	int val1 = voting_image_.at<uchar>(X_[i].pos.y, X_[i].pos.x);
        //	float intensity1 = float(256-val1)/256.f;
        float h = getAdaptiveH(X_[i].pos.x, X_[i].pos.y);
        cv::Point2f xi = X_[i].pos;
        cv::Mat C = computeCovMatrix(features,indices.row(i),
                                        dist.row(i),i);
        /*qDebug()<<"Alpha Cov matrix"<< C.at<float>(0,0)<<","
                    << C.at<float>(0,1)<<","
                    << C.at<float>(1,0)<<","
                    <<C.at<float>(1,1);*/
        for (int j=0; j<dist.cols;++j){
            if (dist.at<float>(i,j) >= 0){
                float d = sqrt(dist.at<float>(i,j));

                int idx = indices.at<int>(i,j);
                cv::Point2f xj(features_.at<float>(idx,0),features_.at<float>(idx,1));

                if (d>1e7){
                    alpha_.at<float>(i,j) = 0.f;
                }else{

                    if (Parameters::getInstance()->skeleton_mvn_kernel){

                      alpha_.at<float>(i,j) = thetaMVN(xj-xi,h,C)/(d + 1e-4);
                   }else{
                    alpha_.at<float>(i,j) = theta(d,h)/(d +
                        1e-4);
                  }
                }
                //	qDebug()	<< "alpha: " << alpha_.at<float>(i,j);

            }

        }
    }
 }


void L1Skeleton::computeBeta(cv::Mat features, SearchResult result){

    cv::Mat &dist = result.distances;
    cv::Mat &indices = result.indices;
    assert(dist.cols > 1);
    beta_ = cv::Mat(dist.rows, dist.cols-1, CV_32FC1, cv::Scalar(0) );
    for (int i=0; i<dist.rows;++i){
        //	int val1 = voting_image_.at<uchar>(X_[i].pos.y, X_[i].pos.x);
        ///	float intensity1 = float(256-val1)/256.f;
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
                    if (Parameters::getInstance()->skeleton_mvn_kernel){
                        beta_.at<float>(i,j-1) = thetaMVN(xj-xi,h,covariances_[i])/(d*d + 1e-4);
                    }else{
                        beta_.at<float>(i,j-1) = theta(d,h)/(d*d + 1e-4);
                    }
                }
            }

        }

    }
    return;
}
void L1Skeleton::processBridgePoints(){
	// merge bridge points that are close to each other
	 
    std::map<int, std::set<int> >::iterator iter;

    //float res =Parameters::getInstance()->skeleton_max_trace_step;
    float r =  2*getTraceStep();
    //float r =  2*res;
    //(h0_/Parameters::getInstance()->skeleton_h0)*res;
	//h0_;//Parameters::getInstance()->skeleton_sample_cell_size*2.5;

	//    qDebug()<<"attempint to merge nearby bridge points (total: " << branchGraph_.size() << ").";
	// iterate over all bridge point
	for (iter=branchGraph_.begin();iter!=branchGraph_.end();++iter){
		// find bridge point in X_. if not found,continue;
		int bridgeId = iter->first;
		auto bridgePt = find_if(X_.begin(),X_.end(),[&bridgeId]
								(const SkeletonPoint &v){
									return v.id == bridgeId;});
		if (bridgePt== X_.end())
			continue;
		// step 1: find neighbors of bridge point!!@@@@
        SearchResult result = findSkeletonNeighbors( bridgePt->pos );
		
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

        if (bridgePts.size()>0){
			bridgePt->pos = avg*float(1.f/bridgePts.size());
			// merge all nearby bridge points 
            for (unsigned int b=0;b<bridgePts.size();++b){
                std::set<int> branchIds1 = iter->second;//branches that contains bridgePt

				// remove bridge point, merge its adacent list 
                std::map<int, std::set<int> >::iterator iter2
					=branchGraph_.find(bridgePts[b].id);

				//      qDebug() <<  "query branch graph by key " << bridgePts[b].id << " -- "<< iter2->first;
//                             <<"," << Common::printIntSet(iter2->second).c_str();

                if (iter2==branchGraph_.end()){
                    qDebug() << "skip unfound branch point " << iter2->first;
                    continue;
                }
                std::set<int> branchIds2 =iter2->second;

				// qDebug() <<"merging branches that contain neighbor bridge point "<< iter2->first
                //         << "["<< Common::printIntSet(branchIds1).c_str()<<"]";
                //qDebug() <<	" with branches that contain "<< iter->first ;
                iter->second.insert(branchIds2.begin(),branchIds2.end());

                iter2->second.clear();
                branchGraph_.erase(iter2);
				// qDebug() << "["<<	Common::printIntSet(branchIds2).c_str()<<"]";
				// qDebug()<<"result: "<< "["<<Common::printIntSet(iter->second).c_str()<<"]";
				
				
			}
		}
		
	}
}

bool L1Skeleton::iterate_step(){
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
    SearchResult result1 = directionalRadiusSearch(inputKNNIndex_,
                                                   features_,
                                                   queries,
                                                   radii);
    computeAlpha(features_,result1);

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
                xi1 += q*alpha_.at<float>(i,j)*pow(density_.at<float>(j),1.0);
                sum1 += alpha_.at<float>(i,j)*pow(density_.at<float>(j),1.0);
            }
        }
		//!!!---!!!! NEXT LINE CAUSE SEG FAULT
		if (result2.indices.cols <= 1){
			continue;
		}
        cv::Mat indices = result2.indices.colRange(1,result2.indices.cols);
        for (int ip=0; ip< indices.cols;++ip){
            int idx = indices.at<int>(i,ip);
            if (idx >=0 && X_[idx].type != SkeletonPointType::NOISE){

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

    bool converged = false; //skeletonConverged(newX);
    iterationPostProcess(newX,converged);

	iter_++;
    return converged;

}
void L1Skeleton::iterationPostProcess(const std::vector<SkeletonPoint> &newX, bool &converged){
    X_=mergeSkeletonPoints(newX);

    // recompute beta and sigma
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
    skeleton_result_ = directionalRadiusSearch(&index3,
                                                   newQueryFeature,
                                                   newQueryFeature,
                                                   newRadii);
    if (skeleton_result_.indices.cols <=1){
        converged = true;
    }else{
        computeSigma(newQueryFeature,skeleton_result_);
        computeBeta(newQueryFeature, skeleton_result_);
    }
}

void L1Skeleton::computeSigma(cv::Mat features, SearchResult result){
    assert(result.indices.cols > 1);
    cv::Mat indices  = result.indices.colRange(1, result.indices.cols);
    cv::Mat dist = result.distances.colRange(1,result.distances.cols);

    sigma_ = cv::Mat(indices.rows, 1, CV_32F, cv::Scalar(0));
    pca_directions_ = cv::Mat(indices.rows,2,CV_32F,cv::Scalar(0));

    for (int i=0; i<indices.rows;++i){
        if (!Parameters::getInstance()->skeleton_position_pca){
            qDebug() << "Error: L1Skeleton::computeSigma requires skeleton_position_pca = true.";
        }
        int k = findStableK(indices.row(i),dist.row(i),i, 3, 15, 3);
        if (k>0){
            std::pair<float,cv::Mat>  sigmaP = sigma(indices.row(i), dist.row(i),i,k);
            sigma_.at<float>(i,0) = sigmaP.first;
            (sigmaP.second).copyTo(pca_directions_.row(i));
           // qDebug() <<"select stable k "<<k;
        } else{
            qDebug() <<"stable k doesn't exists";
        }
    }
}


std::pair<float, cv::Mat> L1Skeleton::sigma(cv::Mat indices,cv::Mat dist,
                                            int i,int k){

    cv::Mat C (2,2,CV_32F,cv::Scalar(0));
    cv::Point2f xi= X_[i].pos;
    float h = getAdaptiveH(xi.x,xi.y);
    for (int ip=0; ip< std::min(k,indices.cols); ++ip){
            //			qDebug() <<i<<","<< ip <<" - " << dist.at<float>(i,ip);
        if ( dist.at<float>(0,ip)>=0 ){
            float d = sqrt(dist.at<float>(0,ip));
            int idx = indices.at<int>(0,ip);
            cv::Point2f xip= X_[idx].pos;
            cv::Mat v = (cv::Mat_<float>(1,2) << (xi-xip).x,
                         (xi-xip).y);

            cv::Mat cov(v.t()*v);
            C+= theta(d,h)*cov;
        }
    }
    if (covariances_.size()<=i){
        qDebug()<<"Error: can not find covariance matrix " << i ;
		cv::Mat cov = (cv::Mat_<float>(2,2) << 0 , 0,0,0);
		return std::make_pair(0,cov); 
	}
    C.copyTo(covariances_[i]);
    /*qDebug()<<"Beta Cov matrix"<< C.at<float>(0,0)<<","
                << C.at<float>(0,1)<<","
                << C.at<float>(1,0)<<","
                <<C.at<float>(1,1);*/
    cv::Mat eigenvalues,eigenvectors;
    eigen(C, eigenvalues,eigenvectors);

    float lambda1 = eigenvalues.at<float>(0,0),
        lambda2 = eigenvalues.at<float>(1,0);

    float sigma;
    if (lambda1+lambda2 < 1e-5){
        sigma = 0;
        //	eigenVeceigenvectors.row(0)
    }else{
        sigma= std::max(lambda1,lambda2)/(lambda1 + lambda2);
    }
    if (lambda1>lambda2){
        return std::make_pair(sigma,eigenvectors.row(0));
    }else{
        return std::make_pair(sigma,eigenvectors.row(1));
    }
}
int L1Skeleton::findStableK(cv::Mat neighbor_ids, cv::Mat dists,int point_id,
                                    int k_min, int k_max, int width){
    float maxSigma = 0;
    float optKl =k_min;
    int kl, ku, k;

    for ( kl = k_min; kl < k_max-width+1; ++kl){
        bool stable = true; float stableSigma = 0.f;
        for ( k=kl;k< kl+width; ++k){
            stableSigma = sigma(neighbor_ids, dists, point_id, k).first;
            if (k <0.85){
                stable =false;
                break;
            }
        }
        if (stable && stableSigma > maxSigma){
            maxSigma = stableSigma;

            optKl = kl;//floor(kl+width/2);
        }
    }

   // return (optKl == k_min || optKl+width == k_max-1)?0:floor(optKl+width/2);
   return floor(optKl+width/2);
    //return 6;
}



cv::Point2f L1Skeleton::addNoise(cv::Point2f p, float maxError){
	// \deprecated this function is not used! 
	int maxErrorScaled = 2*ceil(maxError*1e4);

	cv::Point2f e ((rand()% maxErrorScaled+1)/1e4,
				   (rand()%maxErrorScaled+1)/1e4);
	e.x -= maxError;
	e.y -= maxError;
	qDebug()<<"noise: " <<maxError <<" scaled: " << maxErrorScaled <<" v: "
			<< e.x <<","<<e.y <<endl;
	return p+e;
}



L1Skeleton_::SearchResult L1Skeleton::findSkeletonNeighbors(cv::Point2f p){
    //,
    //													   bool directional){
    //SearchResult result;
    int n = X_.size();
    cv::Mat features(n, 2, CV_32F);
    for (int i=0; i<n;++i){
        features.at<float>(i,0) = X_[i].pos.x;
        features.at<float>(i,1) = X_[i].pos.y;
    }
    float radius = getAdaptiveH(p.x,p.y);

    cvflann::KDTreeIndexParams indexParams;
    cv::flann::GenericIndex<cv::flann::L2<float> > index(features,
                                                         indexParams);
    cv::Mat query = (cv::Mat_<float>(1, 2) << p.x, p.y );
    cv::Mat radii = (cv::Mat_<float>(1,1) << radius);
    //return directionalRadiusSearch(&index, features, query,radii,
    //                               directional);
    return radiusSearch(&index, features, query,radii);
}


L1Skeleton_::SearchResult L1Skeleton::directionalRadiusSearch(
         cv::flann::GenericIndex<cvflann::L2<float> >  *index,
         cv::Mat features, cv::Mat queries,
         cv::Mat radii, bool directional){
    // initialize search result

    SearchResult result =  radiusSearch(index,features,queries,radii);

    if (directional){
        for (int r =0;r<queries.rows;++r){
            cv::Mat indices(1, result.indices.cols, CV_32SC1, result.indices.ptr<int>(r)),
                    dists(1, result.distances.cols, CV_32FC1, result.distances.ptr<int>(r));
            int qx = queries.at<float>(r,0), //query pixel col
                qy = queries.at<float>(r,1); // query pixel row
            // filter neighbors not in the same general direction
            filterSearchResultByDirection(features,
                                          getHeadingDistribution(qx,qy),
                                          indices,dists);
        }
    }
    return result;

}


float L1Skeleton::thetaMVN(cv::Point2f x, float h,cv::Mat C){
    cv::Mat Xt = (cv::Mat_<float>(1,2) << x.x,x.y);
    cv::Mat X;
    cv::transpose(Xt,X);// = (cv::Mat<float>(2,1) << x.x ,x.y );
    cv::Mat product = -Xt*(C.inv()*X);
    return exp(product.at<float>(0)/(2*sqrt(h)));
}


float L1Skeleton_::thetaMVNRotated(cv::Point2f x,float h,cv::Mat C){
    cv::Mat X = (cv::Mat_<float>(2,1) << x.x, x.y );
    cv::Mat R = (cv::Mat_<float>(2,2) << 0,-1, 1, 0);
    X = R*X;

    cv::Mat Xt;
    cv::transpose(X,Xt);

    cv::Mat product = -Xt*(C.inv()*X);
    return exp(product.at<float>(0)/(2*sqrt(h)));
}

void L1Skeleton::visualizeDensity(){
    int height= voting_image_.rows;
    int width = voting_image_.cols;
    double minV,maxV;
    cv::minMaxLoc(density_,&minV, &maxV);

    QImage sigmaImg1(width, height,QImage::Format_RGB32);
     sigmaImg1.fill(QColor(Qt::white).rgb());
     for (size_t i=0;i<X_.size();++i){
         float val1 = 1 - density_.at<float>(i,0);
         osg::Vec4 color1 = ColorMap::getInstance().getContinusColor(val1,minV,maxV,false);
         sigmaImg1.setPixel(X_[i].pos.x, X_[i].pos.y,
                           qRgb( color1[0]*255,color1[1]*255,color1[2]*255));
    }
    double meanV = cv::mean(density_)[0];

    QWidget *window = new QWidget;
    QVBoxLayout *layout = new QVBoxLayout;

    QLabel *header = new QLabel(window);
    header->setText(QString("Density min=%1 max=%2 mean=%3").arg(minV).arg(maxV).arg(meanV));
    QLabel *l1 = new QLabel(window);
    l1->setScaledContents(true);
    l1->setPixmap(QPixmap::fromImage(sigmaImg1));

    layout->addWidget(header);
    layout->addWidget(l1);
    window->setLayout(layout);
    window->show();

}



/*
//uniform sampling while computing average intensity

void L1Skeleton::sampleUniform(){
	int resolution = Parameters::getInstance()->skeleton_sample_cell_size;
	int maxN = Parameters::getInstance()->skeleton_n_samples_per_cell;
	
	int height= ceil(voting_image_.rows/resolution);
	int width = ceil(voting_image_.cols/resolution);
	int index = 0;

	//	float minIntensity = std::numeric_limits<float>::max();
	//float maxIntensity = std::numeric_limits<float>::min();
	density_ = cv::Mat(inputSize_,1,CV_32FC1,cv::Scalar(1));
	for (int r=0;r<height;++r){
		for (int c=0; c<width;++c){

			std::vector<cv::Point2f> pts ;
			float avgIntensity = 0;
			
			for (int rr=r*resolution;rr<(r+1)*resolution;++rr){
				for (int cc=c*resolution;cc< (c+1)*resolution;++cc){
					if (Q_.at<uchar>(rr,cc)==1){
						pts.push_back(cv::Point2f(cc,rr));
						
						cv::Vec3b b = voting_image_.at<uchar>(rr,cc);
						avgIntensity += b[0];
					}
				}
			}
			if (pts.size()>0){

				avgIntensity /= resolution*resolution;
				density_.at<float>(index) = avgIntensity;
				minIntensity = std::min(minIntensity, avgIntensity);
				maxIntensity = std::max(maxIntensity, avgIntensity);
				qDebug() << "computed intensity for point " << index << 
					" in Q " << avgIntensity ;

				I_.push_back(SkeletonPoint{index++,
							cv::Point((c+0.5)*resolution, (r+0.5)*resolution), 
							SkeletonPointType::NONE});

			}
		}
	}

	density_ = (density_- minIntensity)*(1.f / (maxIntensity-minIntensity))+1;

	qDebug() <<"Total # of cells: " << height*width << endl;
	qDebug() <<"Sampled " << I_.size() <<" points out of "
			 << inputSize_+I_.size()
			 << " (%" <<float(I_.size())/ float(inputSize_+I_.size()) *100 
			 <<")" << " max index = " << index <<  endl;
	X_ = I_;
	return;
}
*/
void L1Skeleton_::sample(){
	int resolution = Parameters::getInstance()->skeleton_sample_cell_size;
	int maxN = Parameters::getInstance()->skeleton_n_samples_per_cell;
	
	int height= ceil(voting_image_.rows/resolution);
	int width = ceil(voting_image_.cols/resolution);
	int index = 0;

	for (int r=0;r<height;++r){
		for (int c=0; c<width;++c){

			std::vector<cv::Point2f> pts ;
			double avgIntensity = 0;
			
			for (int rr=r*resolution;rr<(r+1)*resolution;++rr){
				for (int cc=c*resolution;cc< (c+1)*resolution;++cc){
					if (Q_.at<uchar>(rr,cc)==1){
						pts.push_back(cv::Point2f(cc,rr));
						
						cv::Vec3b b = voting_image_.at<uchar>(rr,cc);
						avgIntensity += b[0];
					}
				}
			}
			if (pts.size()>0){

				avgIntensity /= pts.size();

                int n= ceil( float(256-avgIntensity)/256.f * maxN);
				std::vector<cv::Point2f>::iterator last = pts.end();
				if (pts.size() > n )	{		
					// random shuffle point list
					int max  = pts.size()-1;
					while(max !=0){
						int re = rand()%(max);
						cv::Point temp = pts[re];
						pts[re] = pts[max];
						pts[max] = temp;
						max--;
			
					}
					// collect first n points from shuffled list 
					last = pts.begin()+n;
				}
				for (std::vector<cv::Point2f>::iterator it=pts.begin();it!=last;
					 ++it){
					I_.push_back(SkeletonPoint{index,*it,SkeletonPointType::NONE});
					index ++;
				}
			}
		}
	}
	qDebug() <<"Total # of cells: " << height*width << endl;
	qDebug() <<"Sampled " << I_.size() <<" points out of "
			 << inputSize_+I_.size()
			 << " (%" <<float(I_.size())/ float(inputSize_+I_.size()) *100 
			 <<")" << " max index = " << index <<  endl;
	X_ = I_;
	return;
}

std::vector<SkeletonPoint> L1Skeleton::mergeSkeletonPoints(std::vector<SkeletonPoint> skeletonPts){
    //grid filter on skeleton,grid cell size is the original sample_cell_size/2
	std::map<std::pair<int,int>, std::vector<int> > gridFilter;
    float r = Parameters::getInstance()->skeleton_sample_cell_size *0.5f;
    //loop through all non-noise points in X, place them into gridFilter
	//loop through each pair in the map, if size > 2, modify the first element
	//add to new skeleton point list
	for (int i=0; i<skeletonPts.size();++i){
        if (skeletonPts[i].type != SkeletonPointType::NOISE){
            std::pair<int,int> key;
            key.first = floor(skeletonPts[i].pos.x/r);
            key.second = floor(skeletonPts[i].pos.y/r);
            gridFilter[key].push_back(i);
        }
	}
	std::vector<SkeletonPoint> result;
	for (auto it=gridFilter.begin();it!=gridFilter.end();++it){
		int groupSize = it->second.size();
		if (groupSize ==1){
			result.push_back(skeletonPts[it->second[0]]);
		}else if (groupSize >1){
			// compute the average position of t
			// push the first into result, modify its position 
			cv::Point2f avg(0,0);
			bool hasBranchPts = false, hasBridgePts = false;
			int bridgeId=-1, branchId=-1;
			for (int j=0; j < groupSize;++j){
				SkeletonPoint pt = skeletonPts[it->second[j]];
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

			result.push_back(s);
		}
		

	}
	qDebug() << "before: " << skeletonPts.size() << " points, after: " <<
		result.size() << " points.";
	return result;
}



double L1Skeleton::evalEnergyFunction(SearchResult inputSearch, 
									  SearchResult sampleSearch, 
									  float h){
	double totalEnergy=0.f;
	
	cv::Mat alphaSum, betaSum;
	reduce(alpha_, alphaSum , 1, CV_REDUCE_SUM);
	reduce(beta_,betaSum, 1, CV_REDUCE_SUM);

	for (int i=0; i< X_.size();++i){

		for (int j=0; j<inputSearch.distances.cols;++j){
			float d = inputSearch.distances.at<float>(i,j);
			if (d>=0){
				d = sqrt(d);
				totalEnergy += d*theta(d,h);
			}
		}


		cv::Mat dist
			= sampleSearch.distances.colRange(1,sampleSearch.distances.cols);

		float sigma_i = sigma_.at<float>(i,0);
		float mu = Parameters::getInstance()->skeleton_mu;
		float gamma  = mu*sigma_i*alphaSum.at<float>(i,0)/betaSum.at<float>(i,0);

		for (int ip=0; ip< dist.cols;++ip){
			float d = dist.at<float>(i,ip+1e-5);
			if (d> 0 && sigma_i > 0 ){
				d = sqrt(d);
				totalEnergy += gamma * theta(d,h)/(sigma_i*d);
			}
		}

	}
	return totalEnergy;
}

cv::Mat L1Skeleton::getHeadingDistributionNS(int x, int y,bool recursive){

    // compute the non symmetric version of heading distribution
    int nDirections = heading_images_.size() ;
    bool allZero =true;
    cv::Mat headingDist(1,nDirections,CV_32F,cv::Scalar(0)),
            headingDistNormalized;
    for (size_t d=0; d< nDirections;++d){
        headingDist.at<float>(0,d) = heading_images_[d].at<float>(y,x);
        if (headingDist.at<float>(0,d)>1e-5)
            allZero = false;
    }
    cv::normalize(headingDist,headingDistNormalized);

    if (allZero && recursive){
        // use weighted neighbors
        auto nlist = findNearestKInputPoints(x,y,1);
        if (nlist.size()>0){
            IndexVal<float> nearest = nlist.front();
            cv::Point2f p(features_.at<float>(nearest.index,0),
                      features_.at<float>(nearest.index,1));
            headingDistNormalized = getHeadingDistributionNS(p.x, p.y, false);
        }
    }
    return headingDistNormalized;
}

cv::Mat L1Skeleton::getHeadingDistribution(int x, int y,bool recursive){
    int nBiDirections = (heading_images_.size() ==2)? 2:heading_images_.size()/2;

	cv::Mat headingDist(1,nBiDirections,CV_32F,cv::Scalar(0));
    if (nBiDirections==2){
         headingDist.at<float>(0,0) = heading_images_[0].at<float>(y,x);
         headingDist.at<float>(0,1) = heading_images_[1].at<float>(y,x);
           headingDist /=  headingDist.at<float>(0,0)+headingDist.at<float>(0,1);

         if (headingDist.at<float>(0,0)>1e-5 && headingDist.at<float>(0,1)>1e-5
                 && recursive){
             auto nlist = findNearestKInputPoints(x,y,1);
             if (nlist.size()>0){
                 IndexVal<float> nearest = nlist.front();
                 cv::Point2f p(features_.at<float>(nearest.index,0),
                           features_.at<float>(nearest.index,1));
                 headingDist = getHeadingDistribution(p.x, p.y, false);
             }
         }
    }else{
        for (size_t d=0; d<heading_images_.size();++d){
            int modD = d % nBiDirections;
            headingDist.at<float>(0,modD) += heading_images_[d].at<float>(y,x);
        }
        std::ostringstream oss;
        oss << "distrib("<< x<<","<<y<<"):";
        bool allZero =true;
        for (int i=0;i<nBiDirections;++i){
            oss <<  headingDist.at<float>(0,i)<< ", ";
            if (headingDist.at<float>(0,i)>1e-5)
                allZero = false;
        }
        if (allZero && recursive){
            // use weighted neighbors
            auto nlist = findNearestKInputPoints(x,y,1);
            if (nlist.size()>0){
                IndexVal<float> nearest = nlist.front();
                cv::Point2f p(features_.at<float>(nearest.index,0),
                          features_.at<float>(nearest.index,1));
                headingDist = getHeadingDistribution(p.x, p.y, false);
            }
        }

    }
	return headingDist;
}

/*L1Skeleton::SearchResult L1Skeleton::radiusSearch(
    cv::flann::GenericIndex<cvflann::L2<float> > *index,
	cv::Mat features,cv::Mat queries, float radius){
	// compute k nearest neighbors of each xi in Q

	SearchResult result;
	result.indices = cv::Mat(queries.rows, index->size(),CV_32S,cv::Scalar(-1));
	result.distances = cv::Mat(queries.rows, index->size(),CV_32FC1,cv::Scalar(-1)); 


	cvflann::SearchParams searchParams(128);
	for (int r =0;r<queries.rows;++r){
		cv::Mat p(1, queries.cols, CV_32FC1, queries.ptr<float>(r)),
			indices(1, index->size(), CV_32SC1, result.indices.ptr<int>(r)),
			dists(1, index->size(), CV_32FC1, result.distances.ptr<int>(r));
		// radius search
		index->radiusSearch(p,
							indices, 
							dists, 
							radius, searchParams);
	}

	return result;
	}*/
void L1Skeleton::filterSearchResultByDirection(cv::Mat features,
											   cv::Mat headingDistQ,
											   cv::Mat &indices, 
											   cv::Mat &dists){
    if (heading_images_.size()==2){
        filterSearchResultByDirection2(features,headingDistQ, indices,dists);
        return;
    }
	double min, max;
	cv::minMaxLoc(headingDistQ, &min, &max);
	if (max <1e-4){
		return;
	}
	// assume the number of discrete directions is even
	assert( heading_images_.size()%2 == 0);
	int removed =0;
	for (int i = 0; i<indices.cols;++i){
		// for each neighbor pi, find the symmetric chi square error between 
		// heading distribution of pi and query point pr. 
		// chi(p,q) = sum((p-q)^2/p+q)/2

		int neighborIdx = indices.at<int>(0,i); //neighbor index
		if (neighborIdx <0 || neighborIdx > features.rows) 
			continue;
		int px = features.at<float>(neighborIdx,0); // neighbor pixel col
		int py = features.at<float>(neighborIdx,1); // neighbor pixel row
		cv::Mat headingDistP = getHeadingDistribution(px,py);

		cv::Mat diff = headingDistQ - headingDistP;

		cv::minMaxLoc(headingDistP, &min, &max);
		if (max <1e-4){
			// skip if this point doesn't have a direction
			qDebug() <<"zero neighbor heading: "<< px<<","<< py;
			continue;
		}
		cv::Mat tmp;
		cv::multiply(diff,diff,tmp);
			
		tmp = tmp/ (headingDistQ+headingDistP);
		double chiSqr =0.5* cv::sum(tmp)[0];			 
		//			qDebug() <<"headingerror: " <<  chiSqr;
		//remove neighbors with l2 heading error larger than heading_threshold
		if (chiSqr > Parameters::getInstance()->skeleton_heading_thresh){
			indices.at<int>(0,i) = -1;
			dists.at<float>(0,i) = -1;
			removed++;
		}
			
						
	}
	//	qDebug()<<"filtered: " << removed <<" out of max of " << indices.cols;
}
void L1Skeleton::filterSearchResultByDirection2(cv::Mat features,
                                               cv::Mat headingDistQ,
                                               cv::Mat &indices,
                                               cv::Mat &dists){
    cv::Point2f vq (headingDistQ.at<float>(0), headingDistQ.at<float>(1));
    int removed =0;
    for (int i = 0; i<indices.cols;++i){
        int neighborIdx = indices.at<int>(0,i); //neighbor index
        if (neighborIdx <0 || neighborIdx > features.rows)
            continue;
        int px = features.at<float>(neighborIdx,0); // neighbor pixel col
        int py = features.at<float>(neighborIdx,1); // neighbor pixel row


        cv::Mat headingDistP = getHeadingDistributionNS(px,py);
        cv::Point2f vp(headingDistP.at<float>(0),headingDistP.at<float>(1));
       // if (cv::norm(vp)!=1)
        //    qDebug() <<"heading distr vector has non unit length " << cv::norm(vp);
        if (fabs(vq.dot(vp)) < Parameters::getInstance()->skeleton_heading_thresh){
            indices.at<int>(0,i) = -1;
            dists.at<float>(0,i) = -1;
            removed++;
        }

    }
   qDebug()<<"filtered: " << removed <<" out of max of " << indices.cols;
}
														  
/*
L1Skeleton::SearchResult L1Skeleton::directionalRadiusSearch(
	   cv::flann::GenericIndex<cvflann::L2<float> >  *index, 
	   cv::Mat features, cv::Mat queries,
	   float radius, bool directional){
	// initialize search result
	SearchResult result;
	int maxNeighbors = std::min(Parameters::getInstance()->skeleton_max_neighbors, index->size());
	result.indices = cv::Mat(queries.rows, maxNeighbors,CV_32S,cv::Scalar(-1));
	result.distances = cv::Mat(queries.rows, maxNeighbors,CV_32FC1,cv::Scalar(-1)); 

	//	result.indices = cv::Mat(queries.rows, index->size(),CV_32S,cv::Scalar(-1));
	//result.distances = cv::Mat(queries.rows, index->size(),CV_32FC1,cv::Scalar(-1)); 

	cvflann::SearchParams searchParams(64);
	for (int r =0;r<queries.rows;++r){
		cv::Mat p(1, queries.cols, CV_32FC1, queries.ptr<float>(r)),
			indices(1, result.indices.cols, CV_32SC1, result.indices.ptr<int>(r)),
			dists(1, result.distances.cols, CV_32FC1, result.distances.ptr<int>(r));
		// radius search
		index->radiusSearch(p,
							indices, 
							dists, 
							radius*radius, searchParams);

		// filter neighbors not in the same general direction 
		if (directional){
			int qx = queries.at<float>(r,0), //query pixel col
				qy = queries.at<float>(r,1); // query pixel row

			filterSearchResultByDirection(features,
										  getHeadingDistribution(qx, qy),
										  indices,dists);		
		}
	}
	return result;
}

*/


cv::Point2f L1Skeleton::getPrincipalDirection(int i){
	if (i >= 0 && i < pca_directions_.rows){
		//	if (pca_directions_.at<float>(i,0) == 0 || pca_directions_.at
		return cv::Point2f(pca_directions_.at<float>(i,0), 
						   pca_directions_.at<float>(i,1));
	}else{
		qDebug() << "can not find pca direction for point " << i ;
		return cv::Point2f(0,0);
	}
}
float L1Skeleton_::gaussian(float mean,float stdev){
    float var2 = 2*stdev*stdev;
    return 1/sqrt(M_PI*var2) * exp(-mean*mean /var2);
}
/*
double L1Skeleton_::computeSkeletonAvgSigma(int skeletonIdx, int  k){
    // compute skeleton sigma
    // [Note] skeletonIdx is the index of skeleton point in X_. ,
    // typically X_[skeletonIdx].id != skeletonIdx
     if (skeletonIdx <0 || skeletonIdx >= sigma_.rows)
        return -1;
    //neighbors = findNearestKSkeletonPoints( X_[skeletonIdx].pos.x,
    //										X_[skeletonIdx].pos.y, k);
    SearchResult result = findSkeletonNeighbors(X_[skeletonIdx].pos);
    float vote =0.f;int nvotes =0;
    for (int i=0; i<result.indices.cols;++i){
        int neighborIdx = result.indices.at<float>(0,i);
        if (neighborIdx <0)
            break;
        float sigma = sigma_.at<float>(neighborIdx ,0);
        nvotes++;
        vote += theta(result.distances.at<float>(0,i)/2.f,
                      getAdaptiveH(X_[neighborIdx].pos.x,X_[neighborIdx].pos.y))* sigma;
    }
    return vote/float(nvotes);
}*/

double L1Skeleton_::computeSkeletonAvgSigma(int skeletonIdx, int  k){
	// compute skeleton sigma
	// [Note] skeletonIdx is the index of skeleton point in X_. ,
	// typically X_[skeletonIdx].id != skeletonIdx
	std::vector<IndexVal<float> > neighbors;
	if (sigma_.rows == 0)
		return -1;
	neighbors = findNearestKSkeletonPoints( X_[skeletonIdx].pos.x, 
											X_[skeletonIdx].pos.y, k);

    //float h =
	float avgSigma =0;
	std::vector<IndexVal<float> >::iterator it;
    int ct = 0; float hScale = 2;float wSum=0.f;
	for (it=neighbors.begin();it!=neighbors.end();++it){
		if (it->val <= h0_*hScale){
            float r = (it->val>=0)?sqrt(it->val):5*h0_;
            avgSigma += theta(r,h0_) *sigma_.at<float>( it->index, 0);
            wSum+=theta(r,h0_) ;
			ct++;
		}
	}

    if ( wSum>0) avgSigma/= wSum;
	return avgSigma; 
}


std::pair<std::vector<Candidate>,std::vector<Candidate> >
L1Skeleton::getBranchFromSeed(int skeleton_id){
	// compute candidate
	// if skeleton_id is candidate, trace branch
	// [Note]: skeleton_id is the index of skeleton point in X_, 
	// generally X_[skeleton_id].id != skeleton_id!
	qDebug() <<"noise scale " << Parameters::getInstance()->skeleton_noise_scale;
	qDebug() <<"max trace step " << getTraceStep();
	
	std::vector<Candidate> branchCandidates;
	computeBranchCandidate(branchCandidates,true);//exclude branch points
	int pid = X_[skeleton_id].id;
	std::vector<Candidate>::iterator seedIter = 
		find_if(branchCandidates.begin(),branchCandidates.end(),
				[&pid](const Candidate &c){return c.id == pid;});
	if (seedIter==branchCandidates.end()){
		return std::make_pair(std::vector<Candidate>(),std::vector<Candidate>()); //return empty branch if skeleton pt is not candidate
	}
	std::vector<Candidate> noise;
	Candidate seed =*seedIter;
	std::vector<Candidate> branchForward =
		traceBranch(branchCandidates, seed,noise);

 	seedIter = 
		find_if(branchCandidates.begin(),branchCandidates.end(),
				[&pid](const Candidate &c){return c.id == pid;});

	std::vector<Candidate> 	branch = traceBranchReverse(branchCandidates,
														seed,noise);
    branch.pop_back();
    branch.insert(branch.end(),branchForward.begin(),branchForward.end());

	return std::make_pair(branch,noise);
}

void L1Skeleton_::printBranchGraph(){
	qDebug() << "======= branch graph =======";
    for (auto it=branchGraph_.begin(); it != branchGraph_.end(); ++it){
		std::ostringstream oss;
		oss << it->first <<" --> [\n";

        for (auto bit=(it->second).begin(); bit!= (it->second).end() ;++bit){
            std::vector<Candidate> branch = branches_[*bit];
            oss <<*bit <<" (" << branch[0].id;

			for (size_t j=1; j < branch.size(); ++j){
				oss << ", " <<  branch[j].id ;
			}
			oss <<") \n";
		}
		oss <<"]\n";
		qDebug() << oss.str().c_str();
	}
	//	qDebug() << " \n";
}
/*void L1Skeleton::smoothBranches(){
  }

void L1Skeleton::recenterBranches(){
	// for each branch point in a branch, project neighbors to the perpendicular (minor pca), move point to center
	for (size_t i=0; i<branches_.size(); ++i){
		std::vector<Candidate> &branch = branches_[i];	  
		float r =0.5*Parameters::getInstance()->skeleton_sample_cell_size;
	   
		for (size_t j = 0; j <branch.size();++j){
			cv::Mat queries = (cv::Mat_<float>(1,2) << branch[j].pos.x , branch[j].pos.y );
			//	branch[j].pos
			SearchResult result = directionalRadiusSearch(inputKNNIndex_, features_, queries, r,true);
			int pid = branch[j].id;
			
			auto iter = find_if(X_.begin(),X_.end(),[&pid](const SkeletonPoint &v){return v.id == pid;});
			if (iter== X_.end())
				continue;
            //			cv::Point2f dir(pca_directions_.at<float>(iter-X_.begin(),0),
                //			pca_directions_.at<float>(iter-X_.begin(),1));

			cv::Point2f dir;
			if (j>0 && j <branch.size()-1){
				dir = 0.5*(branch[j+1].pos - branch[j-1].pos);
			}
			else {
				dir.x =  pca_directions_.at<float>(iter-X_.begin(),0);
				dir.y = pca_directions_.at<float>(iter-X_.begin(),1);
			}
			cv::Point2f projAxis(-dir.y, dir.x) , projAxisN;
			Common::normalize(projAxis,projAxisN);
			//std::vector<Point2f> neighbors;
			std::vector<float> projDists;
			for (int k=0; k< result.indices.cols;++k){
				// return the resulting 
				int nId = result.indices.at<int>(0,k);
				if (nId != -1){
					cv::Point2f nPos(features_.at<float>(nId,0), 
									 features_.at<float>(nId,1));
					//	neighbors.push_back(cv::Point2f(features_.at<float>(nId,0), 
					//								featuers_.at<float>(nId,1)));
					cv::Point2f vec = nPos-branch[j].pos;
					projDists.push_back(vec.dot(projAxisN));
				}else{
					break;
				}
			}
			
			int m = floor(projDists.size()/2);
			nth_element(projDists.begin(),projDists.begin()+m, 
						projDists.end());
			if (m>=1){

				double median =projDists[m-1];
				qDebug() <<"recentered to " <<  m << " th neighbor out of " 
						 << projDists.size();
				qDebug() <<"displacement:  " <<  median ;
				branch[j].pos = branch[j].pos + projAxisN*median;
				iter->pos = branch[j].pos;
			}
			//for each result, project to vec (-y,x)
			// project results and find center
			//!!!@@@@@@@@@@@@@@

			
		}
	}
}*/

float L1Skeleton_::angleBetween(cv::Point2f v1, cv::Point2f v2){
	float len1 = sqrt(v1.x * v1.x + v1.y * v1.y);
	float len2 = sqrt(v2.x * v2.x + v2.y * v2.y);

	float dot = v1.x * v2.x + v1.y * v2.y;

	float a = dot / (len1 * len2);

	if (a >= 1.0)
		return 0.0;
	else if (a <= -1.0)
		return M_PI;
	else
		return acos(a); // 0..PI
}
void L1Skeleton_::markBridgePoint(int branchId,
								BranchEndpointType endType){
    // Given branch B and its endpoint e,
    // mark the nearest non-branch skeleton point in N_k(e) in the same
    // direction of the branch as a Bridge point
	std::vector<Candidate> &branch = branches_[branchId];	  
	std::vector<Candidate>::iterator lastIter = branch.end()-1, 
		firstIter = branch.begin();
	bool first = endType==BranchEndpointType::FIRST;
    cv::Point2f p = (first)? firstIter->pos:lastIter->pos;

    size_t k = 20;
	std::vector<IndexVal<float> > result = 	
		findNearestKSkeletonPoints(p.x, p.y, std::min(k, X_.size()-1));
	std::sort(result.begin(), result.end());
	cv::Point2f v0,v1;
    //if (branch.size()>= 2){
    v0 = (first)? (branch[1].pos - firstIter->pos):
        (branch[branch.size()-2].pos -lastIter->pos);
    //}else{
    //	v0 = (first)? (firstIter->dir) : -(lastIter->dir);
    //}
    //float branchCosThresh = Parameters::getInstance()->skeleton_branch_cos_thresh;
    std::vector<int> scannedPoints;
    bool foundBridgePt=false;
    for (std::vector<IndexVal<float> >::iterator it=result.begin();
		 it != result.end();++it){
		SkeletonPoint &x = X_[it->index];//ith neighbor
		if (x.type == SkeletonPointType::BRANCH)
			continue; 

		v1 = (first)? (x.pos - firstIter->pos): 
			(x.pos - lastIter->pos);  

        if (cos(angleBetween(v0, v1))< 0 && sqrt(it->val) < getTraceStep()){
		   
			x.type = SkeletonPointType::BRIDGE;

            foundBridgePt = true;
            if (branchGraph_.find(x.id) == branchGraph_.end() ){
				//				qDebug() <<"add new node " << x.id << " to branch graph";
                branchGraph_[x.id] = std::set<int>();
			}
            branchGraph_[x.id].insert(branchId);
            break;
        }else if (X_[it->index].type == SkeletonPointType::NONE){
            scannedPoints.push_back(it->index);
        }
	}
    /** // mark points closer to bridge points as noise
    if (foundBridgePt){
        for (size_t i=0;i<scannedPoints.size();++i){

            //    X_[i].type = SkeletonPointType::NOISE;

        }
        qDebug() <<"[skipped!] mark "<< scannedPoints.size() <<
                    " points as noise on endpoints of branch "<< branchId;

    }
    **/
}

bool L1Skeleton_::hasUnvisited(const std::vector<Candidate> &candidates){
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

std::vector<SkeletonPoint>  L1Skeleton_::getBridgePoints(){
	//!!@@ this function no longer works due to X_ index changes
	std::vector<SkeletonPoint> pts;
    for (auto it = branchGraph_.begin(); it != branchGraph_.end();++it){
		int id = it->first;
		auto bridgePt = find_if(X_.begin(), X_.end(),
								[&id](const SkeletonPoint &v){return v.id == id;});
		if (bridgePt != X_.end()){
            pts.push_back(*bridgePt);
        }
	}
	return pts;
}


Candidate L1Skeleton_::getCandidateInBranch(int branchId, int skeletonId){
	std::vector<Candidate> branch =  branches_[branchId];
	std::vector<Candidate>::iterator cit;

	cit = searchCandidates(branch,skeletonId);
	if (cit!= branch.end()){
        /*qDebug() << "Found skeleton point " << skeletonId << " in branch "
                 << branchId;*/
		return *cit;
	}else{
		return Candidate();//-1,cv::Point2f(0,0),);
	}
}

void L1Skeleton_::markVisited(std::vector<Candidate> &candidates, 
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
void L1Skeleton_::computeProjectionDistance(std::vector<Candidate> &candidates,
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


float L1Skeleton_::getTraceStep(){
	float res =Parameters::getInstance()->map_image_resolution; //skeleton_max_trace_step;
	float step =  (h0_/float(Parameters::getInstance()->skeleton_h0))*res;
	assert(step!=0);
	return step;
}
std::vector<Candidate> L1Skeleton_::traceBranchReverse(
      std::vector<Candidate> &candidates,
	  Candidate seed,
	  std::vector<Candidate> &noise){

	float maxTraceStep =getTraceStep();
	float noiseScale = Parameters::getInstance()->skeleton_noise_scale;
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
        std::vector<Candidate> noiseCandidates;
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
                if ( cv::norm(vec) < maxTraceStep*noiseScale  )
                  noiseCandidates.push_back(*it);
			}
		}
        if (!branchEnded){
            noise.insert(noise.end(),noiseCandidates.begin(),noiseCandidates.end());
        }


    }while(!branchEnded && branch.size() < max_branch_size);

	return branch;
}
std::vector<Candidate> L1Skeleton_::traceBranch(std::vector<Candidate> &candidates,
											   Candidate seed,
											   std::vector<Candidate> &noise){
	//float res =Parameters::getInstance()->skeleton_max_trace_step;

  	float maxTraceStep = getTraceStep();// (h0_/Parameters::getInstance()->skeleton_h0)*res;
	float noiseScale = Parameters::getInstance()->skeleton_noise_scale;
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
		std::vector<Candidate> neighbors(candidates.begin(),candidates.begin()+neighborSize);
		std::vector<Candidate>::iterator it;
		branchEnded = true;

         std::vector<Candidate> noiseCandidates;
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
                if ( cv::norm(vec) < float(maxTraceStep)* noiseScale )
                    noiseCandidates.push_back(*it);
			}
		}
        if (!branchEnded){
            noise.insert(noise.end(),noiseCandidates.begin(),noiseCandidates.end());
        }

    }while(!branchEnded &&  branch.size() <max_branch_size);

	return branch;
}

void L1Skeleton_::processBranch(std::vector<Candidate> &candidates,
							   std::vector<Candidate> branch,
							   std::vector<Candidate> noise){

	//retain branch if its length is at least 5
	if (branch.size() >= 4){
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
		//        qDebug() <<"trace a branch of length " << branch.size() << " with " << noise.size() << " noise points";
		markVisited(candidates,branch);
		removeFromCandidates(candidates,noise);
        removeFromCandidates(candidates,branch);
       //
        //markAsNoise
        for (int i=0; i< (int)noise.size() ; ++i){
            int id = noise[i].id;
            auto vit = find_if(X_.begin(), X_.end(), [&id](const SkeletonPoint &pt)
                               {return pt.id == id;} );
            if (vit->type != SkeletonPointType::BRANCH && vit->type != SkeletonPointType::BRIDGE){
                vit->type = SkeletonPointType::NOISE;
            }
        }
    }
}


std::pair<float, cv::Mat> L1Skeleton::sigmaV(cv::Mat indices,
                                                    cv::Mat dist, int i ){
    int dim = heading_images_.size();
    int n = 2;
    cv::Mat C(std::min(n,indices.cols),dim,CV_32F,cv::Scalar(0));
    cv::Point2f xi=X_[i].pos;
    float h = getAdaptiveH(xi.x,xi.y);

    for (int ip=0; ip< std::min(n, indices.cols) ; ++ip){
        // find velocity vectors of each neighbor and form the C function
        // cv::Mat heading =getHeadingDistribution(xi.x,xi.y);

            if ( dist.at<float>(0,ip)>=0 ){
                float d = sqrt(dist.at<float>(0,ip));
                int idx = indices.at<int>(0,ip);
                cv::Point2f xip= X_[idx].pos;
                cv::Mat v = getHeadingDistributionNS(xip.x,xip.y); //row vector
                cv::Mat cov(v.t()*v);
                C+= theta(d,h)*cov;
            }


    }
    C.copyTo(covariances_[i]);
    cv::Mat eigenvalues,eigenvectors;
    eigen(C, eigenvalues,eigenvectors);

    float lambda1 = eigenvalues.at<float>(0,0),
        lambda2 = eigenvalues.at<float>(1,0);

    float sigma;
    if (lambda1+lambda2 < 1e-5){
        sigma = 0;
    }else{
        sigma = std::max(lambda1,lambda2)/(lambda1 + lambda2);
    }
    if (lambda1>lambda2){
        return std::make_pair(sigma, eigenvectors.row(0));
    }else{
        return std::make_pair(sigma, eigenvectors.row(1));
    }
}


