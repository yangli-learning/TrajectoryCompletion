#include <sstream>
#include <iomanip>

#include <boost/throw_exception.hpp>
#include <QDebug>
#include <sstream>
#include "common.h"

namespace Common {
std::string int2String(int i, int width) {
  std::ostringstream ss;
  ss << std::setw(width) << std::setfill('0') << i;
  return ss.str();
}

float percentile(const cv::Mat &m, float p, int nBins ){
    /*cv::Mat mask(m.size(),CV_8U,cv::Scalar(1));
    float min = 0.01;
    for (int i=0;i<m.rows;++i){
        if (m.at<float>(i,0)==0){
            mask.at<uchar>(i,0) = 0;
        }else{
            min =std::min(min, m.at<float>(i,0));
        }
    }
    float min;*/
    double  minV,maxV;
    cv::minMaxLoc(m, &minV,&maxV);
    float range[] = {(float)minV,(float)maxV}; //min
    const float* ranges[] = { range};
    cv::Mat hist;

    int channels[] = {0};
    int  histSize[] = {nBins};
    //    cv::calcHist(&m,1,channels,mask,hist,1 ,histSize, ranges, true,false);
    cv::calcHist(&m,1,channels,cv::Mat(),hist,1 ,histSize, ranges, true,false);
    // cv::norm(proj_k-shape[sk])std::ostringstream oss;
    //oss<<"percentile ";
    //oss <<hist ;
    //qDebug() << oss.str().c_str();

    float valLoc  = cv::sum(hist).val[0]*p;
    float ct = 0;
    for (int i=0; i < nBins;++i){
        ct+= hist.at<float>(i);
        if (valLoc < ct){
            return ((maxV-minV)/float(nBins))*i;//min
        }
    }
    return 1;
}

void randomK(std::vector<int>& random_k, int k, int N) {
  std::vector<int> random_N(N);
  for (int i = 0; i < N; ++ i)
    random_N[i] = i;

  for (int i = 0; i < k; ++ i) {
    int idx = i + rand()%(N-i);
    std::swap(random_N[i], random_N[idx]);
  }

  random_k.assign(k, 0);
  for (int i = 0; i < k; ++ i)
    random_k[i] = random_N[i];

  return;
}
cv::Point2f unitDirection(float ang){
	
  cv::Point2f dir(0.f, 0.f);

  if (ang >= 0) {
    float rotCCW = M_PI * (360- ang)/180.f;
    dir.x=-sin(rotCCW);
    dir.y = -cos(rotCCW);
  }else{
	  qDebug() <<"error: negative angle " << ang << endl;
  }
  return dir;
}
float cosAngleBetween(const cv::Point2f &v1, const cv::Point2f &v2){
	float len1 = sqrt(v1.x * v1.x + v1.y * v1.y);
	float len2 = sqrt(v2.x * v2.x + v2.y * v2.y);

	float dot = v1.x * v2.x + v1.y * v2.y;
	if (len1*len2 ==0){
		return 0;
	}else{
		return dot / (len1 * len2);
	}
}
	
double polylineLength(const std::vector<cv::Point2f> &poly,
                                       int pos1, int pos2){
    double d = 0.f;
    for (int c =std::max(0,pos1); c < std::min((int)poly.size()-1,pos2);++c){
        d+= cv::norm(poly[c] -poly[c+1]);
    }
    return d;
}

double polylineLength(const std::vector<cv::Point2f> &poly){
    double d= polylineLength(poly,0,poly.size()-1);
	/*if (d<1e-6){
		qDebug() <<"Error: invalid polyline length! polyline : ";
		for (int c= 0;c< poly.size();c++){
			qDebug() << "("<<poly[c].x<<"," <<poly[c].y<<") ";
		}
		qDebug()<< endl;
		}*/
	return d;
}

std::string printIntArray(std::vector<int> array){
	std::ostringstream oss;
	for (size_t j=0;j<array.size();++j){
		oss <<" " << array[j];
	}
	return oss.str();
}
std::string printIntSet(std::set<int> array){
    std::ostringstream oss;
    for (auto it=array.begin(); it != array.end(); ++it){
        oss <<" " << *it;
    }
    return oss.str();


}

bool normalize(const cv::Point2f &vin, cv::Point2f &vout)
{
		float norm = cv::norm(vin);
		if (norm == 0){
			return false;
		}
		vout = vin *float (1.f/norm);
			
		return true;
		
	}

}

float variance1D(std::vector<float> data,float mean )
{
    float sum = 0.f;
    for (size_t i=0;i<data.size();++i){
        sum+= pow(data[i]-mean,2);
    }
    return sum/float(data.size());
}

namespace boost {
void throw_exception( std::exception const & e ) {
  throw e;
}
}
