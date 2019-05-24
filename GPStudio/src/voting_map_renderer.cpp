#include "voting_map_renderer.h"
#include "parameters.h"
#include "l_shape_iterator.h"
#include <QDebug>
#include "bbox_2d.h"
#include "gps_trajectories.h"
#include "gps_trajectory.pb.h"


#if (CV_VERSION_EPOCH  > 2)
#include <opencv2/imgcodecs.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#else
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
//#include <opencv2/highgui/highgui_c.h>
#endif

void VotingMapRenderer::init(){//int width,int height){
    bbox_->expire();
    //bbox_->forceUpdate();
    int width = ceil(bbox_->getWidth()/bbox_->getResolution());
    int height =ceil(bbox_->getHeight()/bbox_->getResolution());

    int nDirections = Parameters::getInstance()->n_directions;
    for (int i=0; i<nDirections;++i){
      cv::Point2f p = LShapeIterator::UnitDirection(i*float(360/nDirections));
      directions_.push_back(p);
      qDebug() << "angle="  << i*2*CV_PI/nDirections <<" vector=(" << p.x <<","<<p.y << ")";
    }
    heading_images_.resize(directions_.size());
    if (width > 3000 || height > 3000) {
        fprintf(stderr, "ERROR! Cannot create map image with %d x %d dimensions!\n", width, height);
        return;
    }

    map_image_ = cv::Mat(height, width, CV_32FC1, cv::Scalar(0));
    map_color_image_ = cv::Mat(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

    for (size_t i=0; i< heading_images_.size();++i){
        heading_images_[i] = cv::Mat(height, width, CV_32FC1, cv::Scalar(0));
    }
    qDebug()<<"++ finished initializing bbox..";
}


void VotingMapRenderer::updateMapColorImage(ToneMapMethod m) {
  switch(m) {
  case NORMALIZE:
    updateMapColorImageNormalize();
    break;
  case GAMMA_COMPRESSION:
    updateMapColorImageGammaCompression();
    break;
  default:
    break;
  }
  qDebug()<< "++ update map color image...";
  return;
}
void VotingMapRenderer::updateMapColorImageGammaCompression(float gamma) {
  double min, max;
  cv::minMaxLoc(map_image_, &min, &max);
  if (min == max) {
    map_color_image_.setTo(cv::Scalar(255, 255, 255));
    return;
  }

  double A = pow(max-min,-gamma);

  for (int y = 0; y < map_image_.rows; ++y) {
    const float* p_map_image = map_image_.ptr<float>(y);
    uchar* p_map_color_image = map_color_image_.ptr<uchar>(y);
    for (int x = 0; x < map_image_.cols; ++x) {
      float value = p_map_image[x];
      uchar gray =255-(A * pow(value-min,gamma))*255;
      p_map_color_image[3 * x + 0] = gray;
      p_map_color_image[3 * x + 1] = gray;
      p_map_color_image[3 * x + 2] = gray;
    }
  }
}

void VotingMapRenderer::updateMapColorImageNormalize(void) {
    updateMapColorImageGammaCompression(1.0f);
}

void VotingMapRenderer::voteMapImage(void){
  bbox_->expire();
    map_image_ = cv::Scalar(0.f);
    for (int i=0; i<heading_images_.size();++i){
        heading_images_[i] =cv::Scalar(0.f);
    }
    //heading_images_.clear();
  qDebug()<<" Vote map image with l-shape.. " << endl;

  GPSTrajectories::GpsTrajs gps_trajs = bbox_->getTrajectories();
  //gps_trajectories_->getGpsTrajs();

  for (size_t i = 0, i_end = gps_trajs.size(); i < i_end; ++i) {
    std::shared_ptr<GpsTraj>& gps_traj = gps_trajs[i];
    if (gps_traj->point_size() < 2)
      continue;
    int ang_prev, ang_curr;
    try {
      ang_prev = (int32_t) (gps_traj->point(0).head());
    } catch(...) {
      continue; //skip unreadable point
    }
    cv::Point point_prev = bbox_->toCvPoint(&(gps_traj->point(0)));
    for (size_t j = 1, j_end = gps_traj->point_size(); j < j_end; ++j) {

      try {
        ang_curr= (int32_t) (gps_traj->point(j).head());
      } catch(...) {
        continue; //skip unreadable points
      }
      cv::Point point_curr = bbox_->toCvPoint(&(gps_traj->point(j)));
      if (point_prev == point_curr)
        continue; //! @todo for simplicity, skip stopping points


      LShapeIterator l_iterator(map_image_, point_prev, point_curr,
                                ang_prev, ang_curr);
      // store the indices the iterator passes through for further
      // processing
      std::vector<cv::Point > iterIndex;

      for (int k = 0, k_end = (j == j_end - 1) ? (l_iterator.count) : (l_iterator.count - 1); k < k_end; ++k, ++l_iterator) {
          if ( insideMap(l_iterator.pos())){
              float* p_value = (float *) (*l_iterator);
              if (p_value){
                  *p_value += 1.0;
                  iterIndex.push_back(l_iterator.pos());
              }else{
                  iterIndex.push_back(cv::Point(-1,-1));
              }
          }else{
              iterIndex.push_back(cv::Point(-1,-1));
          }
      }
      // set heading angles of interpolated samples
      std::vector<float> headings;
      if (Parameters::getInstance()->heading_linear_interpolation){
          headings = interpolateAngles(ang_prev, ang_curr,
                                                        l_iterator.count);
      }else{
          headings = interpolateAnglesLShape(ang_prev, ang_curr,
                                                  l_iterator.count,
                                                  l_iterator.getPivot());
      }
      // compute votes for each direction
      for (size_t id =0;id<iterIndex.size();++id){
          int x = iterIndex[id].x, y = iterIndex[id].y;
          if (x==-1 || y==-1) continue;
           cv::Point2f h = LShapeIterator::UnitDirection(headings[id]);
           if (directions_.size()==2){
               heading_images_[0].at<float>(y,x) += h.x;
               heading_images_[1].at<float>(y,x) += h.y;

          }else{
              for (size_t d = 0; d < directions_.size(); ++d){
                  float dotprod = h.dot(directions_[d]);
                  if (dotprod >0)
                      heading_images_[d].at<float>(y,x) += h.dot(directions_[d]);

              }
          }
      }


      point_prev = point_curr;
      ang_prev = ang_curr;
    }
  }
  // need to normalize heading?

  if (Parameters::getInstance()->heading_blur){
    blurHeadings();
  }
  if (Parameters::getInstance()->skeleton_normalize_heading ){
    normalizeHeadings();
  }
  valid_ = true;
  return;
}

void VotingMapRenderer::exportHeadingField(const std::string& basename){
    qDebug() << "export heading field to files: " << basename.c_str() << endl;
    double min, max;

    for (size_t i=0; i<heading_images_.size();++i){
        int width = heading_images_[i].cols, height = heading_images_[i].rows;
        cv::Mat img(height, width, CV_8UC3, cv::Scalar(255, 255, 255));

        cv::minMaxLoc(heading_images_[i], &min, &max);

        if (min == max) {
            img.setTo(cv::Scalar(255, 255, 255));
            continue;
        }

        for (int y = 0; y <height; ++y) {

            const float* p_row = heading_images_[i].ptr<float>(y);
            uchar* p_color_row = img.ptr<uchar>(y);

            for (int x = 0; x <width; ++x) {

                float value = p_row[x];
                uchar gray = 255-((value - min) / (max-min) * 255);
                p_color_row[3 * x + 0] = gray;
                p_color_row[3 * x + 1] = gray;
                p_color_row[3 * x + 2] = gray;

            }
        }
        std::ostringstream fname;
        fname << basename << "_" <<  floor(i) << ".png";
        qDebug() << "Writing heading field to file " <<  fname.str().c_str()<< endl;

        IplImage* outImg =new IplImage(img);

        cvSaveImage(fname.str().c_str(), outImg);
        delete outImg;

    }
}

void VotingMapRenderer::exportVotingImage(const std::string &fname){
    //updateImpl();
    IplImage* outImg =new IplImage(map_color_image_);
    cvSaveImage(fname.c_str()  ,outImg);
    delete outImg;
    qDebug() << "Writing voting image to " <<  fname.c_str()<< endl;
}

std::vector<float> VotingMapRenderer::getHeadingVectorAtPos(int x, int y){
    std::vector<float> head(directions_.size(),0);
    bool allZero =true;
    for (size_t i=0;i<head.size();++i){
        head[i] =  heading_images_[i].at<float>(y,x) ;
        if (head[i]>1e-4) allZero = false;
    }
    if (allZero){
        for (size_t i=0;i<head.size();++i){
            int ct =0;
            if (y+1 < heading_images_[i].rows){
                head[i]+= heading_images_[i].at<float>(y+1,x);
                ct++;
            }
            if (y-1 >= 0){
                head[i]+=heading_images_[i].at<float>(y-1,x);
                ct++;
            }
            if (x+1 < heading_images_[i].rows){
                head[i]+=heading_images_[i].at<float>(y,x+1);
                ct++;
            }
            if (x-1 >= 0){
                head[i]+=heading_images_[i].at<float>(y,x-1);
                ct++;
            }
            head[i]/=ct;
        }
    }
    return head;
}


std::vector<float>  VotingMapRenderer::interpolateAngles(float a1, float a2, int n){
    std::vector<float> result(n);
    if (n==0) return result;

    float diff = a2-a1+1;
    float step = diff/n;
    for (size_t i=0;i<result.size();++i){
        int mul = floor((a1+i*step) / 360);
        result[i] =  (a1 + i*step) - mul* 360;
    }
}
std::vector<float> VotingMapRenderer::interpolateAnglesLShape(float a1, float a2,
                                                              int n, int mid){
    std::vector<float> result(n);
    if (n==0) return result;

    float diff = a2-a1+1;
    float step = diff/n;

    for (size_t i=0;i<result.size();++i){
        if (i<=mid){
            result[i] = a1;
        }else{
            result[i] = a2;
        }
    }
    return result;
}

void VotingMapRenderer::blurHeadings(void){
    std::vector<cv::Mat>::iterator it = heading_images_.begin();
    for (; it!= heading_images_.end();++it){
        cv::Mat dst;
        //		cv::bilateralFilter(*it,dst,3,6,6);
        cv::GaussianBlur(*it,dst,cv::Size(3,3),0,0);
        *it = dst;
    }
}
void VotingMapRenderer::normalizeHeadings(){
    int height = map_image_.rows, width = map_image_.cols;
    if (heading_images_[0].rows != height || heading_images_[1].cols != width){
        qDebug()<<"Error: problem normalize headings ";
        return;
    }

    cv::Mat zeroMat(map_image_.rows, map_image_.cols,CV_32FC1,cv::Scalar(0));
    for (size_t i=0;i<heading_images_.size();++i){
        heading_images_[i]/= max(map_image_-2,zeroMat);
    }

}

