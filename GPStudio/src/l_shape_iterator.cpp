#define _USE_MATH_DEFINES
#include <cmath>

#include "l_shape_iterator.h"
#include <QtDebug>
LShapeIterator::LShapeIterator(const cv::Mat img,
                               cv::Point pt1, cv::Point pt2,
                               int ang1, int ang2) {
  if (HeadingChanged(ang1,ang2)
      && FindIntersection(pt1,pt2,ang1, ang2)
	  && q_.x >=0 && q_.y >=0 && q_.x < img.cols && q_.y < img.rows) {
    // L-shape
    it1_ = new cv::LineIterator(img, pt1, q_);
    it2_ = new cv::LineIterator(img, q_, pt2);

    count1_ = it1_->count;
    count2_ = it2_->count;

  } else {
    // straight line
    it1_ = new cv::LineIterator(img,pt1,pt2);
    it2_ = NULL;
    count1_ = it1_->count;
    count2_ = 0;
  }
  count = count1_ + count2_;
  ptr = it1_->ptr;  // initiate point
}

LShapeIterator::~LShapeIterator() {

  delete it1_;
  delete it2_;
}

cv::Point LShapeIterator::pos() const {

  if (it1_->pos() == q_ && it2_) {
    return it2_->pos();
  } else {
    return it1_->pos();
  }
}

LShapeIterator&  LShapeIterator::operator ++() {
  // the ptr returned is the value in counter after it has been
  // incremented.

  if (it1_->pos() == q_ && it2_) {
    // reach the end of segment 1
    // increment pointer on segment 2
    ++(*it2_);
    ptr = it2_->ptr;
  } else {
    ++(*it1_);
    ptr  =  it1_->ptr;
  }
  return *this;
}


bool  LShapeIterator::FindIntersection(const cv::Point pt1, const cv::Point pt2, 
					int ang1, int ang2){
	// find intersection between rays (pt1,ang1) and (pt2, ang2)
	cv::Point2f o1(pt1.x, pt1.y);
	cv::Point2f o2(pt2.x, pt2.y);
	cv::Point2f x =  o2-o1;
	cv::Point2f d1 = LShapeIterator::UnitDirection(ang1);
	cv::Point2f d2 = LShapeIterator::UnitDirection(ang2);
	float cross = d1.x*d2.y - d1.y*d2.x;
	
	if (abs(cross) < 1e-8)
		return false;
	double t1 = (x.x * d2.y - x.y * d2.x)/cross;
	
	cv::Point2f in = o1 + d1 * t1;
	q_.x = (int)in.x;
	q_.y = (int)in.y;
	
	return true;
}

bool LShapeIterator::HeadingChanged(int ang1,int ang2){
	cv::Point2f d1 = LShapeIterator::UnitDirection(ang1);
	cv::Point2f d2 =  LShapeIterator::UnitDirection(ang2);
	return fabs(d1.x *d2.x + d1.y*d2.y)<0.25;
	//	double angle = acos(d1.x *d2.x + d1.y*d2.y); 


	//return fabs(cos( angleBetween(d1,d2))) < 0.25;
}

/*
bool LShapeIterator::HeadingChanged(int ang1, int ang2){

  int angDiff = std::min( abs(ang1-ang2), 360 - abs(ang1-ang2));

  return angDiff > 30;
}
*/
cv::Point2f  LShapeIterator::UnitDirection(float ang) {

  cv::Point2f dir(0.f, 0.f);

  if (ang >= 0) {
    float rotCCW = M_PI * (360- ang)/180.f;
    dir.x=-sin(rotCCW);
    dir.y = -cos(rotCCW);
  }/*else{
	  qDebug() <<"error: negative angle " << ang << endl;
	  }*/

  return dir;
}
  
