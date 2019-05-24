#pragma once
#ifndef L_ITER_H
#define L_ITER_H

#include <opencv2/core/core.hpp>

#if (CV_VERSION_EPOCH  > 2)
#include <opencv2/imgcodecs.hpp>
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#else
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp> 
#endif
class LShapeIterator{
 public:

	/**
	   Iterator for pixel on the L-shape line segments that connects
	   pt1 to pt2 .
	   
	   img: 2D point map
	   pt1: starting point
	   pt2: final point
	   ang1: heading angle of pt1 with respect to north
	   ang2: heading angle of pt2 with respect to north
	**/
	LShapeIterator(const cv::Mat img, cv::Point pt1, cv::Point pt2,
				   int ang1, int ang2);

	~LShapeIterator();

	/**
	   Returns pointer to the current line pixel
	**/
	inline uchar* operator *(){

		return ptr;
	}

	/** (gps_traj->point(j).head());
		Increment iterator 
	**/
	LShapeIterator& operator ++();


	/**
	   compute a unit vector in the heading direction
	   negative heading angle results in zero vector
	*/
	static cv::Point2f UnitDirection(float ang);


	/**
	   Returns coordinates of the current pixel
	**/
	cv::Point pos() const;

	uchar* ptr;

	/**
	   Total number of pixel steps 
	**/
    int count;

	int getPivot(){return count1_;}
 protected:
	/**
		Finds the intersection of ray pt1, pt2 with
		angles ang1 and ang2 respectively 
	*/
	bool FindIntersection(const cv::Point pt1, 
						  const cv::Point pt2,
						  int ang1, int ang2);
	/**
	   returns true if heading angles differ >45 degrees 
	 */
	bool HeadingChanged(int ang1, int ang2);

	int count1_, count2_; //length of the first and second line
	
	cv::LineIterator *it1_; // first line seg iterator
	cv::LineIterator *it2_; // second line seg iterator
	cv::Point q_; // position of the inflection point in an L shape

};


#endif
