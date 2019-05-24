#pragma once
#ifndef COMMON_H_
#define COMMON_H_

#define _USE_MATH_DEFINES
#include <cmath>
#include <string>
#include <vector>
#include <set>
#include <limits>
#include <sstream>
#include <osg/Vec3>
#if (CV_VERSION_EPOCH  > 2)
#include <opencv2/core.hpp>
#include <opencv2/imgproc.hpp>
#else
#include <opencv2/core/core.hpp>
#include <opencv2/imgproc/imgproc.hpp>
#endif
enum class ColorMode {
	DEFAULT,
		UNIFORM,
		CURVATURE
		};

enum class PickMode {
	CTRLSHIFT,
		CTRL,
		SHIFT,
		ALT,
		ALTSHIFT
		};

/** int value pair ordered by value **/
template <typename T> struct IndexVal{
	int index;
	T val;
	bool operator < (const IndexVal& rhs) const{
		return val < rhs.val;
	}
};
/** unordered int value pair **/
template <typename T> struct IndexValU{
	int index;
	T val;

};
namespace Common {
	const int inf_i = std::numeric_limits<int>::infinity();
	const float inf = std::numeric_limits<float>::infinity();
    const int palette_size = 10;
    const float palette[10][3] = {
	  {166,206,227},
	  {31,120,180},
	  {178,223,138},
	  {51,160,44},
	  {251,154,153},
	  {227,26,28},
	  {253,191,111},
	  {255,127,0},
	  {202,178,214},
	  {106,61,154}
    };

	std::string int2String(int i, int width);
   // cv::Point2f normalize(const cv::Point2f v);
    double polylineLength(const std::vector<cv::Point2f> &poly);
    double polylineLength(const std::vector<cv::Point2f> &poly, int pos1, int pos2);
	void randomK(std::vector<int>& random_k, int k, int N);
	bool normalize(const cv::Point2f &vin, cv::Point2f &vout);
    float variance1D(std::vector<float> data,float mean=0.f);
    float percentile(const cv::Mat &m, float p, int nBins=200);
	float cosAngleBetween(const cv::Point2f &v1, const cv::Point2f &v2);
	cv::Point2f unitDirection(float ang);
	
    template <typename T>
    T rowVariance(const cv::Mat data, int r, T mean){
        T sum = 0.f;

    //    cv::Mat rData = data.row(r);
        int ct=0;
        for (size_t i=0;i< data.cols;++i){
            if (data.at<T>(r,i) >=0){
                sum+= (data.at<T>(r,i) -mean)*(data.at<T>(r,i) -mean);
                ct++;
            }
        }
        return sum/float(data.cols);
    }

	template <typename Output, typename Input>
		Output Vec3Caster(const Input& input) {
		return Output(input.x(), input.y(), input.z());
	}

	template <typename Output, typename Input>
		Output MatrixCaster(const Input& input) {
		Output output;
		for (int i = 0; i < 4; ++ i)
			for (int j = 0; j < 4; ++ j)
				output(i, j) = input(i, j);
		return output;
	}

	template <typename Output, typename Input>
		Output MatrixTransposeCaster(const Input& input) {
		Output output;
		for (int i = 0; i < 4; ++ i)
			for (int j = 0; j < 4; ++ j)
				output(i, j) = input(j, i);
		return output;
	}

	std::string printIntArray(std::vector<int> array); 
    std::string printIntSet(std::set<int> array);
}

#endif // COMMON_H_
