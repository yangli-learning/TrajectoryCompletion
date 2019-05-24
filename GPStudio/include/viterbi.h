#pragma once
#ifndef VITERBI_H
#define VITERBI_H
#include "common.h"
#include <vector>
#include <opencv2/core/core.hpp>

/**
 * class Viterbi
 * implements a viterbi algorithm for solving HMM problem
 * Currently all candidates have fixed size
 * \@todo allow variable candidate sizes in transition_probability
 */
class Viterbi{
public:
    Viterbi( cv::Mat  emission, std::vector<cv::Mat> transition);

    void init(cv::Mat  emission, std::vector<cv::Mat> transition);
    bool forward();
    std::vector<int> getOptimalPath(){return opt_selection_;}

protected:
   // computeD();
    //computeQ();
    size_t num_states_;
    size_t max_state_size_;
    cv::Mat emit_p_;
    std::vector<cv::Mat> trans_p_;

	std::vector< std::vector< IndexVal<double> > > energy_ ;
	float mu_;
	std::vector<int> opt_selection_;


	void printEnergy();
};

/*class JunctionViterbi:Viterbi{

public:
    JunctionViterbi(){

    }
};*/

#endif
