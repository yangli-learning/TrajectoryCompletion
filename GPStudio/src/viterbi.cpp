#include "viterbi.h"
#include <QDebug>
#include <sstream>
#include "parameters.h"
Viterbi::Viterbi(cv::Mat  emission, std::vector<cv::Mat> transition){
    //em￼ithv =  D[i][k] + mu*Q[i-1][j][k] + D[i-1][j];ttps://slack.com/_p_(emission),trans_p_(transition){
	mu_ = Parameters::getInstance()->skeleton_junction_smoothness;
    init(emission,transition);
	//forward();


}
void Viterbi::init(cv::Mat  emission, std::vector<cv::Mat> transition){
    // normalization emission and transition matrix
	num_states_ = emission.rows;
	max_state_size_=emission.cols;
	  

    emit_p_=cv::Mat(emission.rows,emission.cols,CV_32F ,1e6);
	/* skip normalization for emission cost*/

   for (size_t i=0; i< num_states_;++i){
	   	   size_t c;// = (i==0)? transition[i].rows : transition[i-1].cols;
		   float sum=0.f;
	   for (c=0; c <emission.cols;++c){
		   if (emission.at<float>(i,c)>=1e5){
			   break;
		   }
		   sum += emission.at<float>(i,c);
	   }
	   if (sum==0) continue;
	   //cv::Mat temp,rp;	   
	   //cv::normalize(emission.colRange(0,c).row(i),temp);
	   for (int k=0; k<c;++k){
		   emit_p_.at<float>(i,k) = emission.at<float>(i,k)/fabs(sum);
	   }
   }
   
   std::ostringstream oss;
   oss << emit_p_;
   qDebug() << "====================\nemit_p:";
   qDebug() <<  oss.str().c_str() ;
   cv::Mat temp,rp;	   
   trans_p_.resize(transition.size());
   for (size_t i=0; i < transition.size();++i){
	   
       trans_p_[i] = cv::Mat(transition[i].rows,transition[i].cols,CV_32F );
       for (size_t j=0; j < trans_p_[i].rows;++j){
           cv::normalize( transition[i].row(j),temp);
           rp =  trans_p_[i].rowRange(j,j+1);
           temp.copyTo(rp);
       }
	    oss.str("");
	    oss << trans_p_[i] ;
	    qDebug() << "trans_p[" <<i <<"]:";
        qDebug() <<  oss.str().c_str() ;
   }
}

bool Viterbi::forward(){
	/* Note: need to consider invalid values (making them infinite?)
	 * Need to make sure that optimal values are SMALLER since
	 * we are minimizing sums
	 */
	//    qDebug() <<"compute forward viterbi process.";
    if (num_states_ <=2){
        return false;
    }
    double opt_energy = 0.0;
    energy_.resize(num_states_);
	energy_[0].resize(trans_p_[0].rows);
	// initialize energyValues
	for (int j=0; j< trans_p_[0].rows; j++){
		energy_[0][j].index = -1;
		energy_[0][j].val  =emit_p_.at<float>(0,j);
	}
    // ----------forward step--------------------------

    for (int i=1;i<num_states_;i++){
		
        for (int k=0; k<trans_p_[i-1].cols;k++){
			energy_[i].resize( trans_p_[i-1].cols);//equivalent to tr[i-1].cols
            energy_[i][k].index = -1;
            double minVal = Common::inf;
            for (int j=0; j< trans_p_[i-1].rows;j++){

				// compute v =  D[i][k] + mu*Q[i-1][j][k] + D[i-1][j];
                double val =  energy_[i-1][j].val + mu_* trans_p_[i-1].at<float>(j,k) ;
                
                if (val < minVal){
                    minVal = val;
                    energy_[i][k].index =j;
                }
            }
			//            if ( (*candidate_ids)[i-1].size() >0){
			energy_[i][k].val= emit_p_.at<float>(i,k) + (float)minVal;

        }

    }
	printEnergy();
	
	//------backward step --------------
	opt_selection_= std::vector<int>(num_states_,-1);
	int endpoint = num_states_-1;
	while(energy_[endpoint].size()==0 && endpoint >=0){
	    endpoint--;
	}
	if (endpoint<0)
		return false;
	// find minimizing candidate id and energy value at the last step
	// what is the row dimension of energy?? is it n or n-1??!!
	// what is the comparsion function of indexVal? (needs to be conventional)
	//!
	//! DEBUG by printing stats of energy matrix?
	//!
	auto iter = std::min_element(energy_[endpoint].begin(),
								 energy_[endpoint].end());
	opt_selection_[endpoint] =  iter-energy_[endpoint].begin();//minId;
	opt_energy = iter->val;
	int minId = iter->index;//distance(energy_[num_states_-1].begin(),iter);// iter-energy_[num_states_-1].begin();
	
	/*for (size_t i=0;i< energy_[num_states_-1].size();i++){
		if (  energy_[num_states_-1][i].second < minVal ){
				minId = energy_[endpoint][i].first;
		}
		}*/
	//	opt_energy = minVal; //global variable that may be used later


	for (int i=endpoint-1; i>=0;i--){
		if (minId >=0 && minId < energy_[i].size()){

			opt_selection_[i]= minId;//(minId>=0)?minId:-1;
			minId = energy_[i][minId].index; //update minId by backpointer			
		}else{
			//opt_selection_[i] = -1;
			qDebug()<<"Error: invalid path assignment "
					<< minId << " found from row**" << i;
			while( energy_[i].size()==0 && i >0){
				i--;
			}
			if (energy_[i].size() >0){
				iter = std::min_element( energy_[i].begin(), energy_[i].end());
				opt_selection_[i]= iter-energy_[i].begin();//(minId>=0)?minId:-1;
				minId = iter->index;
			}

		}
	}
	//output the optimal path
	
	std::ostringstream oss;
	oss <<"[Optimal skeleton branch projection]: ";
	for (size_t s=0;s<opt_selection_.size();++s){
		oss << opt_selection_[s]<<" ";
	}
	qDebug()<<oss.str().c_str();
	
	return true;
}

void Viterbi::printEnergy(){
	std::ostringstream oss;
	qDebug()<<"-------------------------------";
	qDebug() << "[Min sum energy table]: ";
	for (size_t i=0;i<energy_.size();++i){
		oss.str("");
		for (size_t j=0;j<energy_[i].size();++j){
			oss << energy_[i][j].val <<"|("
				<< energy_[i][j].index <<") ";
		}
		qDebug()<<i<<": "<< oss.str().c_str();
	}
	qDebug()<<"-------------------------------";
	
}


