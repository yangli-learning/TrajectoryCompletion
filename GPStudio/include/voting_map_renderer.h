#pragma once
#ifndef VOTING_RENDERER_H
#define VOTING_RENDERER_H

#include "common.h"
#include "l_shape_iterator.h"

#include <mutex>
#include <QReadWriteLock>
class BBox2D;
class VotingMapRenderer{

    friend class BBox2D;

public:
    /**
     * @brief The ToneMapMethod enum
     */
    enum ToneMapMethod {NORMALIZE, // equaivalent to gamma =1
                        GAMMA_COMPRESSION //default gamma = 1.5
    };

    /**
     * @brief VotingMapRender constructor
     */
    VotingMapRenderer(BBox2D* bbox):bbox_(bbox),valid_(false){}

    void init();
    bool isValid(){return valid_;}
    void setValid(bool v){valid_=v;}
    //===========compute voting map ======================

    /**
     * @brief voteMapImageLShape vote map by straight or
     * L-shape lines
     */
    void voteMapImage(void);

    //===========update voting map image================

    void updateMapColorImage(ToneMapMethod method);

    void updateMapColorImageNormalize();

    void updateMapColorImageGammaCompression(float gamma=1.5);

    //============import and export======================
    /**
     * @brief export the heading direction of trajectory point as
     *    images.
     */
     void exportHeadingField(const std::string& basename);

     /**
      * @brief exportVotingImage export the voting result as image
      */
     void exportVotingImage(const std::string& fname);

     //============ utility functions====================
     void normalizeHeadings();
     void blurHeadings();
     std::vector<float> getHeadingVectorAtPos(int x, int y);

     /**
      * @brief interpolateAnglesLShape interpolate between angle a1
      * and a2 for n steps. The 1st,2nd,..,(mid-1)th angles equals a1,
      * the rest equals a2
      * @param a1 starting angle
      * @param a2 ending angle
      * @param n number of steps
      * @param mid index where the L shape bents
      * @return list of interpolated angles
      */
     std::vector<float> interpolateAnglesLShape(float a1, float a2,
                                              int n, int mid);
     /**
      * @brief interpolateAngles returns n linearly interpolated angles
      * between a1 and a2
      */
     std::vector<float> interpolateAngles(float a1, float a2, int n);

     bool insideMap(cv::Point pos) {
       return pos.x >=0 && pos.x < map_image_.cols
              && pos.y >= 0 && pos.y < map_image_.rows;
     }


protected:
    /**
     * @brief map_image_ voting map
     */
    cv::Mat map_image_;

    /**
     * @brief map_color_image_ transformed voting map for display
     * and image processing
     */
    cv::Mat map_color_image_;

    /**
     * @brief directions_ unit vectors for a discrete set
     * of angles that partitions [0,2Pi).
     * directions_[0]=(1,0)
     */
    std::vector<cv::Point2f> directions_;

    /**
     * @brief heading_images_
     */
    std::vector<cv::Mat> heading_images_;

    /**
     * @brief bbox_ BBox2D object that contains the voting_map
     */
    BBox2D *bbox_;
    bool valid_;
};

#endif
