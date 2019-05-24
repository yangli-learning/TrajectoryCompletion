#pragma once
#ifndef SCENE_IMAGE_H
#define SCENE_IMAGE_H

#include <memory>

#include "renderable.h"

class BBox2D;
class GpsTraj;
class TrajPoint;

class GPSTrajectories: public Renderable {
 public:
  const	uint32_t MAX_TRAJ = 5000;
  GPSTrajectories(void);
  virtual ~GPSTrajectories(void);

  bool load(const std::string& filename);

  typedef std::vector<std::shared_ptr<GpsTraj> > GpsTrajs;
  const GpsTrajs& getGpsTrajs(void) const {
    return gps_trajs_;
  }
  int size(){return gps_trajs_.size();}
  
  osg::Vec2 getCenter(void) const {
    return (min_corner_ + max_corner_) / 2;
  }
  const osg::Vec2& getMinCorner(void) const {
    return min_corner_;
  }

  const osg::Vec2& getMaxCorner(void) const {
    return max_corner_;
  }

  float getWidth(void) const {
    return max_corner_.x() - min_corner_.x();
  }

  float getHeight(void) const {
    return max_corner_.y() - min_corner_.y();
  }

  void toggleRenderTrajPoints(bool show_traj_points);

protected:
  virtual void updateImpl(void);

private:
  GpsTrajs gps_trajs_;
  osg::Vec2 min_corner_, max_corner_;

  bool show_traj_points_;
  bool loadPBF(const std::string& filename);
  bool loadTXT(const std::string& filename);

 };

#endif // SCENE_IMAGE_H
