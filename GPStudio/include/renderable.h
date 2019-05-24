#pragma once
#ifndef RENDERABLE_H
#define RENDERABLE_H

#include <mutex>
#include <QReadWriteLock>
#include <osg/BoundingBox>
#include <osg/MatrixTransform>

#include "common.h"

#define META_Renderable(name) \
  virtual bool isSameKindAs(const osg::Object* obj) const { return dynamic_cast<const name *>(obj)!=nullptr; } \
  virtual const char* className() const { return #name; } \
  virtual void accept(osg::NodeVisitor& nv) { if (nv.validNodeMask(*this)) { nv.pushOntoNodePath(this); nv.apply(*this); nv.popFromNodePath(); } } \
 
class UpdateCallback;
class SceneWidget;
class Renderable : public osg::MatrixTransform {
 public:
  Renderable(void);
  virtual ~Renderable(void);

  META_Renderable(Renderable);

  void expire(void);
  inline void setSceneWidget(SceneWidget* sw){
	  scene_widget_ = sw;
  }
  inline bool isHidden(void) const {
    return hidden_;
  }
  virtual void toggleHidden(void);

  QReadWriteLock& getReadWriteLock(void) {
    return read_write_lock_;
  }

  osg::BoundingBox getBoundingBox(void);

  virtual void pickEvent(PickMode pick_mode) {}
  virtual void pickEvent(PickMode pick_mode, float x, float y){}

  void testVisibility(const std::vector<osg::Matrix>& camera_path, int width, int height,
                      const std::vector<osg::Vec3>& positions, std::vector<bool>& flags, float fovy=43.0f);
 protected:
  friend class UpdateCallback;
  void update(void);
  friend class ForceUpdateVisitor;
  void forceUpdate(void);
  virtual void updateImpl(void) = 0;

 protected:
  mutable QReadWriteLock              read_write_lock_;
  osg::ref_ptr<osg::MatrixTransform>  content_root_;
  bool                                expired_;
  bool                                hidden_;
  SceneWidget* scene_widget_; 
  static std::mutex										mutex_graphics_context_;
};

#endif // RENDERABLE_H
