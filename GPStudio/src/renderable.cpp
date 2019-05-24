#include <random>

#include <QMenu>
#include <QImage>
#include <QColor>
#include <QCursor>

#include <osg/Point>
#include <osg/Geode>
#include <osg/StateSet>
#include <osg/Material>
#include <osgViewer/Viewer>
#include <osgUtil/UpdateVisitor>
#include <osg/ComputeBoundsVisitor>

#include "osg_utility.h"
#include "update_callback.h"
#include "force_update_visitor.h"

#include "renderable.h"

std::mutex Renderable::mutex_graphics_context_;

Renderable::Renderable(void)
  :read_write_lock_(QReadWriteLock::NonRecursive),
   content_root_(new osg::MatrixTransform),
   expired_(true),
   hidden_(false) {
  addChild(content_root_);

  setUpdateCallback(new UpdateCallback);
  setDataVariance(Object::DYNAMIC);

  return;
}

Renderable::~Renderable(void) {
}

void Renderable::expire(void) {
  QWriteLocker locker(&read_write_lock_);
  expired_ = true;
  //forceUpdate();
  return;
}

void Renderable::update(void) {
  if (!read_write_lock_.tryLockForRead())
    return;

  if (!expired_) {
    read_write_lock_.unlock();
    return;
  }

  expired_ = false;
  content_root_->removeChildren(0, content_root_->getNumChildren());

  if (!hidden_)
    updateImpl();

  read_write_lock_.unlock();
  return;
}

void Renderable::forceUpdate(void) {
  read_write_lock_.lockForRead();

  if (!expired_) {
    read_write_lock_.unlock();
    return;
  }

  expired_ = false;
  content_root_->removeChildren(0, content_root_->getNumChildren());

  if (!hidden_)
    updateImpl();

  read_write_lock_.unlock();
  return;
}

void Renderable::toggleHidden(void) {
  QWriteLocker locker(&read_write_lock_);
  expired_ = true;
  hidden_ = !hidden_;

  return;
}

osg::BoundingBox Renderable::getBoundingBox(void) {
  ForceUpdateVisitor force_update_visitor;
  this->accept(force_update_visitor);

  osg::ComputeBoundsVisitor visitor;
  content_root_->accept(visitor);

  return visitor.getBoundingBox();
}

static void saveColorImage(osg::Image* color_image, osg::Image* depth_image, const std::string& filename) {
  int width = color_image->s();
  int height = color_image->t();
  float* z_buffer = (float*)(depth_image->data());

  QImage q_image(width, height, QImage::Format_ARGB32);
  for (int x = 0; x < width; ++ x) {
    for (int y = 0; y < height; ++ y) {
      float z = z_buffer[y*width+x];
      osg::Vec4 color = color_image->getColor(x, y);
      color = color*255;
      QRgb rgba = QColor(color.r(), color.g(), color.b(), (z==1.0)?(0):(255)).rgba();
      q_image.setPixel(x, height-1-y, rgba);
    }
  }

  q_image.save(filename.c_str());

  return;
}

static void saveDepthImage(osg::Image* depth_image, const std::string& filename) {
  int width = depth_image->s();
  int height = depth_image->t();

  float z_min, z_max;
  z_min = std::numeric_limits<float>::max();
  z_max = std::numeric_limits<float>::lowest();
  float* z_buffer = (float*)(depth_image->data());
  for (size_t x = 0; x < width; ++ x) {
    for (size_t y = 0; y < height; ++ y) {
      float z = z_buffer[y*width+x];
      if (z == 1.0)
        continue;

      z_min = std::min(z_min, z);
      z_max = std::max(z_max, z);
    }
  }

  QImage q_image(width, height, QImage::Format_ARGB32);
  for (size_t x = 0; x < width; ++ x) {
    for (size_t y = 0; y < height; ++ y) {
      float z = z_buffer[y*width+x];
      float value = (z==1.0)?(1.0):(z-z_min)*0.8/(z_max-z_min);
      value *= 255;
      QRgb rgba = QColor(value, value, value, (z==1.0)?(0):(255)).rgba();
      q_image.setPixel(x, height-1-y, rgba);
    }
  }

  q_image.save(filename.c_str());

  return;
}

void Renderable::testVisibility(const std::vector<osg::Matrix>& camera_path, int width, int height,
                                const std::vector<osg::Vec3>& positions, std::vector<bool>& flags, float fovy) {
  osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
  traits->x =0;
  traits->y = 0;
  traits->width = width;
  traits->height = height;
  traits->windowDecoration = false;
  traits->doubleBuffer = false;
  traits->sharedContext = 0;
  traits->pbuffer = true;

  mutex_graphics_context_.lock();
  osg::ref_ptr<osg::GraphicsContext> graphics_context= osg::GraphicsContext::createGraphicsContext(traits.get());
  mutex_graphics_context_.unlock();

  osgViewer::Viewer* viewer = new osgViewer::Viewer;
  osg::ref_ptr<osg::Camera> camera = viewer->getCamera();
  camera->setGraphicsContext(graphics_context);
  camera->setViewport(new osg::Viewport(0, 0, width, height));
  camera->setProjectionMatrixAsPerspective(fovy, 1.0f*width/height, 1.0f, 10000.0f);
  camera->setClearColor(osg::Vec4(1, 1, 1, 1.0));

  osg::ref_ptr<osg::Group> node = new osg::Group;
  osg::ref_ptr<osg::StateSet> state_set = node->getOrCreateStateSet();
  state_set->setAttribute(new osg::Point(4.0f), osg::StateAttribute::OFF);
  state_set->setMode(GL_DEPTH_TEST, osg::StateAttribute::ON);
  node->addChild(this);
  viewer->setSceneData(node);
  viewer->setDataVariance(osg::Object::DYNAMIC);
  viewer->setThreadingModel(osgViewer::Viewer::SingleThreaded);
  viewer->realize();

  osg::ref_ptr<osg::Image>  color_image = new osg::Image;
  osg::ref_ptr<osg::Image>  depth_image = new osg::Image;
  color_image->allocateImage(traits->width, traits->height, 1, GL_RGBA, GL_FLOAT);
  depth_image->allocateImage(traits->width, traits->height, 1, GL_DEPTH_COMPONENT, GL_FLOAT);
  camera->attach(osg::Camera::COLOR_BUFFER, color_image.get());
  camera->attach(osg::Camera::DEPTH_BUFFER, depth_image.get());

  double z_buffer_threshold = 0.01;
  flags.assign(positions.size(), false);
  for (size_t i = 0, i_end = camera_path.size(); i < i_end; ++ i) {
    camera->setViewMatrix(camera_path[i]);
    viewer->frame();

    //saveDepthImage(depth_image.get(), "depth.png");
    //saveColorImage(color_image.get(), depth_image.get(), "color.png");

    osg::Matrix matrix_vpw(camera->getViewMatrix() * camera->getProjectionMatrix());
    matrix_vpw.postMult(camera->getViewport()->computeWindowMatrix());
    osg::Matrix matrix_vpw_inverse;
    matrix_vpw_inverse.invert(matrix_vpw);

    float* z_buffer = (float*)(depth_image->data());
    for (size_t j = 0, j_end = positions.size(); j < j_end; ++ j) {
      if (flags[j])
        continue;

      osg::Vec3 window_xy = positions[j]*matrix_vpw;
      int window_x = (int)(window_xy.x());
      int window_y = (int)(window_xy.y());
      // outside of the current view port
      if (window_x < 0 || window_x >= width || window_y < 0 || window_y >= height)
        continue;

      // no valid depth at the position -> not hidden
      double window_z = z_buffer[window_y*width+window_x];
      if (window_z == 1.0) {
        flags[j] = true;
        continue;
      }

      osg::Vec3 eye, center, up;
      camera_path[i].getLookAt(eye, center, up);
      // the position is epsilon close to the visible point -> not hidden
      osg::Vec3 world = osg::Vec3(window_x, window_y, window_z)*matrix_vpw_inverse;
      double depth_difference = (eye-positions[j]).length() - (eye-world).length();
      if (depth_difference < z_buffer_threshold)
        flags[j] = true;
    }
  }

  mutex_graphics_context_.lock();
  delete viewer;
  mutex_graphics_context_.unlock();

  return;
}
