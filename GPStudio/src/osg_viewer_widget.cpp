#define _USE_MATH_DEFINES
#include <cmath>

#include <QFileDialog>
#include <QTextStream>
#include <QResizeEvent>

#include <osg/Depth>
#include <osg/Point>
#include <osg/LineWidth>
#include <osg/ComputeBoundsVisitor>
#include <osgGA/TrackballManipulator>

#include "information.h"
#include "main_window.h"
#include "osg_utility.h"
#include "axis_indicator.h"
#include "toggle_handler.h"
#include "update_callback.h"
#include "stateset_manipulator.h"
#include "force_update_visitor.h"

#include "osg_viewer_widget.h"

std::mutex OSGViewerWidget::mutex_first_frame_;

OSGViewerWidget::OSGViewerWidget(QWidget * parent, const QGLWidget * shareWidget, Qt::WindowFlags flags, float fovy)
  :AdapterWidget(parent, shareWidget, flags, true),
   scene_root_(new osg::Group),
   other_root_(new osg::Group),
   gl_thread_(this),
   threaded_painter_(this),
   first_frame_(true),
   information_(new Information),
   axis_indicator_(new AxisIndicator(this)),
   up_vector_(0.0f, 1.0f, 0.0f) {
  osg::ref_ptr<osg::Group> root = new osg::Group;
  setSceneData(root);

  root->addChild(other_root_);
  root->addChild(scene_root_);
  scene_root_->getOrCreateStateSet()->setAttribute(new osg::Point(2.0f), osg::StateAttribute::OFF);
  scene_root_->getOrCreateStateSet()->setAttribute(new osg::LineWidth(1.0f), osg::StateAttribute::OFF);
  osg::ref_ptr<osgViewer::ScreenCaptureHandler> sch =
	  new osgViewer::ScreenCaptureHandler(new osgViewer::ScreenCaptureHandler::WriteToFile("capture.png", "png"));
  sch->setKeyEventTakeScreenShot('s');
  addEventHandler(sch);
  addEventHandler(new ToggleHandler(information_, 'i', "Toggle Information."));
  addOtherChild(information_);

  addEventHandler(new ToggleHandler(axis_indicator_, 'x', "Toggle Axis Indicator."));
  addOtherChild(axis_indicator_);

  setCameraManipulator(new osgGA::TrackballManipulator);

  addEventHandler(new osgViewer::HelpHandler);
  //addEventHandler(new osgViewer::StatsHandler);

  double w = width();
  double h = height();
  getCamera()->setViewport(new osg::Viewport(0, 0, w, h));
  getCamera()->setProjectionMatrixAsPerspective(fovy, w/h, 1.0f, 10000.0f);
  getCamera()->setGraphicsContext(getGraphicsWindow());
  getCamera()->setClearColor(osg::Vec4(1, 1, 1, 1.0));

  setThreadingModel(osgViewer::Viewer::SingleThreaded);

  this->doneCurrent();
}

OSGViewerWidget::~OSGViewerWidget() {
  osgViewer::View::EventHandlers handlers = getEventHandlers();
  for (osgViewer::View::EventHandlers::iterator it = handlers.begin(); it != handlers.end(); ++ it)
    this->removeEventHandler(*it);

  stopRendering();
}

void OSGViewerWidget::eventTraversal() {
  double reference_time = _frameStamp->getReferenceTime();
  _frameStamp->setReferenceTime(reference_time*100);

  osgViewer::Viewer::eventTraversal();

  _frameStamp->setReferenceTime(reference_time);

  return;
}

void OSGViewerWidget::increaseLineWidth(void) {
  osg::StateSet* state_Set = scene_root_->getOrCreateStateSet();
  osg::LineWidth* line_width = dynamic_cast<osg::LineWidth*>(state_Set->getAttribute(osg::StateAttribute::LINEWIDTH));
  if (line_width == nullptr)
    return;

  if (line_width->getWidth() >= 16.0) {
    line_width->setWidth(16.0);
    return;
  }

  line_width->setWidth(line_width->getWidth()+1.0);

  return;
}

void OSGViewerWidget::decreaseLineWidth(void) {
  osg::StateSet* state_Set = scene_root_->getOrCreateStateSet();
  osg::LineWidth* line_width = dynamic_cast<osg::LineWidth*>(state_Set->getAttribute(osg::StateAttribute::LINEWIDTH));
  if (line_width == nullptr)
    return;

  if (line_width->getWidth() <= 1.0) {
    line_width->setWidth(1.0);
    return;
  }

  line_width->setWidth(line_width->getWidth()-1.0);

  return;
}

void OSGViewerWidget::increasePointSize(void) {
  osg::StateSet* state_set = scene_root_->getOrCreateStateSet();
  osg::Point* point = dynamic_cast<osg::Point*>(state_set->getAttribute(osg::StateAttribute::POINT));
  if (point == nullptr)
    return;

  if (point->getSize() >= 16.0) {
    point->setSize(16.0);
    return;
  }

  point->setSize(point->getSize()+1.0);

  return;
}

void OSGViewerWidget::decreasePointSize(void) {
  osg::StateSet* state_Set = scene_root_->getOrCreateStateSet();
  osg::Point* point = dynamic_cast<osg::Point*>(state_Set->getAttribute(osg::StateAttribute::POINT));
  if (point == nullptr)
    return;

  if (point->getSize() <= 1.0) {
    point->setSize(1.0);
    return;
  }

  point->setSize(point->getSize()-1.0);

  return;
}

void OSGViewerWidget::startRendering() {
  addEventHandler(new StateSetManipulator(getSceneData()->getOrCreateStateSet()));

  context()->moveToThread(&gl_thread_);
  threaded_painter_.moveToThread(&gl_thread_);
  connect(&gl_thread_, SIGNAL(started()), &threaded_painter_, SLOT(start()));
  gl_thread_.start();
}

void OSGViewerWidget::stopRendering() {
  threaded_painter_.stop();
  gl_thread_.wait();
}

void OSGViewerWidget::paintGL(void) {
  std::lock_guard<std::mutex> lock(mutex_);

  if (first_frame_) {
    first_frame_ = false;
    std::lock_guard<std::mutex> lock(mutex_first_frame_);
    frame();
  } else {
    frame();
  }

  return;
}

void OSGViewerWidget::paintEvent(QPaintEvent * /*event*/) {
  // Handled by the GLThread.
}

void OSGViewerWidget::closeEvent(QCloseEvent *event) {
  stopRendering();
  QGLWidget::closeEvent(event);
}

osg::BoundingSphere OSGViewerWidget::getBoundingSphere(void) const {
  return scene_root_->getBound();
}

osg::BoundingBox OSGViewerWidget::getBoundingBox(void) const {
  ForceUpdateVisitor force_update_visitor;
  scene_root_->accept(force_update_visitor);

  osg::ComputeBoundsVisitor visitor;
  scene_root_->accept(visitor);

  return visitor.getBoundingBox();
}

void OSGViewerWidget::addSceneChild(osg::Node *child) {
  std::lock_guard<std::mutex> lock(mutex_);

  scene_root_->addChild(child);

  return;
}

void OSGViewerWidget::addOtherChild(osg::Node *child) {
  std::lock_guard<std::mutex> lock(mutex_);

  other_root_->addChild(child);

  return;
}

void OSGViewerWidget::removeSceneChild(osg::Node *child) {
  std::lock_guard<std::mutex> lock(mutex_);

  scene_root_->removeChild(child);

  return;
}

void OSGViewerWidget::removeOtherChild(osg::Node *child) {
  std::lock_guard<std::mutex> lock(mutex_);

  other_root_->removeChild(child);

  return;
}

void OSGViewerWidget::removeSceneChildren(void) {
  std::lock_guard<std::mutex> lock(mutex_);

  if (scene_root_->getNumChildren() != 0)
    scene_root_->removeChildren(0, scene_root_->getNumChildren());

  return;
}

void OSGViewerWidget::replaceSceneChild(osg::Node* old_child, osg::Node* new_child) {
  std::lock_guard<std::mutex> lock(mutex_);

  scene_root_->removeChild(old_child);
  scene_root_->addChild(new_child);

  return;
}

bool OSGViewerWidget::hasSceneChild(osg::Node *child) {
  std::lock_guard<std::mutex> lock(mutex_);

  return (scene_root_->getChildIndex(child) != -1);
}

bool OSGViewerWidget::hasOtherChild(osg::Node *child) {
  std::lock_guard<std::mutex> lock(mutex_);

  return (other_root_->getChildIndex(child) != -1);
}

void OSGViewerWidget::slotSaveCamera(void) {
  MainWindow* main_window = MainWindow::getInstance();
  QString filename = QFileDialog::getSaveFileName(main_window, "Save Camera Parameters",
                     main_window->getWorkspace().c_str(), "Camera Paramters(*.txt)");
  if (filename.isEmpty())
    return;

  filename = QFileInfo(filename).path() + "/" + QFileInfo(filename).baseName() +".txt";

  QFile txt_file(filename);
  txt_file.open(QIODevice::WriteOnly | QIODevice::Text);
  QTextStream txt_file_stream(&txt_file);

  const osg::Matrixd& manipulator_matrix = getCameraManipulator()->getMatrix();
  for (size_t i = 0; i < 4; ++ i)
    for (size_t j = 0; j < 4; ++ j)
      txt_file_stream << manipulator_matrix(i, j) << ((j == 3)?("\n"):(" "));

  txt_file_stream << "\n";

  const osg::Matrixd& camera_matrix = getCamera()->getViewMatrix();
  for (size_t i = 0; i < 4; ++ i)
    for (size_t j = 0; j < 4; ++ j)
      txt_file_stream << camera_matrix(i, j) << ((j == 3)?("\n"):(" "));

  return;
}

void OSGViewerWidget::savePOVRayCamera(const std::string& filename) {
  std::ofstream fout(filename);

  osg::Vec3 eye, center, up;
  osg::Matrix matrix = getCamera()->getViewMatrix();
  matrix.getLookAt(eye, center, up);

  double fovy, aspect_ratio, z_near, z_far;
  getCamera()->getProjectionMatrixAsPerspective(fovy, aspect_ratio, z_near, z_far);
  double angle = 2*std::atan(aspect_ratio*std::tan((fovy/2)*M_PI/180))*180/M_PI;

  fout << "\n\n//camera\n"
       << "camera\n{\n"
       << "  perspective\n"
       << "  right -x*image_width/image_height\n"
       << "  location <" << eye[0] << "," << eye[1] << "," << eye[2] << ">\n"
       << "  sky <" << up[0] << "," << up[1] << "," << up[2] << ">\n"
       << "  look_at <" << center[0] << "," << center[1] << "," << center[2] << ">\n"
       << "  angle " << angle << "\n"
       << "}\n";

  double radius = scene_root_->getBound().radius();
  double fade_distance = 4*radius;
  double area_light_size = radius/2;
  double eye_offset = radius/2;
  fout << "\n\n//light\n"
       << "#declare light_value = 1.0;\n"
       << "light_source\n{\n"
       << "  <" << eye[0] << "," << eye[1]+eye_offset << "," << eye[2] << ">\n"
       << "  color rgb <light_value, light_value, light_value>\n"
       << "  fade_distance "<< fade_distance << " fade_power 2\n"
       << "  area_light x*" << area_light_size << ", z*" << area_light_size << ", 8, 8 jitter orient circular adaptive 0\n"
       << ""
       << "}\n";

  return;
}

void OSGViewerWidget::slotLoadCamera(void) {
  MainWindow* main_window = MainWindow::getInstance();
  QString filename = QFileDialog::getOpenFileName(main_window, "Load Camera Parameters",
                     main_window->getWorkspace().c_str(), "Camera Parameters (*.txt)");
  if (filename.isEmpty())
    return;

  QFile txt_file(filename);
  txt_file.open(QIODevice::ReadOnly | QIODevice::Text);
  QTextStream txt_file_stream(&txt_file);

  osg::Matrixd manipulator_matrix;
  for (size_t i = 0; i < 4; ++ i)
    for (size_t j = 0; j < 4; ++ j)
      txt_file_stream >> manipulator_matrix(i, j);
  getCameraManipulator()->setByMatrix(manipulator_matrix);

  osg::Matrixd camera_matrix;
  for (size_t i = 0; i < 4; ++ i)
    for (size_t j = 0; j < 4; ++ j)
      txt_file_stream >> camera_matrix(i, j);
  getCamera()->setViewMatrix(camera_matrix);

  return;
}

void OSGViewerWidget::keyPressEvent(QKeyEvent* event) {
  if (event->key() == Qt::Key_Equal) {
    if (event->modifiers() == Qt::CTRL) {
      increasePointSize();
    } else if (event->modifiers() == Qt::SHIFT) {
      increaseLineWidth();
    }
  } else if (event->key() == Qt::Key_Minus) {
    if (event->modifiers() == Qt::CTRL) {
      decreasePointSize();
    } else if (event->modifiers() == Qt::SHIFT) {
      decreaseLineWidth();
    }
  } else if (event->key() == Qt::Key_8) {
    centerScene();
  }

  AdapterWidget::keyPressEvent(event);

  return;
}

void OSGViewerWidget::prepareContextMenu(QMenu* menu) {
  QMenu* osg_viewer_menu = menu->addMenu("OSG Viewer");
  osg_viewer_menu->addAction("Load Camera", this, SLOT(slotLoadCamera()));
  osg_viewer_menu->addAction("Save Camera", this, SLOT(slotSaveCamera()));
  menu->addMenu(osg_viewer_menu);

  return;
}

void OSGViewerWidget::mouseDoubleClickEvent(QMouseEvent* event) {
  QMenu menu(this);

  prepareContextMenu(&menu);

  menu.exec(event->globalPos());

  AdapterWidget::mousePressEvent(event);

  return;
}

void OSGViewerWidget::updateCamera(const osg::Matrix& matrix) {
  osg::Vec3f eye, center, up;
  matrix.getLookAt(eye, center, up);
  osgGA::CameraManipulator* camera_manipulator = getCameraManipulator();
  camera_manipulator->setHomePosition(eye, center, up);
  camera_manipulator->home(0);

  return;
}
