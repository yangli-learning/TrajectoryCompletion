#pragma once
#ifndef SCENE_WIDGET_H_
#define SCENE_WIDGET_H_

#include "osg_viewer_widget.h"
#include <QListWidgetItem>

class QMenu;
class QMouseEvent;
class QWheelEvent;

class BBox2D;
class GPSTrajectories;

class SceneWidget : public OSGViewerWidget {
  Q_OBJECT

 public:
  SceneWidget(QWidget * parent = 0, const QGLWidget * shareWidget = 0, Qt::WindowFlags f = 0);
  virtual ~SceneWidget();

  virtual QSize sizeHint() const {
    return QSize(256, 256);
  }

  virtual void centerScene(void);
  void centerScene2D(void);
  void removeBBox(BBox2D* bbox);
  void processBBox(BBox2D* bbox);
  void saveSelectedClusters(std::vector<bool> selected);
  void toggleClusterDisplay(std::vector<bool> selected,bool noline);
  
  osg::Node * getSceneChild(int i){
      return scene_root_->getChild(i);
  }

  GPSTrajectories* getGPSTrajectories(void) {
    return gps_trajectories_;
  }

 public slots:
  void slotOpenGPSTrajectories(void);

  void slotSaveBBox(void);
  void slotLoadBBox(void);
  void slotSaveImage(void);
  void slotVoteMapImage(void);
  void slotInitializeSkeleton(void);
  void slotExtractSkeleton(void);
  void slotExtractSkeletonStep(void);
  void slotIncrementH0(void);
  void slotSwitchViewMode(void);
  void slotUpdateMapImages(void);
  void slotToggleRenderTrajPoints(bool);
  void slotToggleRenderMapImage(bool);
  void slotComputeSkeletonGraph(void);
  void slotClusterJunctionTrajectory(void);
  void slotCompleteTrajectories(void);
  void slotToggleCluster(QListWidgetItem *);
  void slotUpdateHeadingPlot(int ,float,float);
  void slotShowAllClusters();
  void slotHideAllClusters();
  void slotFilterClusterBySize(int minSize);
  void slotVisualizeSigma();
  void slotVisualizeDensity();
  void slotToggleTrajHeading(int checkStatus);
  void slotShowDenseTraj(int trajId);
  void slotFilterJunctionEdge();
 // void slotSaveSelectedClusters(void);
 signals:
  void existBBoxesToBeDeleted();
  void selectListToBeUpdated(int childId);
  void pointSelected(int idx,float x,float y);
  void nTrajToBeUpdated(int);
 protected:
  virtual void prepareContextMenu(QMenu* menu);
  virtual void mousePressEvent(QMouseEvent* event);
  virtual void mouseMoveEvent(QMouseEvent* event);
  virtual void mouseReleaseEvent(QMouseEvent* event);
  virtual void wheelEvent(QWheelEvent * event);

 private slots:
  void slotDeleteBBoxes(void);

//  void slotUpdateSelectList(int childId);

 private:
  Q_DISABLE_COPY(SceneWidget)

  osg::ref_ptr<GPSTrajectories>               gps_trajectories_;
  osg::Matrix                                 manipulator_matrix_;
  osg::Matrix                                 camera_matrix_;
  osg::ref_ptr<BBox2D>                        active_box_;
  std::vector<osg::ref_ptr<BBox2D>>           bboxes_to_be_deleted_;


  void updateInformation(void);
  bool isIn2DViewMode(void);
  osg::Vec3 screenToWorld(QPointF point_screen);
  void deleteAllBBoxes(void);
  void zoomin(void);
  void zoomout(void);
  void pan(int delta_x, int delta_y);
  int mouse_x_prev_, mouse_y_prev_;
  float distance_eye_;
  float look_at_x_;
  float look_at_y_;

};

#endif // SCENE_WIDGET_H_
