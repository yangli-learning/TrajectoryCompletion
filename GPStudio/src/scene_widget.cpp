#define _USE_MATH_DEFINES
#include <cmath>

#include <thread>
#include <string>
#include <sstream>

#include <QMenu>
#include <QFile>
#include <QFileDialog>
#include <QMetaObject>
#include <QElapsedTimer>
#include <QDebug>

#include <osg/LightModel>
#include <osgGA/TrackballManipulator>

#include <boost/filesystem/path.hpp>

#include "common.h"
#include "bbox_2d.h"
#include "color_map.h"
#include "parameters.h"
#include "information.h"
#include "main_window.h"
#include "pick_handler.h"
#include "gps_trajectories.h"
#include "l_shape_iterator.h"
#include "scene_widget.h"
#include "l1skeleton.h" 
#include <QInputDialog>
SceneWidget::SceneWidget(QWidget * parent, 
						 const QGLWidget * shareWidget, Qt::WindowFlags f) :
    OSGViewerWidget(parent, shareWidget, f), 
	gps_trajectories_(new GPSTrajectories()), active_box_(new BBox2D(gps_trajectories_)) {
  osg::ref_ptr<osg::LightModel> light_model = new osg::LightModel();
  light_model->setTwoSided(true);
  scene_root_->getOrCreateStateSet()->setAttributeAndModes(light_model.get(), osg::StateAttribute::ON);
  active_box_->setSceneWidget(this);
  addSceneChild(gps_trajectories_);
  addSceneChild(active_box_);
  addEventHandler(new PickHandler());

  centerScene();

  connect(this, SIGNAL(existBBoxesToBeDeleted()), this, SLOT(slotDeleteBBoxes()));
  connect(this, SIGNAL(pointSelected(int,float,float)),
		  this,SLOT(slotUpdateHeadingPlot(int,float,float)));
  //connect(this, SIGNAL(pointSelected(int,float,float)),
   //       parent, SLOT(slotUpdateSigmaLabel(in)))
  //connect(this, SIGNAL(selectListToBeUpdated(int)), this,SLOT(slotUpdateSelectedList(int)));
  return;
}

SceneWidget::~SceneWidget() {

}

void SceneWidget::centerScene(void) {
  std::lock_guard < std::mutex > lock(mutex_);

  osg::BoundingBox bounding_box = getBoundingBox();
  if (bounding_box.xMax() == -std::numeric_limits<float>::max()) {
    bounding_box._min = osg::Vec3(-1.0, -1.0, -1.0);
    bounding_box._max = osg::Vec3(1.0, 1.0, 1.0);
  }

  osg::Vec3 center = (bounding_box._min + bounding_box._max) / 2;
  double distance = (bounding_box._max - bounding_box._min).length();
  osg::Vec3 eye_offset = osg::Vec3(0, 0, 1) * distance;

  osgGA::CameraManipulator* camera_manipulator = getCameraManipulator();
  if (camera_manipulator != nullptr) {
    camera_manipulator->setHomePosition(center + eye_offset, center, up_vector_);
    camera_manipulator->home(0);
  }

  return;
}

void SceneWidget::centerScene2D(void) {
  std::lock_guard < std::mutex > lock(mutex_);

  osg::BoundingBox bounding_box = getBoundingBox();
  if (bounding_box.xMax() == -std::numeric_limits<float>::max()) {
    bounding_box._min = osg::Vec3(-1.0, -1.0, -1.0);
    bounding_box._max = osg::Vec3(1.0, 1.0, 1.0);
  }
  osg::Vec3 center = (bounding_box._min + bounding_box._max) / 2;

  double fovy, aspect_ratio, z_near, z_far;
  getCamera()->getProjectionMatrixAsPerspective(fovy, aspect_ratio, z_near, z_far);
  fovy = fovy * M_PI / 180;
  double fovx = 2 * std::atan(aspect_ratio * std::tan(fovy / 2));

  double distancey = (gps_trajectories_->getHeight() / 2) / std::tan(fovy / 2);
  double distancex = (gps_trajectories_->getWidth() / 2) / std::tan(fovx / 2);
  distance_eye_ = std::max(distancey, distancex);
  look_at_x_ = center.x();
  look_at_y_ = center.y();
  osg::Vec3 eye_offset = osg::Vec3(0, 0, 1) * distance_eye_;
  getCamera()->setViewMatrixAsLookAt(center + eye_offset, center, up_vector_);

  return;
}

void SceneWidget::zoomout(void) {
  std::lock_guard < std::mutex > lock(mutex_);

  osg::Vec3 center(look_at_x_, look_at_y_, 0.0f);
  distance_eye_ *= 1.1;
  osg::Vec3 eye_offset = osg::Vec3(0, 0, 1) * distance_eye_;
  getCamera()->setViewMatrixAsLookAt(center + eye_offset, center, up_vector_);

  return;
}

void SceneWidget::zoomin(void) {
  std::lock_guard < std::mutex > lock(mutex_);

  osg::Vec3 center(look_at_x_, look_at_y_, 0.0f);
  distance_eye_ /= 1.1;
  osg::Vec3 eye_offset = osg::Vec3(0, 0, 1) * distance_eye_;
  getCamera()->setViewMatrixAsLookAt(center + eye_offset, center, up_vector_);

  return;
}

void SceneWidget::pan(int delta_x, int delta_y) {
  std::lock_guard < std::mutex > lock(mutex_);

  double fovy, aspect_ratio, z_near, z_far;
  getCamera()->getProjectionMatrixAsPerspective(fovy, aspect_ratio, z_near, z_far);
  fovy = fovy * M_PI / 180;
  double fovx = 2 * std::atan(aspect_ratio * std::tan(fovy / 2));

  float height = 2 * distance_eye_ * tan(fovy / 2);
  float width = 2 * distance_eye_ * tan(fovx / 2);
  look_at_x_ += delta_x * (width / this->width());
  look_at_y_ += delta_y * (height / this->height());
  osg::Vec3 center(look_at_x_, look_at_y_, 0.0f);
  osg::Vec3 eye_offset = osg::Vec3(0, 0, 1) * distance_eye_;
  getCamera()->setViewMatrixAsLookAt(center + eye_offset, center, up_vector_);

  return;
}

void SceneWidget::prepareContextMenu(QMenu* menu) {
  menu->addAction("Open GPS Trajectories", this, SLOT(slotOpenGPSTrajectories()));
  //QMenu* scene_widget_menu = menu->addMenu("Scene Widget");
  //menu->addMenu(scene_widget_menu);

  OSGViewerWidget::prepareContextMenu(menu);

  return;
}

osg::Vec3 SceneWidget::screenToWorld(QPointF point_screen) {
  osg::Matrix VPW = getCamera()->getViewMatrix() * getCamera()->getProjectionMatrix() * getCamera()->getViewport()->computeWindowMatrix();
  osg::Matrix VPW_inverse;
  VPW_inverse.invert(VPW);

  osg::Vec3 p1 = osg::Vec3(point_screen.x(), (height() - point_screen.y()), 0.0f) * VPW_inverse;
  osg::Vec3 p2 = osg::Vec3(point_screen.x(), (height() - point_screen.y()), 1.0f) * VPW_inverse;

  return osg::Vec3((p1.x() * p2.z() - p2.x() * p1.z()) / (p2.z() - p1.z()), (p1.y() * p2.z() - p2.y() * p1.z()) / (p2.z() - p1.z()), 0.0f);
}

void SceneWidget::mousePressEvent(QMouseEvent* event) {
  if (isIn2DViewMode() && event->modifiers() == Qt::ControlModifier) {
    active_box_->addCorner(screenToWorld(event->localPos()));
  }

  if (isIn2DViewMode() && event->button() == Qt::MidButton) {
    mouse_x_prev_ = event->x();
    mouse_y_prev_ = event->y();
  }

  OSGViewerWidget::mousePressEvent(event);
  return;
}
void SceneWidget::mouseMoveEvent(QMouseEvent* event) {
  if (isIn2DViewMode() && event->modifiers() == Qt::ControlModifier && active_box_->initialized()) {
    active_box_->addCorner(screenToWorld(event->localPos()));
  }

  if (isIn2DViewMode() && event->buttons() == Qt::MidButton) {
    int delta_x = event->x() - mouse_x_prev_;
    int delta_y = event->y() - mouse_y_prev_;
    pan(-delta_x, delta_y);
    mouse_x_prev_ = event->x();
    mouse_y_prev_ = event->y();
  }

  OSGViewerWidget::mouseMoveEvent(event);
  return;
}
void SceneWidget::mouseReleaseEvent(QMouseEvent* event) {
  if (isIn2DViewMode() && event->modifiers() == Qt::ControlModifier) {
    if (active_box_->isDegenerated()) {
      active_box_->clearCorners();
    } else {
      active_box_->setBuilt(true);
      processBBox(active_box_);
      qDebug()<< "++ finished calling process bbox...";
      active_box_ = new BBox2D(gps_trajectories_);
      qDebug()<< "++ finished calling bbox constructor...";
	  active_box_->setSceneWidget(this);
	  addSceneChild(active_box_);
    }
  }

  OSGViewerWidget::mouseReleaseEvent(event);
  return;
}

void SceneWidget::wheelEvent(QWheelEvent * event) {
  if (isIn2DViewMode()) {
    if (event->angleDelta().y() > 0)
      zoomout();
    else if (event->angleDelta().y() < 0)
      zoomin();
  }
  OSGViewerWidget::wheelEvent(event);
  return;
}


void SceneWidget::slotInitializeSkeleton(void){
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr)
        bbox->initializeSkeleton();
    }
  }	
  return;
}
void SceneWidget::slotIncrementH0(void){
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr)
        bbox->incrementH0();
    }
  }	
  return;
}

void SceneWidget::slotShowDenseTraj(int trajId){
    int num_children = scene_root_->getNumChildren();
    for (int i = 0; i < num_children; ++i) {
      osg::Node* child = scene_root_->getChild(i);
      if (child != active_box_) {
        BBox2D* bbox = dynamic_cast<BBox2D*>(child);
        if (bbox != nullptr)
          bbox->completeSingleTrajectory(trajId);
      }
    }
    return;
}
 
void SceneWidget::slotUpdateHeadingPlot(int idx,float px,float py){
	qDebug()<< "updating heading plot at " <<px <<","<<py;
	std::vector<float> head;
    cv::Point2f principalDir; float h =0.0;
    cv::Mat cov;
    L1Skeleton_* sk =0; //L1Skeleton* l1skeleton=0;
    bool useHeading = Parameters::getInstance()->skeleton_use_heading;
	int num_children = scene_root_->getNumChildren();
	for (int i = 0; i < num_children; ++i) {
		osg::Node* child = scene_root_->getChild(i);
		if (child != active_box_) {
			BBox2D* bbox = dynamic_cast<BBox2D*>(child);

			if (bbox != nullptr){
                VotingMapRenderer *voting_map =bbox->getVotingMap();
                head= voting_map->getHeadingVectorAtPos(px,py);
                 sk  = bbox->getSkeleton();
				 principalDir = sk->getPrincipalDirection(idx);
				 h = sk->getAdaptiveH(px,py);
                 //l1skeleton = dynamic_cast<L1Skeleton*>(sk);
                // if ( l1skeleton){
                cov = sk->getSkeletonCovariance(idx);
					 //					 cov = cv::Mat(dynamic_cast<L1Skeleton*>(sk)->
                 //}
                std::ostringstream oss;
                oss<<cov;
                qDebug()<<"Covariance matrix =" << oss.str().c_str();
				 break;
			}
		}
	}	
	if (sk==nullptr){
		return;
	}
	int nDirections = head.size();
    MainWindow* main_window = MainWindow::getInstance();
    emit main_window->sigmaUpdated(sk->getSkeletonLinearity(idx)); //update sigma label

	QCustomPlot* heading_plot = main_window->getHeadingPlot();
	heading_plot->legend->setVisible(true);
    heading_plot->legend->setIconSize(5,2);
	int nCleared = heading_plot->clearPlottables();
	nCleared+= heading_plot->clearItems();
	heading_plot->clearGraphs();

    // make color map of gaussian /multivariate gaussian weight ------------
    QCPColorMap *cmap=new QCPColorMap(heading_plot->xAxis2,heading_plot->yAxis);
    cmap->data()->setSize(50,50);
    cmap->data()->setRange(QCPRange(-1,1), QCPRange(-1, 1));
    float xVal,yVal,zVal;
    for (int x=0; x<50; ++x){
        for (int y=0; y<50; ++y){
            xVal = h*(x-25.0)/25.f ; yVal = h*(y-25.0)/25.f ;
            if (Parameters::getInstance()->skeleton_mvn_kernel){
                if ( cov.cols>0 && cov.rows>0){
                   // zVal = sk->thetaMVN(cv::Point2f(xVal, yVal),
                    //                                           h, cov);
                    zVal = sk->thetaMVN(cv::Point2f(xVal, yVal),h,cov);
                    cmap->data()->setCell(x, y, zVal);
                }
            }else{
                zVal = sk->theta(sqrt(xVal*xVal+yVal*yVal),h);
                cmap->data()->setCell(x, y, zVal);
            }
        }
    }
    cmap->setGradient(QCPColorGradient::gpJet);
    cmap->rescaleDataRange(true);
    cmap->setName(QString("h=%1").arg(h) );
     QCPRange range = cmap->dataRange();
    range.lower=-0.2;
     cmap->setDataRange(range);
    heading_plot->addPlottable(cmap);

    if (nDirections ==2){
		QCPItemLine *arrow1 = new QCPItemLine(heading_plot);
		//	arrow1->setName("Average heading vector");
		arrow1->setPen(QPen(Qt::blue));
		heading_plot->addItem(arrow1);
		arrow1->start->setCoords(0,0);
		arrow1->end->setCoords(head[0], head[1]);
		arrow1->setHead(QCPLineEnding::esSpikeArrow);

		heading_plot->addGraph();
		heading_plot->graph(0)->setName("Average heading");
		heading_plot->graph(0)->setPen(QPen(Qt::blue));
		heading_plot->addGraph();
		heading_plot->graph(1)->setName("PCA direction");
		heading_plot->graph(1)->setPen(QPen(Qt::black));

	}else if(nDirections >2){
		QVector<double> t(nDirections+1),x(nDirections+1),y(nDirections+1),
			x0(nDirections+1),y0(nDirections+1);
		int i;
         //normalize head vector for display
        float maxVote = (Parameters::getInstance()->skeleton_normalize_heading )?
                    1.f: *(std::max_element(head.begin(),head.end()));

		for (i=0;i<nDirections;++i){
            cv::Point2f p = LShapeIterator::UnitDirection(( nDirections-i)*float(360/nDirections));
            if (useHeading){
                if (i<nDirections/2){
                    x[i] = p.x *(head[i]+head[i+int(nDirections/2)])/(2*maxVote);
                    y[i] = p.y * (head[i]+head[i+ int(nDirections/2)])/(2*maxVote) ;
                }else {
                    x[i] = p.x *(head[i]+head[i%int(nDirections/2)])/(2*maxVote);
                    y[i] = p.y *(head[i]+head[i%int(nDirections/2)])/(2*maxVote) ;
                }
            }else{
                x[i] =p.x*head[i];
                y[i]=p.y*head[i];
            }
            x0[i] = p.x; // boundary x
            y0[i] = p.y; // boundary y
			t[i]= i;
		}

		t[i] = nDirections;
		cv::Point2f p =  LShapeIterator::UnitDirection(0);
        x[i] = x[0];
        y[i] = y[0];
		x0[i]=p.x;
		y0[i]=p.y;

        QCPCurve *g =  new QCPCurve(heading_plot->xAxis2,heading_plot->yAxis);
		g->setBrush(QBrush(QColor(0, 0, 255, 20)));
        g->setName(QString("(%1,%2)").arg(px).arg(py));
		heading_plot->addPlottable(g);

        QCPCurve *g0 = new QCPCurve(heading_plot->xAxis2,heading_plot->yAxis);
        g0->setPen(QPen(Qt::red));
        g->setData(t,x,y);
        g0->setData(t,x0,y0);
        heading_plot->addPlottable(g0);
	}

    connect(heading_plot->xAxis2, SIGNAL(rangeChanged(QCPRange)),
            heading_plot->xAxis, SLOT(setRange(QCPRange)));
    connect(heading_plot->yAxis, SIGNAL(rangeChanged(QCPRange)),
            heading_plot->yAxis2, SLOT(setRange(QCPRange)));


	QCPItemLine *arrow = new QCPItemLine(heading_plot);
	//	arrow->setName("PCA direction");
	heading_plot->addItem(arrow);
	arrow->start->setCoords(0,0);
	arrow->end->setCoords(principalDir.x, principalDir.y);
	arrow->setHead(QCPLineEnding::esSpikeArrow);

	//	g->rescaleAxes();
	heading_plot->replot();
	heading_plot->setInteractions(QCP::iRangeDrag | QCP::iRangeZoom | QCP::iSelectPlottables);
}
void SceneWidget::slotExtractSkeleton(void){
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr){
        bbox->extractSkeleton();
      }
    }
  }	
  return;
}

void SceneWidget::slotComputeSkeletonGraph(void){
   
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr)
        bbox->computeSkeletonGraph();
    }
  }	
  return;	
	
}
void SceneWidget::slotClusterJunctionTrajectory(void){

	int num_children = scene_root_->getNumChildren();
	for (int i = 0; i < num_children; ++i) {
		osg::Node* child = scene_root_->getChild(i);
		if (child != active_box_) {
			BBox2D* bbox = dynamic_cast<BBox2D*>(child);
            if (bbox != nullptr){
				bbox->clusterJunctionTrajectory();
                emit selectListToBeUpdated(i);
            }
		}
	}	
	return;
}
void SceneWidget::slotCompleteTrajectories(void){

	int num_children = scene_root_->getNumChildren();
	for (int i = 0; i < num_children; ++i) {
		osg::Node* child = scene_root_->getChild(i);
		if (child != active_box_) {
			BBox2D* bbox = dynamic_cast<BBox2D*>(child);
            if (bbox != nullptr){
				bbox->completeTrajectories();
				//               emit selectListToBeUpdated(i);
            }
		}
	}	
	return;
}

void SceneWidget::slotExtractSkeletonStep(void){
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr)
        bbox->extractSkeletonStep();
    }
  }	
  return;
}

void SceneWidget::slotVoteMapImage(void) {
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr){
        VotingMapRenderer *voting_map = bbox->getVotingMap();
        voting_map->voteMapImage();

		std::ostringstream oss;
		oss <<"votingImg" << i<<".png";
        voting_map->exportVotingImage(oss.str());

		oss.str("");
		oss << "heading" << i ;
        voting_map->exportHeadingField(oss.str());
	  }
    }
  }

  return;
}

void SceneWidget::slotDeleteBBoxes(void) {
  while (!bboxes_to_be_deleted_.empty()) {
    removeSceneChild(bboxes_to_be_deleted_.back());
    bboxes_to_be_deleted_.pop_back();
  }
  return;
}
void SceneWidget::slotToggleCluster(QListWidgetItem* item){
    QListWidget * lw = item->listWidget();
    int clusterId = lw->row(item);

  //  qDebug() <<"show cluster " <<clusterId;
    int num_children = scene_root_->getNumChildren();
    for (int i = 0; i < num_children; ++i) {
      osg::Node* child = scene_root_->getChild(i);
      if (child != active_box_) {
        BBox2D* bbox = dynamic_cast<BBox2D*>(child);
        if (bbox != nullptr){
          item->setCheckState( (item->checkState()==Qt::Checked)?
                                        Qt::Unchecked: Qt::Checked);
          bbox->toggleTrajectoryCluster(clusterId,
										item->checkState()==Qt::Checked );//voteMapImageLShape();
        }
      }
    }
    return;
}

void SceneWidget::slotToggleTrajHeading(int checkStatus){

    bool checked = checkStatus == Qt::Checked;
    qDebug() << "Display trajectory heading " << checked;
     MainWindow* main_window = MainWindow::getInstance();
     QListWidget *selector = main_window->getClusterSelector();
     std::vector<bool> checkedIds;
     for (int i=0; i< selector->count(); ++i){
         if (selector->item(i)->checkState() == Qt::Checked){
             checkedIds.push_back(true);
         }else{
			    checkedIds.push_back(false);
		 }
     }

    int num_children = scene_root_->getNumChildren();
    for (int i = 0; i < num_children; ++i) {
      osg::Node* child = scene_root_->getChild(i);
      if (child != active_box_) {
        BBox2D* bbox = dynamic_cast<BBox2D*>(child);
        if (bbox != nullptr){
          bbox->displayTrajectoryHeading(checkedIds, checked);
        }
      }
    }

}

void SceneWidget::slotVisualizeSigma(){
    int num_children = scene_root_->getNumChildren();
    for (int i = 0; i < num_children; ++i) {
      osg::Node* child = scene_root_->getChild(i);
      if (child != active_box_) {
        BBox2D* bbox = dynamic_cast<BBox2D*>(child);
        if (bbox != nullptr){
          bbox->getSkeleton()->visualizeSigma();
        }
      }
    }

    return;
}
void SceneWidget::slotVisualizeDensity(){
    int num_children = scene_root_->getNumChildren();
    for (int i = 0; i < num_children; ++i) {
      osg::Node* child = scene_root_->getChild(i);
      if (child != active_box_) {
        BBox2D* bbox = dynamic_cast<BBox2D*>(child);
        if (bbox != nullptr){
          bbox->getSkeleton()->visualizeDensity();
        }
      }
    }
    return;
}

void SceneWidget::removeBBox(BBox2D* bbox) {
  bboxes_to_be_deleted_.push_back(bbox);
  emit existBBoxesToBeDeleted();

  return;
}

void SceneWidget::processBBox(BBox2D* bbox) {
  // the is executed after the bbox is drawn.
	//bbox->voteMapImageLShape();
	qDebug()<< "Created bounding box...";
  return;
}

void SceneWidget::updateInformation(void) {
  std::string information;
  information_->setText(information);

  if (information_->isHidden())
    information_->toggleHidden();

  return;
}

void SceneWidget::deleteAllBBoxes(void) {
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr)
        bboxes_to_be_deleted_.push_back(bbox);
    }
  }
  if (!bboxes_to_be_deleted_.empty())
    emit existBBoxesToBeDeleted();

  return;
}

void SceneWidget::slotOpenGPSTrajectories(void) {
  MainWindow* main_window = MainWindow::getInstance();

#ifdef _DEBUG
#define Q_FILE_DIALOG_WORKAROUND , nullptr, QFileDialog::DontUseNativeDialog
#else
#define Q_FILE_DIALOG_WORKAROUND
#endif
  QString filename = QFileDialog::getOpenFileName(main_window, "Open GPS Trajectories", main_window->getWorkspace().c_str(),
      "GPS Trajectories (*.pbf *.txt)" Q_FILE_DIALOG_WORKAROUND);
  if (filename.isEmpty())
    return;

  if (gps_trajectories_->load(filename.toStdString())) {
    centerScene();
    deleteAllBBoxes();
    if (!isIn2DViewMode()) {
      slotSwitchViewMode();
    } else {
      centerScene2D();
    }
	emit nTrajToBeUpdated(gps_trajectories_->size());
  }

  return;
}

bool SceneWidget::isIn2DViewMode(void) {
  return (getCameraManipulator() == nullptr);
}

void SceneWidget::slotSwitchViewMode(void) {
  if (isIn2DViewMode()) {
    setMouseTracking(false);

    setCameraManipulator(new osgGA::TrackballManipulator);
    getCameraManipulator()->setByMatrix(manipulator_matrix_);
    getCamera()->setViewMatrix(camera_matrix_);
  } else {
    setMouseTracking(true);

    manipulator_matrix_ = getCameraManipulator()->getMatrix();
    camera_matrix_ = getCamera()->getViewMatrix();
    setCameraManipulator(nullptr);

    centerScene2D();
  }
  return;
}

void SceneWidget::slotToggleRenderTrajPoints(bool checked) {
  gps_trajectories_->toggleRenderTrajPoints(checked);
  return;
}

void SceneWidget::slotToggleRenderMapImage(bool checked) {
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr)
        bbox->toggleRenderMapImage(checked);
    }
  }

  return;
}

void SceneWidget::slotUpdateMapImages(void) {
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr){
          VotingMapRenderer *voting_map=bbox->getVotingMap();
          voting_map->init();
      }
    }
  }
  return;
}
void SceneWidget::slotLoadBBox(void){
	MainWindow* main_window = MainWindow::getInstance();

#ifdef _DEBUG
#define Q_FILE_DIALOG_WORKAROUND , nullptr, QFileDialog::DontUseNativeDialog
#else
#define Q_FILE_DIALOG_WORKAROUND
#endif
	QString fname = 
		QFileDialog::getOpenFileName(this,tr("Bounding box base name"),
									 main_window->getWorkspace().c_str(),
									 tr("XML files (*.xml *.yml *.yaml)")
									 );
	
	if (isIn2DViewMode() ) {
		if (active_box_->importBBox(fname.toStdString()) ){
			active_box_ = new BBox2D(gps_trajectories_);
			active_box_->setSceneWidget(this);
			addSceneChild(active_box_);
		}

	}

	return;
}
void SceneWidget::slotSaveImage(){
	
	MainWindow* main_window = MainWindow::getInstance();


#ifdef _DEBUG
#define Q_FILE_DIALOG_WORKAROUND , nullptr, QFileDialog::DontUseNativeDialog
#else
#define Q_FILE_DIALOG_WORKAROUND
#endif
	QString filename = 
		QFileDialog::getSaveFileName(this,tr("Save image as"),
									 main_window->getWorkspace().c_str());
	//!!!!!!!!!!!!!!!!!!!
	/*

	int x,y,width,height;
	auto camera = getCamera();
	x = camera.getViewport()->x();
	y = camera.getViewport()->y();
	width = camera.getViewport()->width();
	height = camera.getViewport()->height();

	osg::ref_ptr<osg::Image> image = new osg::Image;
	image->readPixels(x,y,width,height,GL_RGB,GL_UNSIGNED_BYTE);
	if (osgDB::writeImageFile(*image, filename)){
		qDebug() <<"Saving screen as "<<filename<<"...";
	}else{
		qDebug() <<"Error saving screen to image file!";
		}*/
	return;
		}
void SceneWidget::slotSaveBBox(void) {
	MainWindow* main_window = MainWindow::getInstance();


#ifdef _DEBUG
#define Q_FILE_DIALOG_WORKAROUND , nullptr, QFileDialog::DontUseNativeDialog
#else
#define Q_FILE_DIALOG_WORKAROUND
#endif
	QString basename = 
		QFileDialog::getSaveFileName(this,tr("Bounding box basename"),
									 main_window->getWorkspace().c_str());
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
	  QString fname = basename+QString::number(i)+QString(".yaml");
      if (bbox != nullptr){
		  bbox->exportBBox(fname.toStdString());
		  main_window->showStatus("Loaded bounding box");
	  }
    }
  }
  return;
}
/*
void SceneWidget::toggleClusterDisplay(std::vector<bool> selected,bool noline){

}*/

void SceneWidget::saveSelectedClusters(std::vector<bool> selected){
    MainWindow* main_window = MainWindow::getInstance();
    qDebug() <<"scenewidget::saveSelectedClusters...";

    #ifdef _DEBUG
    #define Q_FILE_DIALOG_WORKAROUND , nullptr, QFileDialog::DontUseNativeDialog
    #else
    #define Q_FILE_DIALOG_WORKAROUND
    #endif

    QString basename=
        QFileDialog::getSaveFileName(this,tr("Cluster file basename"),
                                     main_window->getWorkspace().c_str());
    qDebug() << "basename is " <<basename;
    int num_children = scene_root_->getNumChildren();
    for (int i = 0; i < num_children; ++i) {
      osg::Node* child = scene_root_->getChild(i);
      if (child != active_box_) {
        BBox2D* bbox = dynamic_cast<BBox2D*>(child);
        QString fname = basename+QString::number(i)+QString(".dat");
        QString fname2 = basename+QString::number(i)+QString(".index.dat");
        if (bbox != nullptr){
            bbox->exportClusters(selected, fname.toStdString(),fname2.toStdString());
            main_window->showStatus("Saved cluster file");
        }
      }
    }
}

void SceneWidget::slotShowAllClusters(){
    MainWindow* main_window = MainWindow::getInstance();
	QListWidget *selector = main_window->getClusterSelector();
	
	for (int i=0; i< selector->count(); ++i){
	    selector->item(i)->setCheckState(Qt::Checked);
    }

	int num_children = scene_root_->getNumChildren();
	for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr)
		  bbox->displayTrajectoryClusters(true);
    }
  }
 
}

void SceneWidget::slotHideAllClusters(){
    MainWindow* main_window = MainWindow::getInstance();
	QListWidget *selector = main_window->getClusterSelector();
	
	for (int i=0; i< selector->count(); ++i){
	    selector->item(i)->setCheckState(Qt::Unchecked);
    }
  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr)
		  bbox->displayTrajectoryClusters(false);
    }
  }
  return;	
}

void SceneWidget::slotFilterClusterBySize(int minSize){
    MainWindow* main_window = MainWindow::getInstance();
    QListWidget *selector = main_window->getClusterSelector();


  int num_children = scene_root_->getNumChildren();
  for (int i = 0; i < num_children; ++i) {
    osg::Node* child = scene_root_->getChild(i);
    if (child != active_box_) {
      BBox2D* bbox = dynamic_cast<BBox2D*>(child);
      if (bbox != nullptr){
          std::vector<bool> row = bbox->filterClustersBySize(minSize);
          // update selector checkbox! selector->
          for (int r=0;r<selector->count();++r){
              selector->item(r)->setCheckState( row[r]?Qt::Checked:Qt::Unchecked);
          }

      }
    }
  }
  return;
}

void SceneWidget::slotFilterJunctionEdge(){
	//	MainWindow* main_window = MainWindow::getInstance();
	double thresh; bool ok;
	thresh =QInputDialog::getDouble(this, tr("Filter junction edge by flow size"), tr("minimum flow size"),
									10,0,500,0.5,&ok);
	if (ok)
		qDebug() <<"Filter junction edge by minimum flow size " << thresh;
	
    MainWindow* main_window = MainWindow::getInstance();
	
	int num_children = scene_root_->getNumChildren();
	for (int i = 0; i < num_children; ++i) {
		osg::Node* child = scene_root_->getChild(i);
		if (child != active_box_) {
			BBox2D* bbox = dynamic_cast<BBox2D*>(child);
			if (bbox != nullptr){
				main_window->slotShowStatus(QString("Filter junction edges with flow less than %1").arg(thresh),0);
				bbox->computeSkeletonGraph(thresh);
				bbox->clusterJunctionTrajectory();
                emit selectListToBeUpdated(i);
			}
		}
	}
}
