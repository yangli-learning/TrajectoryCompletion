#include <iostream>

#include <QToolTip>
#include <QKeyEvent>
#include <QSettings>
#include <QGridLayout>
#include <QDockWidget>
#include <QFileDialog>
#include <QMessageBox>
#include <QTextStream>
#include <QApplication>

#include <boost/filesystem.hpp>

#include "config.h"
#include "singleton.h"
#include "parameters.h"
#include "scene_widget.h"
#include "task_dispatcher.h"
#include "gps_trajectories.h"
#include "bbox_2d.h"
#include "main_window.h"

MainWindow* MainWindow::getInstance(void) {
  // The MainWindow can not be deleted by the normal Singleton in a crash free way.
  // Instead, setAttribute(Qt::WA_DeleteOnClose) for it, and make it deleted there.
  return SingletonWithNullDeleter<MainWindow>::instance();
}

MainWindow::MainWindow(QWidget * parent, Qt::WindowFlags flags)
  :QMainWindow(parent, flags),
   workspace_("."),
   scene_widget_(nullptr),
   allow_close_(true) {
  ui_.setupUi(this);
  installEventFilter(this);
  setAttribute(Qt::WA_DeleteOnClose);

  loadSettings();

  ui_.dockWidgetParameterEditor->hide();
  ui_.toolBar->addAction(ui_.dockWidgetParameterEditor->toggleViewAction());
  connect(this, SIGNAL(parametersLoaded()), this, SLOT(slotUpdateParameterEditor()));
  connect(ui_.pushButtonUpdateParameters, SIGNAL(pressed()), this, SLOT(slotSaveParameters()));
  loadParameters();

  scene_widget_ = new SceneWidget(this);

  setCentralWidget(scene_widget_);
  scene_widget_->startRendering();

  heading_plot_ = new QCustomPlot(this);//ui_.verticalLayout1);
  heading_plot_->setMinimumHeight(250);
  heading_plot_->xAxis2->setRange(-1,1);
  heading_plot_->xAxis->setRange(-1,1);
  heading_plot_->yAxis->setRange(-1,1);

  QLabel *heading_label = new QLabel(this);
  heading_label->setText("Heading Angle Visualization");
  ui_.verticalLayout1->addWidget(heading_label);
  ui_.verticalLayout1->addWidget(heading_plot_);
  heading_plot_->yAxis->setRangeReversed(true);
  heading_plot_->xAxis->setVisible(false);
  heading_plot_->xAxis2->setVisible(true);
  heading_plot_->show();

  //QLabel *sigma_label = new QLabel(this);
  //sigma_label->setText("Sigma at selected point: ");
  //ui_.verticalLayout1->addWidget(sigma_label);

  connect(this, SIGNAL(showInformationRequested(const QString&)), this, SLOT(slotShowInformation(const QString&)));
  connect(this, SIGNAL(showStatusRequested(const QString&, int)), this, SLOT(slotShowStatus(const QString&, int)));
  connect(scene_widget_,SIGNAL(selectListToBeUpdated(int)), this,SLOT(slotUpdateSelectList(int)) );
  // File
  connect(ui_.actionSetWorkspace, SIGNAL(triggered()), this, SLOT(slotSetWorkspace()));
  connect(ui_.actionReloadParameters, SIGNAL(triggered()), this, SLOT(slotLoadParameters()));
  connect(ui_.actionSaveBBox, SIGNAL(triggered()), scene_widget_, SLOT(slotSaveBBox()));
  connect(ui_.actionLoadBBox, SIGNAL(triggered()), scene_widget_, SLOT(slotLoadBBox()));
  connect(ui_.actionSaveImage,SIGNAL(triggered()),scene_widget_,SLOT(slotSaveImage()));
  // Algorithm
  connect(ui_.actionVoteMapImage, SIGNAL(triggered()), scene_widget_, SLOT(slotVoteMapImage()));
  connect(ui_.actionInitializeSkeleton, SIGNAL(triggered()), scene_widget_, SLOT(slotInitializeSkeleton()));
  connect(ui_.actionExtractSkeleton, SIGNAL(triggered()), scene_widget_, SLOT(slotExtractSkeleton()));
  connect(ui_.actionExtractSkeletonStep, SIGNAL(triggered()), scene_widget_, SLOT(slotExtractSkeletonStep()));
  connect(ui_.actionIncrementH0, SIGNAL(triggered()),scene_widget_,SLOT(slotIncrementH0()));
  connect(ui_.actionComputeSkeletonGraph, SIGNAL(triggered()), scene_widget_, 
          SLOT(slotComputeSkeletonGraph()));
  
  connect(ui_.actionClusterJunctionTrajectory,SIGNAL(triggered()),scene_widget_,SLOT(slotClusterJunctionTrajectory()));

  connect(ui_.actionCompleteTrajectories,SIGNAL(triggered()), scene_widget_,SLOT(slotCompleteTrajectories()));
  
  connect(ui_.clusterSelector, SIGNAL(itemClicked(QListWidgetItem*)),scene_widget_,
          SLOT(slotToggleCluster(QListWidgetItem*)));
  connect(ui_.actionSaveSelectedClusters,SIGNAL(triggered()),
          this,SLOT(slotSaveSelectedClusters()));
  connect(ui_.buttonSaveCluster, SIGNAL(clicked()),this,SLOT(slotSaveSelectedClusters()));
  connect(ui_.buttonHideAllClusters,SIGNAL(clicked()),scene_widget_,SLOT(slotHideAllClusters()));
  connect(ui_.buttonShowAllClusters, SIGNAL(clicked()),scene_widget_,SLOT(slotShowAllClusters()));
  // sigma label and visualization
  connect(ui_.buttonDisplaySigma, SIGNAL(clicked()),scene_widget_,SLOT(slotVisualizeSigma()));
  connect(this, SIGNAL(sigmaUpdated(float)),this,SLOT(slotUpdateSigmaLabel(float)));
  //visualize density
  connect(ui_.buttonDisplayDensity,SIGNAL(clicked()),scene_widget_,SLOT(slotVisualizeDensity()));
  connect(ui_.checkBoxTrajHeading,SIGNAL(stateChanged(int)),scene_widget_,SLOT(slotToggleTrajHeading(int)));
  // slider
  connect(ui_.sliderClusterSize,SIGNAL(sliderMoved(int)),this,SLOT(slotUpdateSliderLabel(int)));
  connect(ui_.sliderClusterSize,SIGNAL(sliderMoved(int)),scene_widget_,SLOT(slotFilterClusterBySize(int)));

  // Visualization
  connect(ui_.actionSwitchViewMode, SIGNAL(triggered()), scene_widget_, SLOT(slotSwitchViewMode()));
  connect(ui_.actionToggleTrajPoints, SIGNAL(toggled(bool)), scene_widget_, SLOT(slotToggleRenderTrajPoints(bool)));
  connect(ui_.actionToggleMapImage, SIGNAL(toggled(bool)), scene_widget_, SLOT(slotToggleRenderMapImage(bool)));
 connect(ui_.buttonShowDenseTraj,SIGNAL(clicked()),this,SLOT(slotCompleteSelectedTraj()));
 connect(scene_widget_,SIGNAL(nTrajToBeUpdated(int)),this,SLOT(slotUpdateTrajSpinBox(int)));

 //filter junction edge by flow
 connect(ui_.actionFilterJunctionEdge,SIGNAL(triggered()), scene_widget_,SLOT(slotFilterJunctionEdge()));
  // option pannel
  

  return;
}
MainWindow::~MainWindow() {
	// delete scene_widget_;
	//	delete heading_plot_;
}

void MainWindow::slotShowInformation(const QString& information) {
  QToolTip::showText(QCursor::pos(), information);
}

void MainWindow::showInformation(const std::string& information) {
  emit showInformationRequested(information.c_str());
}

void MainWindow::slotShowStatus(const QString& status, int timeout) {
  ui_.statusBar->showMessage(status, timeout);
}

void MainWindow::showStatus(const std::string& status, int timeout) {
	ui_.statusBar->clearMessage();
  emit showStatusRequested(status.c_str(), timeout);
}

void MainWindow::setAllowClose(bool allow_close) {
  allow_close_ = allow_close;

  return;
}

bool MainWindow::eventFilter(QObject *watched, QEvent *event) {
  if (!allow_close_) {
    if (event->type() == QEvent::Close) {
      event->ignore();
      showMinimized();
      return true;
    }
  }

  return QMainWindow::eventFilter(watched, event);
}

void MainWindow::closeEvent(QCloseEvent *event) {
  saveSettings();

  QMainWindow::closeEvent(event);

  return;
}

void MainWindow::keyPressEvent(QKeyEvent* event) {
  if (event->key() == Qt::Key_Down) {
    emit keyDownPressed();
  }

  QMainWindow::keyPressEvent(event);

  return;
}

void MainWindow::loadParameters(void) {
  std::string filename(PARAMETER_FILENAME);
  if (!boost::filesystem::exists(filename)) {
    std::cout << "Can not locate parameter file for GPStudio!" << std::endl;
    std::cout << "Expected location: \"" << filename << "\"" << std::endl;
    return;
  }

  Parameters::getInstance()->setParameterFile(filename);
  slotLoadParameters();

  return;
}

void MainWindow::slotUpdateParameterEditor(void) {
  std::string filename = Parameters::getInstance()->getFilename();
  if (filename.empty()) {
    ui_.plainTextEditParameterEditor->setPlainText("Error: parameter file is not set!");
    return;
  }

  QFile file(filename.c_str());
  file.open(QIODevice::ReadOnly);
  QTextStream text(&file);
  ui_.plainTextEditParameterEditor->setPlainText(text.readAll());

  return;
}

void MainWindow::slotSaveParameters(void) {
  std::string filename = Parameters::getInstance()->getFilename();
  if (filename.empty()) {
    ui_.plainTextEditParameterEditor->setPlainText("Error: parameter file is not set!\n"
        +ui_.plainTextEditParameterEditor->toPlainText());
    return;    
  }

  QFile file(filename.c_str());
  file.open(QIODevice::WriteOnly);
  QTextStream text(&file);
  text << ui_.plainTextEditParameterEditor->toPlainText();
  file.close();

  ui_.statusBar->showMessage("Parameter Saved!");

  float map_image_resolution = Parameters::getInstance()->map_image_resolution;
  slotLoadParameters(false);
  if(map_image_resolution != Parameters::getInstance()->map_image_resolution) {
    scene_widget_->slotUpdateMapImages();
  }

  return;
}

void MainWindow::slotLoadParameters(bool show_message) {
  Parameters::getInstance()->readParameters();
  emit parametersLoaded();

  if (show_message)
    ui_.statusBar->showMessage("Parameter Loaded!");

  return;
}

bool MainWindow::slotSetWorkspace(void) {
  QString directory = QFileDialog::getExistingDirectory(this, tr("Set Workspace"), workspace_.c_str(),
                      QFileDialog::ShowDirsOnly);

  if (directory.isEmpty())
    return false;

  workspace_ = directory.toStdString();

  return true;
}

void MainWindow::loadSettings(void) {
  QSettings settings("GPStudio", "GPStudio");

  workspace_ = settings.value("workspace", ".").toString().toStdString();

  return;
}

void MainWindow::saveSettings(void) {
  QSettings settings("GPStudio", "GPStudio");

  QString workspace(workspace_.c_str());
  settings.setValue("workspace", workspace);

  return;
}

bool MainWindow::slotShowYesNoMessageBox(const std::string& text, const std::string& informative_text) {
  QMessageBox msg_box;
  msg_box.setText(text.c_str());
  msg_box.setInformativeText(informative_text.c_str());
  msg_box.setStandardButtons(QMessageBox::Yes | QMessageBox::No);
  msg_box.setDefaultButton(QMessageBox::Yes);
  int ret = msg_box.exec();

  return (ret == QMessageBox::Yes);
}

void MainWindow::slotCleanUpTaskDispatcher(void) {
  TaskDispatcher::getInstance()->cleanUp();
  return;
}


void MainWindow::slotUpdateSelectList(int childId){
	qDebug()<< "updaing cluster select list...";
    osg::Node* child = scene_widget_->getSceneChild(childId);
	BBox2D* bbox = dynamic_cast<BBox2D*>(child);
	if (bbox != nullptr){
        
		JunctionNetwork* network  = bbox->getJunctionNetwork();
		ui_.clusterSelector->show();
		ui_.clusterSelector->clear();
		if (network!=nullptr){
			int nClusters = network->getNumClusters();
            int maxSize = network->getMaxClusterSize();
			qDebug() << "adding " <<nClusters <<" items...";
            QStringList strlist =  network->clusterToQStringList();
			ui_.clusterSelector->addItems(strlist);
			for (int k=0; k<ui_.clusterSelector->count();++k){
				ui_.clusterSelector->item(k)->setCheckState( Qt::Unchecked);
			}
            ui_.sliderClusterSize->setRange(1,maxSize);

			qDebug() <<"read " << ui_.clusterSelector->count() << " rows ";
			qDebug() <<"currently selected " << ui_.clusterSelector->selectedItems().size() << " items";
		}else{
			qDebug()<<"Error: Junction network does not exist!";
		}
        
	}
}

void MainWindow::slotUpdateSigmaLabel(float value){
    QString txt = QString("Sigma of selected point: %1").arg(value);
    ui_.sigmaLabel->setText(txt);
}

void MainWindow::slotUpdateSliderLabel(int value){
	QString txt = QString("Display clusters of minimum size %1").arg(value);
	ui_.clusterFilterLabel->setText(txt);
}

void MainWindow::slotSaveSelectedClusters(void){
    std::vector<bool> selected(ui_.clusterSelector->count(),false) ;
    int ct=0;
    for (int k=0; k<ui_.clusterSelector->count();++k){
        QListWidgetItem *item = ui_.clusterSelector->item(k);
        if (item->checkState() == Qt::Checked){
            selected[k] = true;ct++;
        }
    }
    qDebug() <<"saving " << ct <<
               " out of "  << selected.size() << " selected clusters ..";
    scene_widget_->saveSelectedClusters(selected);

}


void MainWindow::slotCompleteSelectedTraj(void){
	int  trajId = ui_.trajIdSpinBox->value();
	scene_widget_->slotShowDenseTraj(trajId);
}

void MainWindow::slotUpdateTrajSpinBox(int max){
	ui_.trajIdSpinBox->setMaximum(max);
}
