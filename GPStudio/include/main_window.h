#pragma once
#ifndef MainWindow_H
#define MainWindow_H

#include <mutex>
#include <string>
#include <cassert>

#include <QMainWindow>

#include "ui_main_window.h"
#include "qcustomplot.h"

class SceneWidget;
class QPlainTextEdit;

class MainWindow: public QMainWindow {
Q_OBJECT

public:
  MainWindow(QWidget * parent = 0, Qt::WindowFlags flags = 0);
  virtual ~MainWindow(); 

  static MainWindow* getInstance(void);

  const std::string& getWorkspace(void) const {
    return workspace_;
  }
  SceneWidget* getSceneWidget(void) {
    return scene_widget_;
  }

  void showInformation(const std::string& information);
  void showStatus(const std::string& status, int timeout = 0);

  void setAllowClose(bool allow_close);
  QCustomPlot* getHeadingPlot(){
    return heading_plot_;
  }
  QListWidget* getClusterSelector(){
	  return ui_.clusterSelector;
  }
public slots:
  bool slotShowYesNoMessageBox(const std::string& text, const std::string& informative_text);
  void slotCleanUpTaskDispatcher(void);
  void slotUpdateSelectList(int);
  void slotSaveSelectedClusters(void);
  void slotUpdateSliderLabel(int);
  void slotUpdateSigmaLabel(float);
  void slotCompleteSelectedTraj(void);
  void slotUpdateTrajSpinBox(int);
  void slotShowStatus(const QString& status, int timeout);
 signals:
  void keyDownPressed(void);
  void parametersLoaded(void);
  void showInformationRequested(const QString& information);
  void showStatusRequested(const QString& status, int timeout);
  void sigmaUpdated(float);
protected:
  virtual void closeEvent(QCloseEvent *event);
  virtual void keyPressEvent(QKeyEvent* event);
  bool eventFilter(QObject *watched, QEvent *event);

private slots:
  bool slotSetWorkspace(void);
  void slotShowInformation(const QString& information);

  void slotLoadParameters(bool show_message = true);
  void slotUpdateParameterEditor(void);
  void slotSaveParameters(void);

private:
  Q_DISABLE_COPY(MainWindow)

  void loadSettings(void);
  void saveSettings(void);
  void saveStatusLog(void);
  void loadParameters(void);

  Ui::MainWindowClass ui_;

  std::string workspace_;
  SceneWidget* scene_widget_;
  QCustomPlot * heading_plot_;

  bool allow_close_;
};

#endif // MainWindow_H
