#include <QResizeEvent>

#include "osg_viewer_widget.h"


ThreadedPainter::ThreadedPainter(OSGViewerWidget *widget)
  : gl_widget_(widget),
    do_rendering_(true) {
}

ThreadedPainter::~ThreadedPainter(void) {

}

void ThreadedPainter::start() {
  gl_widget_->AdapterWidget::makeCurrent();

  do_rendering_ = true;
  startTimer(20);
}

void ThreadedPainter::stop() {
  std::lock_guard<std::mutex> lock(mutex_);
  do_rendering_ = false;
}

void ThreadedPainter::timerEvent(QTimerEvent *event) {
  std::lock_guard<std::mutex> lock(mutex_);
  if (!do_rendering_) {
    killTimer(event->timerId());
    QThread::currentThread()->quit();
    return;
  }

  gl_widget_->updateGL();
}