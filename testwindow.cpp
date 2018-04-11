#include "testwindow.h"
#include <QPainter>
#include <QTimerEvent>
#include <QKeyEvent>
TestWindow::TestWindow(QWidget *parent)
: QOpenGLWidget(parent)
//: QWidget(parent)
{
    QSurfaceFormat fmt = format();
    fmt.setVersion(3,3);
    fmt.setProfile(QSurfaceFormat::CoreProfile);
    setFormat(fmt);
    setFixedSize(640,640);
//    if (!ex) ex = new PainterExample();
}
TestWindow::~TestWindow()
{
    delete ex;
}

/*void TestWindow::paintEvent(QPaintEvent *event) {
    QPainter p;
    QSize s = size();
    p.begin(this);
    p.fillRect(0, 0, s.width(), s.height(), Qt::black);
    if (!ex->paintOnQPainter(&p, s, wireframe)) {
        if (timer<0) timer = startTimer(33);
    }
    p.end();
}*/

void TestWindow::initializeGL()
{
    qDebug()<<"CI="<<context();
    if (!ex) ex = new PainterExample();
//    context()->makeCurrent(this);
}

void TestWindow::paintGL() {
    if (!ex) return;
    QSize s = size();
    if (!ex->paintOnQPainter(nullptr, s, wireframe)) {
        if (timer<0) timer = startTimer(33);
    }
//    context()->makeCurrent(this);
}
void TestWindow::timerEvent(QTimerEvent *event) {
    if (event->timerId() != timer) {
        QWidget::timerEvent(event);
    } else {
        killTimer(timer);
        timer = -1;
        update();
    }
}
void TestWindow::keyPressEvent(QKeyEvent *event)
{
    if (event->key() == Qt::Key_Space) {
        wireframe = !wireframe;
        update();
    }
}
