#include "testwindow.h"
#include <QPainter>
#include <QTimerEvent>
#include <QKeyEvent>
TestWindow::TestWindow(QWidget *parent) : QWidget(parent)
{
    setFixedSize(640,640);
}
void TestWindow::paintEvent(QPaintEvent *event) {
    QPainter p;
    QSize s = size();
    p.begin(this);
    p.fillRect(0, 0, s.width(), s.height(), Qt::black);
    if (!ex.paint(&p, s, wireframe)) {
        if (timer<0) timer = startTimer(33);
    }
    p.end();
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
