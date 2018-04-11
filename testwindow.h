#ifndef TESTWINDOW_H
#define TESTWINDOW_H

#include <QOpenGLWidget>
#include <QWidget>
#include "painter_example.h"

class TestWindow
        : public QOpenGLWidget
        //: public QWidget
{
    Q_OBJECT
public:
    explicit TestWindow(QWidget *parent = 0);
    ~TestWindow();

private:
    PainterExample * ex = nullptr;
    int timer           = -1;
    bool wireframe      = false;

    void initializeGL();
    void paintGL();
    void timerEvent(QTimerEvent *event);
    void keyPressEvent(QKeyEvent *event);
   // void paintEvent(QPaintEvent *event);
};

#endif // TESTWINDOW_H
