#ifndef TESTWINDOW_H
#define TESTWINDOW_H

#include <QWidget>
#include "painter_example.h"
class TestWindow : public QWidget
{
    Q_OBJECT
public:
    explicit TestWindow(QWidget *parent = 0);

    void paintEvent(QPaintEvent *event);
    void timerEvent(QTimerEvent *event);
private:
    PainterExample ex;
    int timer = -1;
    bool wireframe = false;
    void keyPressEvent(QKeyEvent *event);
};

#endif // TESTWINDOW_H
