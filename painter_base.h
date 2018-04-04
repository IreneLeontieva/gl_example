#pragma once
#include <QOffscreenSurface>
#include <QOpenGLContext>
#include <QOpenGLFunctions_3_0>
#include <QOpenGLFunctions_3_1>
#include <QOpenGLFunctions_3_2_Core>
#include <QOpenGLFunctions_3_3_Core>

typedef QOpenGLFunctions_3_3_Core GLFuncs;

struct PainterBase {
public:
    PainterBase();
    virtual ~PainterBase();

    //opengl context and some objects
    GLFuncs * gl() const { return mGL; }
    bool beginPaint();
    void endPaint();
private:
    QOffscreenSurface         * mOffscreen = nullptr;
    QOpenGLContext            * mContext   = nullptr;
    GLFuncs                   * mGL = nullptr;
    int                         mBindCount = 0;
    Q_DISABLE_COPY(PainterBase)
};

