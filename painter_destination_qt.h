#pragma once

#include "painter_source.h"
#include <QImage>
#include <QPainter>

class PainterDestinationQt : public PainterSource {
public:
    PainterDestinationQt();
    virtual ~PainterDestinationQt();

    bool paintOnQPainter(QPainter* p, const QSize& size, bool wireframe);
private:
    GLuint mFramebuffer = 0;
    GLuint mBackBuffer  = 0;
    GLuint mFrontBuffer = 0;
    GLuint mDepthBuffer = 0;
    GLsync mSync        = NULL;
    QSize  mSize;
};
