#pragma once

#include "painter_source.h"
#include <QImage>

class PainterDestination : public PainterSource {
public:
    PainterDestination();
    virtual ~PainterDestination();

    bool paintOnQPainter(QPainter* p, const QSize& size, bool wireframe);
    bool beginPaint(const QSize& size);
    void endPaint();
    QImage image() const;
private:
    GLuint mFramebuffer = 0;
    GLuint mColorBuffer = 0;
    GLuint mDepthBuffer = 0;
    GLuint mFront       = 0;
    GLuint mBack        = 0;
    GLsync mSync        = NULL;
    void  *mPinnedMemory= nullptr;
    QSize  mSize;

};
