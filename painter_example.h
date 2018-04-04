#pragma once

#include "painter_base.h"
#include "painter_destination.h"
#include "painter_source.h"
#include <QOpenGLTexture>
#include <QScopedPointer>

class PainterExample
        : public PainterDestination
{
public:
    PainterExample();
    ~PainterExample();

    bool paint(QPainter * p, const QSize& size, bool wireframe);
private:
    static QOpenGLTexture * createGradient(const QVector<float>& f);
    typedef QScopedPointer<QOpenGLTexture> PTEX;
    PTEX mKittenTexture;
    PTEX mHeartTexture;
    PTEX mWaterTexture;
    GLuint mGradTex = 0;
};
