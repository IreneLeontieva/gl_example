#include "painter_base.h"
#include <stdexcept>
#include <QDebug>


#define N "\n"
static const char * __shader_vertex =
"fuck you nvidia!";


PainterBase::PainterBase()
{
    qDebug()<<"PainterBase constructor started";
    try {
        mOffscreen = new QOffscreenSurface();
        auto fmt = mOffscreen->format();
        fmt.setVersion(3,3);
        fmt.setRenderableType(QSurfaceFormat::OpenGL);
        fmt.setProfile(QSurfaceFormat::CoreProfile);
        mOffscreen->setFormat(fmt);
        mOffscreen->create();
        if (!mOffscreen->isValid())
            throw std::runtime_error("failed to create offscreen surface");
        mContext = new QOpenGLContext();
        mContext->setFormat(mOffscreen->format());
        if (!mContext->create())
            throw std::runtime_error("failed to create opengl context");
        qDebug()<<mContext->format();


        if (!mContext->makeCurrent(mOffscreen))
            throw std::runtime_error("failed to make opengl context current");
        mGL = mContext->versionFunctions<GLFuncs>();
        if (!mGL->initializeOpenGLFunctions())
            throw std::runtime_error("failed to make initialize GL functions");
        mContext->doneCurrent();
    } catch(...) {
        if (mContext)
            mContext->doneCurrent();
        delete mContext;
        delete mOffscreen;
        throw;
    }
    qDebug()<<"PainterBase constructed";
}
bool PainterBase::beginPaint()
{
    if (++mBindCount == 1) {
        return mContext->makeCurrent(mOffscreen);
    }
    return true;
}
void PainterBase::endPaint()
{
    if (--mBindCount == 0) {
        mContext->doneCurrent();
    }
}

PainterBase::~PainterBase()
{
    qDebug()<<"deleting PainterBase";
    if (mContext)
        mContext->doneCurrent();
    delete mContext;
    delete mOffscreen;
}

