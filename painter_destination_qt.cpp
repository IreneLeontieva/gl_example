#include "painter_destination_qt.h"
#include <stdexcept>
#include <QDebug>
#include <GL/glu.h>

PainterDestinationQt::PainterDestinationQt()
{
    qDebug()<<"PainterDestinationQt constructor started";
    if (!PainterBase::beginPaint())
        throw std::runtime_error("cannot bind context");
    auto $ = PainterBase::gl();

    try {
        GLenum e;

        while($->glGetError());
        $->glGenFramebuffers(1, &mFramebuffer);
        $->glGenTextures(1, &mFrontBuffer);
        $->glGenTextures(1, &mBackBuffer);
        $->glGenTextures(1, &mDepthBuffer);

        e = $->glGetError();
        if (e) {
            qDebug()<<">"<<((const char*)gluErrorString(e));
            throw std::runtime_error("pending GL errors(fb objects)");
        }
    } catch(...) {
        $->glFinish();
        if (mSync) $->glDeleteSync(mSync);
        if (mFramebuffer) $->glDeleteFramebuffers(1, &mFramebuffer);
        if (mFrontBuffer) $->glDeleteTextures(1, &mFrontBuffer);
        if (mBackBuffer)  $->glDeleteTextures(1, &mBackBuffer);
        if (mDepthBuffer) $->glDeleteTextures(1, &mDepthBuffer);
        throw;
    }
    PainterBase::endPaint();
    qDebug()<<"PainterDestinationQt constructed";
}
PainterDestinationQt::~PainterDestinationQt()
{
    PainterBase::beginPaint();
    auto $ = PainterBase::gl();
    $->glFinish();
    if (mSync) $->glDeleteSync(mSync);
    if (mFramebuffer) $->glDeleteFramebuffers(1, &mFramebuffer);
    if (mFrontBuffer) $->glDeleteTextures(1, &mFrontBuffer);
    if (mBackBuffer)  $->glDeleteTextures(1, &mBackBuffer);
    if (mDepthBuffer) $->glDeleteTextures(1, &mDepthBuffer);
    PainterBase::endPaint();
}

bool PainterDestinationQt::paintOnQPainter(QPainter * p, const QSize &size, bool wireframe)
{
    if (size.isEmpty())
        return false;

    auto $ = PainterBase::gl();

    bool ret = false;
    //get current framebuffer
    GLint dst_buf = 0;
    $->glGetIntegerv(GL_DRAW_FRAMEBUFFER_BINDING, &dst_buf);

    if (p) p->beginNativePainting();
    do {
        auto originalContext = QOpenGLContext::currentContext();
        QSurface * originalSurface = originalContext->surface();
        //first, we bind "our context" because framebuffers cannot be shared
        bool nothing_to_draw_yet = false;
        if (!PainterBase::beginPaint())
            break;
        do {
            if (mSync) {
                auto e = $->glClientWaitSync(mSync, GL_SYNC_FLUSH_COMMANDS_BIT, 10000000);
                if (e == GL_TIMEOUT_EXPIRED) {
                    qDebug("wont paint1");
                    break;
                }
                if (e == GL_WAIT_FAILED)
                    $->glFinish();
                $->glDeleteSync(mSync);
                mSync = nullptr;
            }

            //in case when we have to resize a buffer, we won't a have a suitable texture to draw
            $->glBindFramebuffer(GL_FRAMEBUFFER, mFramebuffer);
            if (size != mSize) {
                while($->glGetError());
                $->glFramebufferTexture2D(
                            GL_DRAW_FRAMEBUFFER,
                            GL_COLOR_ATTACHMENT0,
                            GL_TEXTURE_2D,
                            0, 0);
                $->glFramebufferTexture2D(
                            GL_DRAW_FRAMEBUFFER,
                            GL_DEPTH_STENCIL_ATTACHMENT,
                            GL_TEXTURE_2D,
                            0, 0);

                $->glBindTexture(
                            GL_TEXTURE_2D,
                            mFrontBuffer);
                $->glTexImage2D(
                            GL_TEXTURE_2D,
                            0,
                            GL_RGBA,
                            size.width(),
                            size.height(),
                            0,
                            GL_RGBA,
                            GL_UNSIGNED_INT_8_8_8_8_REV,
                            NULL);
                $->glTexParameteri(
                            GL_TEXTURE_2D,
                            GL_TEXTURE_MAG_FILTER,
                            GL_LINEAR);
                $->glTexParameteri(
                            GL_TEXTURE_2D,
                            GL_TEXTURE_MIN_FILTER,
                            GL_LINEAR);

                $->glBindTexture(
                            GL_TEXTURE_2D,
                            mBackBuffer);
                $->glTexImage2D(
                            GL_TEXTURE_2D,
                            0,
                            GL_RGBA,
                            size.width(),
                            size.height(),
                            0,
                            GL_RGBA,
                            GL_UNSIGNED_INT_8_8_8_8_REV,
                            NULL);
                $->glTexParameteri(
                            GL_TEXTURE_2D,
                            GL_TEXTURE_MAG_FILTER,
                            GL_LINEAR);
                $->glTexParameteri(
                            GL_TEXTURE_2D,
                            GL_TEXTURE_MIN_FILTER,
                            GL_LINEAR);

                $->glBindTexture(
                            GL_TEXTURE_2D,
                            mDepthBuffer);
                $->glTexImage2D(
                            GL_TEXTURE_2D,
                            0,
                            GL_DEPTH24_STENCIL8,
                            size.width(),
                            size.height(),
                            0,
                            GL_DEPTH_STENCIL,
                            GL_UNSIGNED_INT_24_8,
                            NULL);
                $->glTexParameteri(
                            GL_TEXTURE_2D,
                            GL_TEXTURE_MAG_FILTER,
                            GL_LINEAR);
                $->glTexParameteri(
                            GL_TEXTURE_2D,
                            GL_TEXTURE_MIN_FILTER,
                            GL_LINEAR);
                $->glFramebufferTexture2D(
                            GL_DRAW_FRAMEBUFFER,
                            GL_DEPTH_STENCIL_ATTACHMENT,
                            GL_TEXTURE_2D,
                            mDepthBuffer, 0);
                $->glBindTexture(GL_TEXTURE_2D, 0);
                GLenum e = $->glGetError();
                if (e) {
                    qDebug()<<">"<<((const char*)gluErrorString(e));
                    throw std::runtime_error("pending GL errors(pixel buffer)");
                }
                mSize = size;
                nothing_to_draw_yet = true;
            } else {
                std::swap(mFrontBuffer, mBackBuffer);
            }
            /*
             * clear buffer
             */
            GLenum e = $->glGetError();
            if (e) {
                qDebug()<<">a"<<((const char*)gluErrorString(e));
            }
            $->glFramebufferTexture2D(
                        GL_FRAMEBUFFER,
                        GL_COLOR_ATTACHMENT0,
                        GL_TEXTURE_2D,
                        mBackBuffer, 0);
            $->glViewport(0, 0, size.width(), size.height());
            $->glClearColor(1, 1, 1, 1);
            $->glClearDepth(0);
            $->glClearStencil(0);
            $->glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
            PainterSource::paint(size, wireframe);
            $->glFlush();
        } while(0);
        PainterBase::endPaint();
        //second, we bind original context
        originalContext->makeCurrent(originalSurface);
        $->glBindFramebuffer(GL_FRAMEBUFFER,dst_buf);
        if (!nothing_to_draw_yet) {
            GLuint fbo = 0;
            $->glGenFramebuffers(1, &fbo);
            $->glBindFramebuffer(GL_READ_FRAMEBUFFER,fbo);
            $->glFramebufferTexture2D(
                        GL_READ_FRAMEBUFFER,
                        GL_COLOR_ATTACHMENT0,
                        GL_TEXTURE_2D,
                        mFrontBuffer, 0);
            $->glBlitFramebuffer(0, 0, size.width(), size.height(),
                                 0, 0, size.width(), size.height(),
                                 GL_COLOR_BUFFER_BIT, GL_NEAREST);
            $->glBindFramebuffer(GL_READ_FRAMEBUFFER,0);
            $->glDeleteFramebuffers(1, &fbo);
            ret = true;
        } else {
            $->glClearColor(1, 1, 1, 1);
            $->glClearDepth(0);
            $->glClearStencil(0);
            $->glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);
        }
    } while(0);
    if (p) p->endNativePainting();
    return ret;
}

