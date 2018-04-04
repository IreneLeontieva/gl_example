#include "painter_destination.h"
#include <stdexcept>
#include <QDebug>
#include <GL/glu.h>

PainterDestination::PainterDestination()
{
    qDebug()<<"PainterDestination constructor started";
    if (!PainterBase::beginPaint())
        throw std::runtime_error("cannot bind context");
    auto $ = PainterBase::gl();

    try {
        GLenum e;

        while($->glGetError());
        $->glGenFramebuffers(1, &mFramebuffer);
        $->glGenTextures(1, &mColorBuffer);
        $->glGenTextures(1, &mDepthBuffer);
        $->glGenBuffers(1, &mFront);
        $->glGenBuffers(1, &mBack);

        e = $->glGetError();
        if (e) {
            qDebug()<<">"<<((const char*)gluErrorString(e));
            throw std::runtime_error("pending GL errors(fb objects)");
        }
    } catch(...) {
        $->glFinish();
        if (mSync) $->glDeleteSync(mSync);
        if (mFramebuffer) $->glDeleteFramebuffers(1, &mFramebuffer);
        if (mColorBuffer) $->glDeleteTextures(1, &mColorBuffer);
        if (mDepthBuffer) $->glDeleteTextures(1, &mDepthBuffer);
        if (mFront) $->glDeleteBuffers(1, &mFront);
        if (mBack) $->glDeleteBuffers(1, &mBack);
        throw;
    }
    PainterBase::endPaint();
    qDebug()<<"PainterDestination constructed";
}
PainterDestination::~PainterDestination()
{
    PainterBase::beginPaint();
    auto $ = PainterBase::gl();
    $->glFinish();
    if (mSync) $->glDeleteSync(mSync);
    if (mFramebuffer) $->glDeleteFramebuffers(1, &mFramebuffer);
    if (mColorBuffer) $->glDeleteTextures(1, &mColorBuffer);
    if (mDepthBuffer) $->glDeleteTextures(1, &mDepthBuffer);
    if (mPinnedMemory) {
        $->glBindBuffer(GL_PIXEL_PACK_BUFFER, mFront);
        $->glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
        $->glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    }
    if (mFront) $->glDeleteBuffers(1, &mFront);
    if (mBack) $->glDeleteBuffers(1, &mBack);
    PainterBase::endPaint();
}

bool PainterDestination::beginPaint(const QSize &size)
{
    if (!PainterBase::beginPaint())
        return false;
    if (size.isEmpty())
        return false;

    auto $ = PainterBase::gl();
    if (mSync) {
        auto e = $->glClientWaitSync(mSync, GL_SYNC_FLUSH_COMMANDS_BIT, 10000000);
        if (e == GL_TIMEOUT_EXPIRED) {
            PainterBase::endPaint();
            return false;
        }

        if (e == GL_WAIT_FAILED)
            $->glFinish();
        $->glDeleteSync(mSync);
        mSync = nullptr;
    }

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
                    mColorBuffer);
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
        $->glFramebufferTexture2D(
                    GL_DRAW_FRAMEBUFFER,
                    GL_COLOR_ATTACHMENT0,
                    GL_TEXTURE_2D,
                    mColorBuffer, 0);
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
            throw std::runtime_error("pending GL errors(framebuffer)");
        }
        $->glBindBuffer(GL_PIXEL_PACK_BUFFER, mFront);
        $->glBufferData(GL_PIXEL_PACK_BUFFER,
                        size.width()*size.height()*4+256,
                        NULL,
                        GL_STREAM_READ);
        $->glBindBuffer(GL_PIXEL_PACK_BUFFER, mBack);
        $->glBufferData(GL_PIXEL_PACK_BUFFER,
                        size.width()*size.height()*4+256,
                        NULL,
                        GL_STREAM_READ);
        $->glBindBuffer(GL_PIXEL_PACK_BUFFER,
                        0);
        e = $->glGetError();
        if (e) {
            qDebug()<<">"<<((const char*)gluErrorString(e));
            throw std::runtime_error("pending GL errors(pixel buffer)");
        }
        e = $->glCheckFramebufferStatus(GL_FRAMEBUFFER);

        mSize = size;
    } else {
        if (mPinnedMemory) {
            $->glBindBuffer(GL_PIXEL_PACK_BUFFER, mFront);
            $->glUnmapBuffer(GL_PIXEL_PACK_BUFFER);
            $->glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
            mPinnedMemory = nullptr;
        }

        std::swap(mFront, mBack);
        $->glBindBuffer(GL_PIXEL_PACK_BUFFER, mFront);
        mPinnedMemory = $->glMapBuffer(GL_PIXEL_PACK_BUFFER, GL_READ_ONLY);
        $->glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    }
    /*
     * clear buffer
     */
    $->glViewport(0, 0, size.width(), size.height());
    $->glClearColor(1, 1, 1, 1);
    $->glClearDepth(0);
    $->glClearStencil(0);
    $->glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT|GL_STENCIL_BUFFER_BIT);

    return true;
}
void PainterDestination::endPaint()
{
    /*
     * done with painting, schedule content copying to PBO
     */
    auto $ = PainterBase::gl();
    $->glBindFramebuffer(GL_FRAMEBUFFER, 0);
    $->glFlush();
    $->glBindBuffer(GL_PIXEL_PACK_BUFFER, mBack);
    $->glBindFramebuffer(GL_READ_FRAMEBUFFER, mFramebuffer);
    $->glReadBuffer(GL_COLOR_ATTACHMENT0);
    $->glReadPixels(0, 0, mSize.width(), mSize.height(), GL_BGRA, GL_UNSIGNED_INT_8_8_8_8_REV, 0);
    $->glBindFramebuffer(GL_READ_FRAMEBUFFER, 0);
    $->glBindBuffer(GL_PIXEL_PACK_BUFFER, 0);
    $->glFlush();
    mSync = $->glFenceSync(GL_SYNC_GPU_COMMANDS_COMPLETE, 0);
    PainterBase::endPaint();
}
QImage PainterDestination::image() const {
    if (!mPinnedMemory)
        return QImage();

    return QImage((const uchar*)mPinnedMemory,
                  mSize.width(),
                  mSize.height(),
                  4*mSize.width(),
                  QImage::Format_ARGB32);
}
