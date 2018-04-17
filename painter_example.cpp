#include "painter_example.h"
#include <QPainter>
#include <QDebug>
#include <stdexcept>
#include <QLinearGradient>

static const float line[] = {
    20, 20, 100, 20
};
static const float join1[] = {
    20, 50, 60, 80, 100, 50, 140, 80
};
static const float join2[] = {
    20, 100, 60, 280, 100, 100, 140, 280
};
static const float shape[] = {
    -1, 4, 0, 2, 1, 5,
    2, 2,  5, 4, 4, 1,
    6, 0, 4, -1, 5, -3,
    3, -2, 3, -5, 1, -3,
    -1, -6, -2, -2,
    -5, -3, -4, -1, -6, 0,
    -4, 1, -4, 3, -2, 2
};


static const float LetterH[] = {
    37,34,
    43,59,
    46,100,
    43,167,
    37,170,
    31,163,
    32,154,
    54,121,
    63,110,
    80,100,
    91,97,
    100,83,
    99,57,
    92,45,
    83,55,
    83,88,
    85,133,
    91,175
};
static const float LetterE[] = {
    112,128,
    127,135,
    146,134,
    161,126,
    162,117,
    149,103,
    132,113,
    124,132,
    125,147,
    141,169,
    160,169,
    168,164
};
static const float LetterL[] = {
    172,46,
    180,46,
    188,62,
    190,80,
    192,157,
    195,164,
    210,169,
};
static const float LetterL2[] = {
    218,44,
    226,53,
    230,64,
    228,102,
    228,125,
    234,162,
    242,170,
    255,170
};
static const float LetterO[] = {
    277,124,
    283,156,
    297,163,
    326,143,
    320,125,
    302,113,
    283,117,
};


unsigned char sgrad[40] = {
    0,  0,  0, 60,
    0, 20, 20, 60,
    0, 80, 40, 60,
   40, 80, 40, 70,
   50, 80, 20, 70,

   80, 80, 20, 80,
  180,120, 60, 80,
  250,250,120,180,
  120,250,250,190,
   30, 30, 30, 20,
};
#define NOF(xx) (sizeof(xx)/sizeof(xx[0])/2)
PainterExample::PainterExample()
{
    qDebug()<<"PainterExample constructor started";
    if (!PainterBase::beginPaint())
        throw std::runtime_error("cannot bind context");
    auto $ = PainterBase::gl();

    try {
        mKittenTexture.reset(new QOpenGLTexture(QImage(":/icons/kitten.png")));
        mHeartTexture.reset(new QOpenGLTexture(QImage(":/icons/heart.png")));
        mWaterTexture.reset(new QOpenGLTexture(QImage(":/icons/water.png")));
        mKittenTexture->setWrapMode(QOpenGLTexture::Repeat);
        mHeartTexture->setWrapMode(QOpenGLTexture::Repeat);
        mWaterTexture->setWrapMode(QOpenGLTexture::Repeat);
        if (!mKittenTexture->create())
            throw std::runtime_error("failed to create kitten texture");
        if (!mHeartTexture->create())
            throw std::runtime_error("failed to create heart texture");
        if (!mWaterTexture->create())
            throw std::runtime_error("failed to create water texture");

        $->glGenTextures(1, &mGradTex);
        $->glBindTexture(GL_TEXTURE_1D, mGradTex);
        $->glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, 10, 0, GL_RGBA, GL_UNSIGNED_BYTE, sgrad);
        $->glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
        $->glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        $->glBindTexture(GL_TEXTURE_1D, 0);
    } catch(...) {
        if (mGradTex) $->glDeleteTextures(1, &mGradTex);
        qDebug()<<"failed to create PainterExample";
        throw;
    }

    PainterSource::setPath(shape, NOF(shape), true);

    PainterSource::setTransform(200, 100, 40, 10 );
    PainterSource::setTexture(mKittenTexture->textureId(), 0.0, 0.0, 0.13);
    PainterSource::fillPath();

    PainterSource::setPenSize(0.3, 0.0, true);
    PainterSource::setTexture(0, 0.0, 0.0, 0.0);
    PainterSource::setRGBA(1.0, 0.7, 0.2, 1.0);
    PainterSource::strokePath();
    /*PainterSource::setPath(line, NOF(line), false);
    PainterSource::strokePath();
    PainterSource::setPath(join1, NOF(join1), false);
    PainterSource::strokePath();
    PainterSource::setPath(join2, NOF(join2), false);
    PainterSource::strokePath();*/

    PainterSource::setTransform(0, 400, 0, 1);
    PainterSource::setRGBA(1.0, 0.0, 0.2, 1.0);
    PainterSource::setPenSize(15.0, 0.0, true);
    PainterSource::setPath(LetterH, NOF(LetterH), false);
    PainterSource::strokePath();
    PainterSource::setTexture(mWaterTexture->textureId(), 0.0, 0.0, 0.13);
    PainterSource::setPath(LetterL, NOF(LetterL), false);
    PainterSource::setPath(LetterE, NOF(LetterE), false);
    PainterSource::strokePath();
    PainterSource::setTexture(mHeartTexture->textureId(), 0.0, 0.0, 0.13);
    PainterSource::setPath(LetterL, NOF(LetterL), false);
    PainterSource::strokePath();
    PainterSource::setLinearGradient(mGradTex, 100.0, 50.0, 100.0, 151.0);
    PainterSource::setPath(LetterL2, NOF(LetterL2), false);
    PainterSource::strokePath();
    PainterSource::setRadialGradient(mGradTex, 320.0, 150.0, 50.0);
    PainterSource::setPath(LetterO, NOF(LetterO), true);
    PainterSource::strokePath();
    PainterBase::endPaint();
    qDebug()<<"PainterExample constructed";

}
PainterExample::~PainterExample()
{
    PainterBase::beginPaint();
    auto $ = PainterBase::gl();
    if (mGradTex) $->glDeleteTextures(1, &mGradTex);
    PainterBase::endPaint();
}
bool PainterExample::paint(QPainter * p, const QSize& size, bool wireframe)
{

    /*if (!PainterDestination::beginPaint(size))
        return false;

    PainterSource::paint(size, wireframe);
    PainterDestination::endPaint();

    QImage image = PainterDestination::image();
    if (!image.isNull())
        p->drawImage(0,0,image);
    return !image.isNull();*/

    return paintOnQPainter(p, size, wireframe);
}

