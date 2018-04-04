#pragma once

#include "painter_base.h"
#include <QByteArray>
#include <QOpenGLShaderProgram>
#include <GL/glu.h>

struct VertexData {
    float x,   y,   z,   _0;
    float e1x, e1y, e1z, e1w;
    float e2x, e2y, e2z, e2w;
    float e3x, e3y, e3z, e3w;

    float fill_xx, fill_xy, fill_x0, fill_G1;
    float fill_yx, fill_yy, fill_y0, fill_G2;
    float fill_r,  fill_g,  fill_b,  fill_a;

    float mat_xx, mat_xy, mat_x0, _u0;
    float mat_yx, mat_yy, mat_y0, _y0;
};

class ShaderSource : public PainterBase
{
public:
    enum : unsigned { PRIM_RESTART = 0x7fffffffu };
    ShaderSource ();
    ~ShaderSource();
    //paints on provided framebuffer
    bool paint(const QSize &size);
    //clears caches primities
    void clear();

    //filling the buffers
    //source color settings
    void setRGB(float r, float g, float b);
    void setRGBA(float r, float g, float b, float a);
    void setA(float a);
    void setTexture(GLuint texture);
    void setGradient(GLuint gtexture,
                     float start,
                     float stop,
                     bool radial);
    void setPenSize(float pen, float miter, bool rounded);
    void setPath(double * coords, int ncoords, bool closed);
    void stokePath();
    void fillPath();
private:
    GLUtesselator        * mTesselator = nullptr;
    QOpenGLShaderProgram * mProgram    = nullptr;
    GLint                  mPaintDPS  = -1;
    GLuint                 mVAO = 0;
    GLuint                 mVBO = 0;
    GLuint                 mIBO = 0;
    enum { VBO_SIZE = 65536, IBO_SIZE = 65536 };
    QByteArray             mVertices;
    QByteArray             mIndices;
    //tesselated path
    VertexData             mVertexSample;
    QByteArray             mTessBuffer;
    int                    mPathSize = 0;
    int                    mPathExt  = 0;
    bool                   mPathClosed = false;
    //aux data for tesselator
    GLenum                 mTessType;
    int                    mTessVertexBase;
    int                    mTessVertexCount;
    //save these position for case of error
    int                    mTessFirstVertex;
    int                    mTessFirstIndex;
    bool                   mTessBroken;

    void insertRCS(int prev);

    //cached values
    GLuint                 mTex2 = 0;
    GLuint                 mTex1 = 0;
    float                  mPenW = 1.0;
    float                  mMiter= 1.0;
    bool                   mRounded = false;
    bool                   mGradRadial= false;
    struct RenderCommand {
        int    mFirstIndex;
        int    mNumIndices;
        GLuint mTex1, mTex2;
    };
    QVector<RenderCommand> mRCS;



    static void tessBegin(GLenum type, void *self);
    static void tessVertex(void * vertex_data, void * self);
    static void tessEnd(void * self);
    static void tessCombineData(GLdouble coords[3], void *vertex_data[4],
                                GLfloat weight[4], void **outData,
                                void * self);
    static void tessError(GLenum errno, void *self);


};


