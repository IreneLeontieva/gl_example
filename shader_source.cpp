#include "shader_source.h"
#include <stdexcept>
#include <cmath>
#include <assert.h>
#include <QDebug>
static const char * __vertex;
static const char * __fragment;

ShaderSource::ShaderSource()
    : PainterBase()
{
    auto $ = PainterBase::gl();
    try {
        while($->glGetError());

        mProgram = new QOpenGLShaderProgram();
        bool ok;
        ok = mProgram->addShaderFromSourceCode(QOpenGLShader::Vertex,
                                               QByteArray(__vertex));
        if (!ok) {
            qDebug()<<mProgram->log();
            throw std::runtime_error("failed to compile vertex shader");
        }
        ok = mProgram->addShaderFromSourceCode(QOpenGLShader::Fragment,
                                               QByteArray(__fragment));
        if (!ok) {
            qDebug()<<mProgram->log();
            throw std::runtime_error("failed to compile fragment shader");
        }
        ok = mProgram->link();
        if (!ok) {
            qDebug()<<mProgram->log();
            throw std::runtime_error("failed to link shader");
        }
        mProgram->bind();
        mProgram->setUniformValue("source1", 0);
        mProgram->setUniformValue("source2", 1);
        mPaintDPS  = mProgram->uniformLocation("display");
        mProgram->release();
        if (mPaintDPS < 0)
            throw std::runtime_error("uniforms not found");

        GLenum e = $->glGetError();
        if (e) {
            qDebug()<<">"<<gluErrorString(e);
            throw std::runtime_error("pending GL errors");
        }

        $->glGenVertexArrays(1, &mVAO);
        $->glGenBuffers(1, &mVBO);
        $->glGenBuffers(1, &mIBO);
        e = $->glGetError();
        if (e) {
            qDebug()<<">"<<gluErrorString(e);
            throw std::runtime_error("pending GL errors");
        }

        const char * attrs[] = {
            "vertex",
            "edge1",
            "edge2",
            "edge3",
            "fill1",
            "fill2",
            "fill3",
            "trX",
            "trY",
            NULL
        };
        $->glBindVertexArray(mVAO);
        $->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mIBO);
        $->glBufferData(GL_ELEMENT_ARRAY_BUFFER, IBO_SIZE, NULL, GL_STREAM_DRAW);
        $->glBindBuffer(GL_ARRAY_BUFFER, mVBO);
        $->glBufferData(GL_ARRAY_BUFFER, VBO_SIZE, NULL, GL_STREAM_DRAW);
        try {
            for(int j = 0; attrs[j]; ++j) {
                auto a0 = mProgram->attributeLocation(attrs[j]);
                if (a0<0) {
                    qDebug()<<"attribute "<<attrs[j]<<" not found";
                    throw std::runtime_error("program attribute not found");
                }
                $->glVertexAttribPointer(a0, 4, GL_FLOAT, GL_FALSE,
                                         16*9, (GLvoid*)(long)(16*j));
                $->glEnableVertexAttribArray(a0);
            }
        }catch(...) {
            $->glBindVertexArray(0);
            $->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
            $->glBindBuffer(GL_ARRAY_BUFFER, 0);
            throw;
        }
        $->glBindVertexArray(0);
        $->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
        $->glBindBuffer(GL_ARRAY_BUFFER, 0);

        mVertices.reserve(VBO_SIZE);
        mIndices.reserve(IBO_SIZE);
        mTessBuffer.reserve(VBO_SIZE);

        mTesselator = gluNewTess();
        if (!mTesselator)
            throw std::bad_alloc();
        gluTessCallback(mTesselator, GLU_TESS_BEGIN_DATA, (_GLUfuncptr)tessBegin);
        gluTessCallback(mTesselator, GLU_TESS_VERTEX_DATA, (_GLUfuncptr)tessVertex);
        gluTessCallback(mTesselator, GLU_TESS_COMBINE_DATA, (_GLUfuncptr)tessCombineData);
        gluTessCallback(mTesselator, GLU_TESS_END_DATA, (_GLUfuncptr)tessEnd);
        gluTessCallback(mTesselator, GLU_TESS_ERROR_DATA, (_GLUfuncptr)tessError);

        memset(&mVertexSample, 0, sizeof(mVertexSample));
        //clear color
        mVertexSample.fill_a = 1;
        //clear transform
        mVertexSample.mat_xx = 1;
        mVertexSample.mat_yy = 1;

    } catch(...) {
        if (mVAO) $->glDeleteVertexArrays(1, &mVAO);
        if (mVBO) $->glDeleteBuffers(1, &mVBO);
        if (mIBO) $->glDeleteBuffers(1, &mIBO);
        if (mTesselator)
            gluDeleteTess(mTesselator);
        delete mProgram;
        throw;
    }
}
ShaderSource::~ShaderSource()
{
    auto $ = PainterBase::gl();
    if (mVAO) $->glDeleteVertexArrays(1, &mVAO);
    if (mVBO) $->glDeleteBuffers(1, &mVBO);
    if (mIBO) $->glDeleteBuffers(1, &mIBO);
    if (mTesselator)
        gluDeleteTess(mTesselator);
    delete mProgram;
}
bool ShaderSource::paint(const QSize& size) {
    if (size.isEmpty())
        return false;

    if (mIndices.size() && mRCS.size())
    {
        auto $ = PainterBase::gl();
        if (!mProgram->bind())
            return false;
        mProgram->setUniformValue(mPaintDPS,
                                  1.0f/float(size.width()),
                                  1.0f/float(size.height()));
        $->glBindVertexArray(mVAO);
        for(int i = 0; i < mRCS.size(); ++i) {
            const RenderCommand& rc = mRCS.at(i);
            $->glActiveTexture(GL_TEXTURE0);
            $->glBindTexture(GL_TEXTURE_1D, rc.mTex1);
            $->glActiveTexture(GL_TEXTURE1);
            $->glBindTexture(GL_TEXTURE_2D, rc.mTex2);
            $->glDrawElementsBaseVertex(GL_TRIANGLE_STRIP,
                                        rc.mNumIndices,
                                        GL_UNSIGNED_INT,
                                        NULL,
                                        rc.mFirstIndex);
        }
        $->glBindVertexArray(0);
        mProgram->release();
    }
    return true;
}
void ShaderSource::clear() {
    mIndices.clear();
    mVertices.clear();
    mRCS.clear();
}
void ShaderSource::setPath(double *coords, int ncoords, bool closed) {
    mTessBuffer.clear();
    mPathSize   = ncoords;
    mPathExt    = ncoords;
    mPathClosed = closed;

    mTessBuffer.resize(3*ncoords*2);
    mTessBuffer.detach();
    double * dst = (double*)mTessBuffer.data();
    for(int i = 0; i < ncoords; ++i) {
        dst[3*i]   = coords[2*i];
        dst[3*i+1] = coords[2*i+1];
        dst[3*i]   = 0.0;
    }
}


static inline double TURN(double dx1, double dy1, double dx2, double dy2, double miter) {
    double d = dx2*dy2-dx1*dy1;
    if (d == 0.0) d = 1.0;
    return std::copysign(2.0+miter, d);
}
static inline void EXTENT(double &ox, double &oy,
                          double dx1, double dy1,
                          double dx2, double dy2,
                          double width)
{
    /*
     * let A = (dx1,dy1);
     * let B = (dx2,dy2);
     *
     * let A'= (-dy1,dx1);
     * let B'= (-dy2,dx2);
     *
     * assuming |A|=1, |B|=1,
     * |A'|=1;
     * |B'|=1;
     *
     * we need to find T such:
     * TA' = width;
     * TB' = width;
     *
     * we assume:
     * T = k(A'+B')
     *
     * kA'(A'+B')=width;
     *   -> k(1+A'B')=width;
     * kB'(A'+B')=width;
     *   -> k(1+A'B')=width;
     *
     * k = width/(1+A'B')
     *
     * worst case scenario is U-turn in path
     */
    double c = dx1*dx2+dy1*dy2;
    //FIXME: handle this case better
    c = std::max(c, -0.999999);
    double k = width/(1.0+c);
    double avx = dx1 + dx2;
    double avy = dy1 + dy2;
    ox = -avy*k;
    oy =  avx*k;
}
#define EMIT_VERTEX() \
    mVertices.append((const char*)&mVertexSample, sizeof(mVertexSample));

void ShaderSource::stokePath() {
    if (mPathSize < 2) return;
    assert(mPathExt >= mPathSize);
    assert(uint(mTessBuffer.size()) >= mPathExt*3*sizeof(double));

    //save info for case of error
    mTessBroken     = false;
    mTessFirstIndex = mIndices.size();
    mTessFirstVertex= mVertices.size();
    bool  reallyClosed = mPathClosed & (mPathSize > 2);

    const double * src = (const double*)mTessBuffer.data();
    const double * A, * B, * C, * D, * M;

    A = nullptr; B = src; C = src+3; D = nullptr;
    M = src + mPathSize*3;

    if (reallyClosed)
        A = M-3;
    if (mPathSize > 2)
        D = src+6;

    //create initial 2 vertices as for disjoint line
    double dax,day,dal;
    double ddx,ddy,ddl;
    double dx = C[0] - B[0];
    double dy = C[1] - B[1];
    double dl  = sqrt(dx*dx+dy*dy);
    if (dl < 0.00001) {
        dx = 1.0;
        dy = 0.0;
    } else {
        //normalize vector
        dl = 1.0/dl; dx *= dl; dy *= dl;
    }
    dax = ddx = dx;
    day = ddy = dy;
    mVertexSample.e1x = -dx;
    mVertexSample.e1y = -dy;
    mVertexSample.e1z = B[0]*dx+B[1]*dy;
    mVertexSample.e1w = 0.0;
    mVertexSample.e2x = -dy;
    mVertexSample.e2y = dx;
    mVertexSample.e2z = B[0]*dy-B[1]*dx;
    mVertexSample.e2w =-1.0;
    if (!mRounded) mVertexSample.e2w = mMiter;
    mVertexSample.e3x = dx;
    mVertexSample.e3y = dy;
    mVertexSample.e3z = -C[0]*dx-C[1]*dy;
    mVertexSample.e3w = 0.0;
    if (A) {
        dax = B[0] - A[0];
        day = B[1] - A[1];
        dal  = sqrt(dax*dax+day*day);
        if (dal < 0.00001) {
            dax = 1.0;
            day = 0.0;
        } else {
            //normalize vector
            dal = 1.0/dal; dax *= dal; day *= dal;
        }
        mVertexSample.e1x = -day;
        mVertexSample.e1y =  dax;
        mVertexSample.e1z = B[0]*day-B[1]*dax;
        mVertexSample.e1w = TURN(dax,day,dx,dy,mMiter);
    }

    bool first = true;
    do {
        if (D) {
            ddx = D[0] - C[0];
            ddy = D[1] - C[1];
            ddl  = sqrt(ddx*ddx+ddy*ddy);
            if (ddl < 0.00001) {
                ddx = 1.0;
                ddy = 0.0;
            } else {
                //normalize vector
                ddl = 1.0/ddl; ddx *= ddl; ddy *= ddl;
            }
            mVertexSample.e3x = -ddy;
            mVertexSample.e3y =  ddx;
            mVertexSample.e3z = C[0]*ddy-C[1]*ddx;
            mVertexSample.e3w = TURN(dx,dy,ddx,ddy,mMiter);
        }
        //emit vertex
        double offset_x, offset_y, k;
        if (first) {
            EXTENT(offset_x, offset_y, dax, day, dx, dy, mPenW);
            k = reallyClosed ? mPenW : 0.0;
            mVertexSample.x = B[0] + offset_x - dx*k;
            mVertexSample.x = B[1] + offset_y - dy*k;
            EMIT_VERTEX();
            mVertexSample.x = B[0] - offset_x - dx*k;
            mVertexSample.x = B[1] - offset_y - dy*k;
            EMIT_VERTEX();
            first = false;
        }
        EXTENT(offset_x, offset_y, dx, dy, ddx, ddy, mPenW);
        k = D ? 0.0 : mPenW;
        mVertexSample.x = C[0] + offset_x + dx*k;
        mVertexSample.x = C[1] + offset_y + dy*k;
        EMIT_VERTEX();
        mVertexSample.x = C[0] - offset_x + dx*k;
        mVertexSample.x = C[1] - offset_y + dy*k;
        EMIT_VERTEX();

        //if we are rasterizing a line, just stop here
        if (!D || C == src) break;

        //iterate to next vertex
        A = B;
        B = C;
        C = D;
        //check if we
        //copy parameters, dont recalculate them
        mVertexSample.e1x = mVertexSample.e2x;
        mVertexSample.e1y = mVertexSample.e2y;
        mVertexSample.e1z = mVertexSample.e2z;
        mVertexSample.e1w = mVertexSample.e3w;//that's not typo!
        mVertexSample.e2x = mVertexSample.e3x;
        mVertexSample.e2y = mVertexSample.e3y;
        mVertexSample.e2z = mVertexSample.e3z;
        //check if there is a next vertex
        D  += 3;
        if (D < M) continue;

        if (!reallyClosed) {
            mVertexSample.e3x = ddx;
            mVertexSample.e3y = ddy;
            mVertexSample.e3z = -C[0]*ddx-C[1]*ddy;
            mVertexSample.e3w = 0.0;
            D = nullptr;
            continue;
        }
        D = src;
    } while(true);
    //OMG, I did it! I did it! I did it! I did it! I did it!
    //Damn line rasterizer!
    if (mTessBroken) {
        mIndices.resize(mTessFirstIndex);
        mVertices.resize(mTessFirstVertex);
    } else {
        insertRCS(mTessFirstIndex);
    }
}
void ShaderSource::fillPath() {
    if (!mPathClosed) return;
    if (mPathSize < 3) return;
    assert(mPathExt >= mPathSize);
    assert(uint(mTessBuffer.size()) >= mPathExt*3*sizeof(double));
    //save info for case of error
    mTessBroken     = false;
    mTessFirstIndex = mIndices.size();
    mTessFirstVertex= mVertices.size();
    //disable stroking
    mVertexSample.e1x = mVertexSample.e1y = mVertexSample.e1z = mVertexSample.e1w = 0.0f;
    mVertexSample.e2x = mVertexSample.e2y = mVertexSample.e2z = mVertexSample.e2w = 0.0f;
    mVertexSample.e3x = mVertexSample.e3y = mVertexSample.e3z = mVertexSample.e3w = 0.0f;

    double * dst = (double*)mTessBuffer.data();
    gluTessBeginPolygon(mTesselator, this);
    gluTessBeginContour(mTesselator);
    for(int i = 0; i < mPathSize; ++i) {
        gluTessVertex(mTesselator,
                      dst+i*3,
                      this);
    }
    gluTessEndContour(mTesselator);
    gluTessEndPolygon(mTesselator);

    if (mTessBroken) {
        mIndices.resize(mTessFirstIndex);
        mVertices.resize(mTessFirstVertex);
    } else {
        insertRCS(mTessFirstIndex);
    }
}
//====utility functions=============================
void ShaderSource::insertRCS(int prev) {
    do {//check if RCS can be combined
        if (mRCS.isEmpty())
            break;
        auto last = mRCS.last();
        if (last.mTex1 != mTex1 ||
            last.mTex2 != mTex2)
            break;
        mRCS.last().mNumIndices += (mIndices.size() - prev)/4;
        return;
    } while(0);
    RenderCommand rc;
    rc.mFirstIndex = prev/4;
    rc.mNumIndices += (mIndices.size() - prev)/4;
    rc.mTex1 = mTex1;
    rc.mTex2 = mTex2;
}
//====tesselator=================================
void ShaderSource::tessBegin(GLenum type, void *self) {
    ShaderSource * p = (ShaderSource*)self;
    if (p->mTessBroken) return;
    p->mTessType        = type;
    p->mTessVertexBase  = p->mVertices.size()/
                          sizeof(VertexData);
    p->mTessVertexCount = 0;
}
void ShaderSource::tessVertex(void *vertex_data, void *self) {
    ShaderSource * p = (ShaderSource*)self;
    if (p->mTessBroken) return;

    double * v = (double*)vertex_data;
    p->mVertexSample.x = v[0];
    p->mVertexSample.y = v[1];
    p->mVertices.append((const char*)&p->mVertexSample,
                        sizeof(VertexData));
    ++p->mTessVertexCount;
}
void ShaderSource::tessEnd(void *self) {
    ShaderSource * p = (ShaderSource*)self;
    if (p->mTessBroken) return;

    int i, j, k;
    int n = p->mTessVertexBase;
    int more = 1;
    if (p->mTessType == GL_TRIANGLES)
        more = 4*p->mTessVertexCount/3;

    int * buf = (int*)malloc(sizeof(int)*(p->mTessVertexCount+more));
    if (!buf) {
        p->mTessBroken = true;
        return;
    }

    switch(p->mTessType) {
    case GL_TRIANGLES:
        for(i = 0; i+2 < p->mTessVertexCount; i+=3) {
            buf[4*i]   = n+3*i;
            buf[4*i+1] = n+3*i+1;
            buf[4*i+2] = n+3*i+2;
            buf[4*i+3] = PRIM_RESTART;
        }
        break;
    case GL_TRIANGLE_FAN:
        buf[0] = n;
        i = 1; j = 1; k = p->mTessVertexCount-1;
        do {
            buf[i++] = j++;
            buf[i++] = k--;
        } while(j < k);
        if (j == k)
            buf[i++] = j;
        buf[i] = PRIM_RESTART;
        break;
    case GL_TRIANGLE_STRIP:
        for(i = 0; i < p->mTessVertexCount; i++)
            buf[i] = n+i;
        buf[n+p->mTessVertexCount] = PRIM_RESTART;
        break;
    default:
        free(buf);
        return;
    }
    p->mIndices.append((const char*)buf,
                       sizeof(int)*(p->mTessVertexCount+more));
    free(buf);
}
void ShaderSource::tessError(GLenum *__errno_location(), void *self) {
    Q_UNUSED(__errno_location);
    ShaderSource * p = (ShaderSource*)self;
    p->mTessBroken = true;
    p->mVertices.resize(p->mTessFirstVertex);
    p->mIndices.resize(p->mTessFirstIndex);
}
void ShaderSource::tessCombineData(GLdouble coords[],
                              void *vertex_data[],
                              GLfloat weight[],
                              void **outData,
                              void *self)
{
    Q_UNUSED(vertex_data);
    Q_UNUSED(weight);
    ShaderSource * p = (ShaderSource*)self;
    double * cd = (double*)p->mTessBuffer.data();
    int n = p->mPathExt;
    assert((3*n+3)*sizeof(double) < uint(p->mTessBuffer.size()));
    cd[n+0] = coords[0];
    cd[n+1] = coords[1];
    cd[n+2] = coords[2];
    *outData = cd+n;
    p->mPathExt++;
}

void ShaderSource::setTexture(GLuint texture) {

}
void ShaderSource::setA(float a)
{
    mVertexSample.fill_a = qBound(0.0f, a, 1.0f);
}
void ShaderSource::setRGB(float r, float g, float b)
{
    mVertexSample.fill_r = qBound(0.0f, r, 1.0f);
    mVertexSample.fill_g = qBound(0.0f, g, 1.0f);
    mVertexSample.fill_b = qBound(0.0f, b, 1.0f);
}
void ShaderSource::setRGBA(float r, float g, float b, float a)
{
    mVertexSample.fill_r = qBound(0.0f, r, 1.0f);
    mVertexSample.fill_g = qBound(0.0f, g, 1.0f);
    mVertexSample.fill_b = qBound(0.0f, b, 1.0f);
    mVertexSample.fill_a = qBound(0.0f, a, 1.0f);
}




}setRGB(float r, float g, float b);
void setRGBA(float r, float g, float b, float a);
void setA(float a);

void setTexture(GLuint texture);
void setGradient(GLuint gtexture,
                 float start,
                 float stop,
                 bool radial);


void tessBegin(GLenum type, void *self);
void tessVertex(void * vertex_data, void * self);
void tessEnd(void * self);
void tessCombineData(GLdouble coords[3], void *vertex_data[4],
                            GLfloat weight[4], void **outData,
                            void * self);
void tessError(GLenum errno, void *self);


/*
 * edge1, edge2 and edge3 are the edge parameters
 * for each edge E{1,2,3},
 *    E.xy is 2D vector perpendicular to the edge
 *         it's length should be 1/2*line_width
 *         ti should look to the left of the vector
 *
 *    E.z  is a position factor, it's value is
 *         E.z = -P.x*E.x-P.y*E.y for each P
 *         at the edge. Start or end point of
 *         the edge can be taken
 *
 *    E1.w joint factors. set to zero to draw line ending
 *    E3.w if non-zero, sign shows turn direction of join
 *         <0 if turn was left, >0 if turn was right
 *         magnitude = 2.0 + miter factor
 *
 *    E2.w set to <zero to draw rounded lines
 *         set to >zero to draw miter joined or bewel lines
 *
 * fill1 and fill2 define a rule to obtain source color
 *    fill1.xyz and fill2.xyz make a matrix to transform
 *    pixel coordinates to texel coordinates in
 *    source 2D-texture.
 *
 *    also, s coordinate in gradient texture
 *    we will use ecomb.w to select between
 *    linear(ecomb.w>=0) and radial gradients(ecomb.w<0).
 *    clamp(s*ecomb.w,0,1) is a gradient parameter
 *    and fill1.w-fill2.w define the color interval
 *    in the gradient texture
 */
#define N "\n"
static const char * __vertex = {
    "#version 150"                      N
    "in vec4 vertex;"                   N
    "in vec4 edge1;"                    N
    "in vec4 edge2;"                    N
    "in vec4 edge3;"                    N
    "in vec4 fill1;"                    N
    "in vec4 fill2;"                    N
    "in vec4 fill3;"                    N
    "in vec4 trX;"                      N
    "in vec4 trY;"                      N
    "uniform vec2 display;"             N
                                        N
    "varying out vec2 PP;"              N
    "flat out vec4 E1;"                 N
    "flat out vec4 E2;"                 N
    "flat out vec4 E3;"                 N
    "flat out vec4 F1;"                 N
    "flat out vec4 F2;"                 N
    "flat out vec4 F3;"                 N
                                        N
    "void main {"                       N
    "  float x = dot(vertex.xy, trX.xy)+trX.z;" N
    "  float y = dot(vertex.xy, trY.xy)+trY.z;" N
    "  PP = vertex.xy;"                 N
    "  x = 2.0*display.x*x-1.0;"        N
    "  y = 2.0*display.y*y-1.0;"        N
    "  gl_Position = vec4(x,y,vertex.z,1);" N
    "  E1 = edge1;"                     N
    "  E2 = edge2;"                     N
    "  E3 = edge3;"                     N
    "  F1 = fill1;"                     N
    "  F2 = fill2;"                     N
    "  F3 = fill3;"                     N
    "}"                                 N

};
static const char * __fragment = {
   "#version 150"                       N
    "varying in vec2 PP;"               N
    "flat in vec4 E1;"                  N
    "flat in vec4 E2;"                  N
    "flat in vec4 E3;"                  N
    "flat in vec4 F1;"                  N
    "flat in vec4 F2;"                  N
    "flat in vec4 F3;"                  N
    "uniform sampler2D source2;"        N
    "uniform sampler1D source1;"        N
                                        N
    "out color;"                        N
                                        N
    "void main {"                       N

    //calculate edge parameters
    "vec3 ps=vec3(E1.z,E2.z,E3.z);"          N
    "     ps+=vec3(E1.x,E2.x,E3.x)*PP.x;"    N
    "     ps+=vec3(E1.y,E2.y,E3.y)*PP.y;"    N
    "vec2 tf=vec3(E1.w,E3.w);"               N

    //E1.w=0 means edge is over at start
    //E3.w=0 means edge is over at end
    //Ex.w in range [1..3] are translated to [-1,1]
    //and are interpreted as line center for joints

    //this is a special parameter set for line endings
    "vec4 eg;"                               N
    "     eg.x=max(ps.x,0.0);"               N
    "     eg.x=abs(ps.y);"                   N
    "     eg.x=abs(ps.y);"                   N
    "     eg.x=max(ps.z,0.0);"               N

    //this a parameter set for line joints
    "vec4 sg=sign(tf.xxyy);"                 N
    "vec4 pk=ps.xyyz*sg;"                    N
    "vec4 pf=abs(tf.xxyy)"                   N
    "     pf=(pk-pf+vec4(2.0)) /"            N
    "        (vec4(3.0)-pf);"                N

    //now we can choose shading function
    "     pf = mix(eg,pf,abs(sg));"          N

    //E2.w will select between rounded joints
    //and miter joints
    "     float jf = step(0.0,E2.w);"        N
    "     pf*= mix(pf, vec4(jf), jf);"       N
    "     float A = clamp(F3.a*"             N
    "                     max(pf.x+pf.y,"    N
    "                         pf.z+pf.w),"   N
    "                     0.0,1.0);"         N

    //now it's a time to apply fill rule
    //first we determine the fill direction
    "  float fs = dot(PP,F1.xy)+F1.z;"       N
    "  float ft = dot(PP,F2.xy)+F2.z;"       N

    //read from source texture
    //will return RGBA(0,0,0,1) if there is no texture
    "  vec4  tx = texture(source2,vec2(fs,ft));" N

    //calculate gradient parameter
    "  float fg = sqrt(fs*fs+ft*ft);"        N

    //read from source texture
    //will return RGBA(0,0,0,1) if there is no texture
    "  vec4  tx2= texture(source1,fg));"     N

    //calculate the resulting color
    "  vec4  out = vec4(tx.rgb+"             N
    "                   tx2.rgb+"            N
    "                   F3.rgb,"             N
    "                   tx.a*tx2.a*A);"      N
    "  color = clamp(out,0.0,1.0);"          N
    "}"                                      N
};
