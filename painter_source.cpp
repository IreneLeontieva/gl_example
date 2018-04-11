#include "painter_source.h"
#include <stdexcept>
#include <QOpenGLShader>
#include <cmath>
#include <assert.h>
#include <GL/glu.h>
#include <QDebug>

static QString __vertex();
static QString __fragment();

#ifdef _WIN32
typedef void (CALLBACK *GLU_CALLBACK_T)();
#else
#define GLU_CALLBACK_T _GLUfuncptr
#endif

#define N "\n"
static const char * __shader_vertex =
/*01*/    "#version 150"                              N
/*02*/    "in vec4 vertex;"                           N
/*03*/    "in vec4 edge1;"                            N
/*04*/    "in vec4 edge2;"                            N
/*05*/    "in vec4 edge3;"                            N
/*06*/    "in vec4 fill1;"                            N
/*07*/    "in vec4 fill2;"                            N
/*08*/    "in vec4 fill3;"                            N
/*09*/    "in vec4 trX;"                              N
/*10*/    "in vec4 trY;"                              N
/*11*/    "uniform vec4 display;"                     N
/*12*/    "out vec2 PP;"                              N
/*13*/    "flat out vec4 E1;"                         N
/*14*/    "flat out vec4 E2;"                         N
/*15*/    "flat out vec4 E3;"                         N
/*17*/    "flat out vec4 F1;"                         N
/*18*/    "flat out vec4 F2;"                         N
/*19*/    "flat out vec4 F3;"                         N
/*20*/    "void main() {"                             N
/*21*/    "  float x = dot(vertex.xy, trX.xy)+trX.z;" N
/*22*/    "  float y = dot(vertex.xy, trY.xy)+trY.z;" N
/*23*/    "  PP = vertex.xy;"                         N
/*24*/    "  x = display.x*x+display.z;"              N
/*25*/    "  y = display.y*y+display.w;"              N
/*26*/    "  gl_Position = vec4(x,y,vertex.z,1.0);"   N
/*27*/    "  E1 = edge1;"                             N
/*28*/    "  E2 = edge2;"                             N
/*29*/    "  E3 = edge3;"                             N
/*30*/    "  F1 = fill1;"                             N
/*31*/    "  F2 = fill2;"                             N
/*32*/    "  F3 = fill3;"                             N
/*33*/    "}";

static const char * __shader_fragment =
/*01*/    "#version 150"                              N
/*02*/    "in vec2 PP;"                               N
/*03*/    "flat in vec4 E1;"                          N
/*04*/    "flat in vec4 E2;"                          N
/*05*/    "flat in vec4 E3;"                          N
/*06*/    "flat in vec4 F1;"                          N
/*07*/    "flat in vec4 F2;"                          N
/*08*/    "flat in vec4 F3;"                          N
/*09*/    "uniform sampler2D source2;"                N
/*10*/    "uniform sampler1D source1;"                N
/*11*/                                                N
/*12*/    "out vec4 color;"                           N
/*13*/                                                N
/*14*/    "void main() {"                             N
//calculate edge parameters                           N
/*15*/    "vec3 ps=vec3(E1.z,E2.z,E3.z);"             N
/*16*/    "     ps+=vec3(E1.x,E2.x,E3.x)*PP.x;"       N
/*17*/    "     ps+=vec3(E1.y,E2.y,E3.y)*PP.y;"       N
/*18*/    "vec2 tf=vec2(E1.w,E3.w);"                  N
//E1.w=0 means edge is over at start                  N
//E3.w=0 means edge is over at end                    N
//Ex.w in range [1..3] are translated to [-1,1]       N
//and are interpreted as line center for joints       N
//this is a special parameter set for line endings    N
/*19*/    "vec4 eg=vec4(max(ps.x,0.0),"               N
/*20*/    "             abs(ps.y),"                   N
/*21*/    "             abs(ps.y),"                   N
/*22*/    "             max(ps.z,0.0));"              N
//this a parameter set for line joints
/*23*/    "vec4 sg=sign(tf.xxyy);"                    N

/*24*/    "vec4 pk=ps.xyyz*sg;"                       N
/*25*/    "vec4 pf=abs(tf.xxyy)-vec4(2.0);"           N
/*26*/    "     pf=(pk-pf) /"                         N
/*27*/    "        (vec4(1.0)-pf);"                   N
          "     pf=max(pf, vec4(0.0));"
//now we can choose shading function                  N
/*28*/    "     pf = mix(eg,pf,abs(sg));"             N
//E2.w will select between rounded joints             N
//and miter joints                                    N
/*29*/    "     float jf = E2.w;"                     N
          "     float A1 = max(pf.x+pf.y,pf.z+pf.w);" N
          "     pf *= pf;"                            N
          "     float A2 = max(sqrt(pf.x+pf.y),sqrt(pf.z+pf.w));" N
/*30*/    "     float A= mix(A2, A1*jf, step(0.0, jf));" N
          "     A = step(A, 1.0);"
//now it's a time to apply fill rule                  N
//first we determine the fill direction               N
/*35*/    "  float fs = dot(PP,F1.xy)+F1.z;"          N
/*36*/    "  float ft = dot(PP,F2.xy)+F2.z;"          N
//read from source texture                            N
//will return RGBA(0,0,0,1) if there is no texture    N
/*37*/    "  vec4  tx = texture(source2,vec2(fs,ft));" N
//calculate gradient parameter                        N
/*38*/    "  float fg = sqrt(fs*fs+ft*ft);"           N
          "  fg = mix(fs,fg,F1.w);"
//read from source texture                            N
//will return RGBA(0,0,0,1) if there is no texture    N
/*39*/    "  vec4  tx2= texture(source1,fg);"         N
//calculate the resulting color                       N
/*40*/    "  vec4  o = vec4(tx.rgb+"                  N
/*41*/    "                 tx2.rgb+"                 N
/*42*/    "                 F3.rgb,"                   N
/*43*/    "                 tx.a*tx2.a*A*F3.a);"      N
/*44*/    "  color = clamp(o,0.0,1.0);"               N
/*45*/    "}"  ;


PainterSource::PainterSource()
{
    qDebug()<<"PainterSource constructor started";
    if (!PainterBase::beginPaint())
        throw std::runtime_error("cannot bind context");
    auto $ = PainterBase::gl();
    try {
        while($->glGetError());

        mProgram = new QOpenGLShaderProgram();
        bool ok;
        ok = mProgram->addShaderFromSourceCode(QOpenGLShader::Vertex,
                                               __vertex());
        if (!ok) {
            qDebug()<<mProgram->log();
            throw std::runtime_error("failed to compile vertex shader");
        }
        ok = mProgram->addShaderFromSourceCode(QOpenGLShader::Fragment,
                                               __fragment());
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
            qDebug()<<">"<<((const char*)gluErrorString(e));
            throw std::runtime_error("pending GL errors(shaders)");
        }

        $->glPrimitiveRestartIndex(PRIM_RESTART);

        $->glGenVertexArrays(1, &mVAO);
        $->glGenBuffers(1, &mVBO);
        $->glGenBuffers(1, &mIBO);
        e = $->glGetError();
        if (e) {
            qDebug()<<">"<<((const char*)gluErrorString(e));
            throw std::runtime_error("pending GL errors(gen objects)");
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
                    throw std::runtime_error("program attribute not found");
                }
                $->glVertexAttribPointer(a0, 4, GL_FLOAT, GL_FALSE,
                                         16*9, (GLvoid*)(long)(16*j));
                $->glEnableVertexAttribArray(a0);
            }
            e = $->glGetError();
            if (e) {
                qDebug()<<">"<<((const char*)gluErrorString(e));
                throw std::runtime_error("pending GL errors(vao)");
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
        gluTessProperty(mTesselator, GLU_TESS_WINDING_RULE, GLU_TESS_WINDING_NONZERO);
        gluTessCallback(mTesselator, GLU_TESS_BEGIN_DATA,
                        (GLU_CALLBACK_T)tessBegin);
        gluTessCallback(mTesselator, GLU_TESS_VERTEX_DATA,
                        (GLU_CALLBACK_T)tessVertex);
        gluTessCallback(mTesselator, GLU_TESS_COMBINE_DATA,
                        (GLU_CALLBACK_T)tessCombineData);
        gluTessCallback(mTesselator, GLU_TESS_END_DATA,
                        (GLU_CALLBACK_T)tessEnd);
        gluTessCallback(mTesselator, GLU_TESS_ERROR_DATA,
                        (GLU_CALLBACK_T)tessError);

        memset(&mVertexSample, 0, sizeof(mVertexSample));
        //clear color
        mVertexSample.fill_a = 1;
        //clear transform
        mVertexSample.mat_xx = 1;
        mVertexSample.mat_yy = -1.0;
        mVertexSample.mat_y0 = 300;
        mVertexSample._y0 =1.0;
    } catch(...) {
        if (mVAO) $->glDeleteVertexArrays(1, &mVAO);
        if (mVBO) $->glDeleteBuffers(1, &mVBO);
        if (mIBO) $->glDeleteBuffers(1, &mIBO);
        if (mTesselator)
            gluDeleteTess(mTesselator);
        delete mProgram;
        throw;
    }
    PainterBase::endPaint();
    qDebug()<<"PainterSource constructed";
}
PainterSource::~PainterSource()
{
    PainterBase::beginPaint();
    auto $ = PainterBase::gl();
    if (mVAO) $->glDeleteVertexArrays(1, &mVAO);
    if (mVBO) $->glDeleteBuffers(1, &mVBO);
    if (mIBO) $->glDeleteBuffers(1, &mIBO);
    if (mTesselator)
        gluDeleteTess(mTesselator);
    delete mProgram;
    PainterBase::endPaint();
}
bool PainterSource::paint(const QSize& size, bool wireframe) {
    if (size.isEmpty())
        return false;

    if (mIndices.size() && mVertices.size() && mRCS.size())
    {
        auto $ = PainterBase::gl();
        if (mUpdated) {
            $->glBindBuffer(GL_ARRAY_BUFFER, mVBO);
            $->glBufferData(GL_ARRAY_BUFFER, mVertices.size(), mVertices.constData(), GL_STREAM_DRAW);
            $->glBindBuffer(GL_ARRAY_BUFFER, 0);
            $->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mIBO);
            $->glBufferData(GL_ELEMENT_ARRAY_BUFFER, mIndices.size(), mIndices.constData(), GL_STREAM_DRAW);
            $->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);
            mUpdated = false;
        }

        if (!mProgram->bind())
            return false;
        $->glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        $->glBlendEquation(GL_FUNC_ADD);
        $->glDepthFunc(GL_GEQUAL);

        if (!wireframe) {
            $->glEnable(GL_BLEND);
            $->glEnable(GL_DEPTH_TEST);
        }
        $->glEnable(GL_PRIMITIVE_RESTART);
        $->glPolygonMode(GL_FRONT_AND_BACK, wireframe ?
                             GL_LINE : GL_FILL);

        mProgram->setUniformValue(mPaintDPS,
                                  0.005f, 0.005f,
                                  -1.0f,
                                  -1.0f);
                                  //200.0f/float(size.width()),
                                  //200.0f/float(size.height()));
        $->glBindVertexArray(mVAO);
        for(int i = 0; i < mRCS.size(); ++i) {
            const RenderCommand& rc = mRCS.at(i);
            $->glActiveTexture(GL_TEXTURE0);
            $->glBindTexture(GL_TEXTURE_1D, rc.mTex1);
            $->glActiveTexture(GL_TEXTURE1);
            $->glBindTexture(GL_TEXTURE_2D, rc.mTex2);
            $->glDrawElements(GL_TRIANGLE_STRIP,
                              rc.mNumIndices,
                              GL_UNSIGNED_INT,
                              (GLvoid*)(long)(rc.mFirstIndex*4));
        }
        $->glBindVertexArray(0);
        $->glDisable(GL_BLEND);
        $->glDisable(GL_PRIMITIVE_RESTART);
        $->glDisable(GL_DEPTH_TEST);
        mProgram->release();
    }
    return true;
}
void PainterSource::clear() {
    mIndices.clear();
    mVertices.clear();
    mRCS.clear();
    mDepthIndex = 0.0f;
    mUpdated = true;
}
void PainterSource::setPath(const double *coords, int ncoords, bool closed) {
    mTessBuffer.clear();
    mPathSize   = ncoords;
    mPathExt    = ncoords;
    mPathClosed = closed;

    mTessBuffer.resize(3*ncoords*2*sizeof(double));
    mTessBuffer.detach();
    double * dst = (double*)mTessBuffer.data();
    for(int i = 0; i < ncoords; ++i) {
        dst[3*i]   = coords[2*i];
        dst[3*i+1] = coords[2*i+1];
        dst[3*i+2] = coords[2*i];
    }
}


static inline double TURN(double dx1, double dy1, double dx2, double dy2, double miter) {
    double d = dx2*dy1-dx1*dy2;
    if (d == 0.0) d = 1.0;
    return std::copysign(2.0+miter, d);
}
static inline bool EXTENT(double &ox, double &oy,
                          double dx1, double dy1,
                          double dx2, double dy2)
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
     * we detect it by negatic c which is cos(phi)
     */
    double c = dx1*dx2+dy1*dy2;
    double f = c < 0.0 ? -1.0 : 1.0;
    //FIXME: handle this case better
    double k = 1.0/(1.0+fabs(c));
    double avx = dx1 + dx2*f;
    double avy = dy1 + dy2*f;
    ox = avx*k;
    oy = avy*k;
    return (c < 0.0);
}
#define EMIT_VERTEX(v, n) \
    {\
        int val = n;\
        mVertices.append((const char*)&v, sizeof(v)); \
        mIndices.append((const char*)&val, sizeof(val));\
    }
#define DEPTH_INCREMENT (1.0f/float(1<<22))

void PainterSource::strokePath() {
    if (mPathSize < 2) return;
    assert(mPathExt >= mPathSize);
    assert(uint(mTessBuffer.size()) >= mPathExt*3*sizeof(double));

    mVertexSample.z = mDepthIndex;
    mDepthIndex += DEPTH_INCREMENT;



    //save info for case of error
    mTessBroken     = false;
    mTessFirstIndex = mIndices.size();
    mTessFirstVertex= mVertices.size();
    bool  reallyClosed = mPathClosed & (mPathSize > 2);
    int   indexCount   = mVertices.size()/sizeof(mVertexSample);

    const double * src = (const double*)mTessBuffer.data();
    const double * A, * B, * C, * D, * M;

    A = nullptr; B = src; C = src+3; D = nullptr;
    M = src + mPathSize*3;

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
    float wc = 1.0f/qMax(1.0f, mPenW);
    mVertexSample.e1x = -dx*wc;
    mVertexSample.e1y = -dy*wc;
    mVertexSample.e1z = (B[0]*dx+B[1]*dy)*wc;
    mVertexSample.e1w = 0.0;
    mVertexSample.e2x = -dy*wc;
    mVertexSample.e2y = dx*wc;
    mVertexSample.e2z = (B[0]*dy-B[1]*dx)*wc;
    mVertexSample.e2w = -1.0;
    if (!mRounded) mVertexSample.e2w = mMiter;
    mVertexSample.e3x = dx*wc;
    mVertexSample.e3y = dy*wc;
    mVertexSample.e3z = (-C[0]*dx-C[1]*dy)*wc;
    mVertexSample.e3w = 0.0;
    if (reallyClosed) {
        A = M-3;
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
        mVertexSample.e1x = -day*wc;
        mVertexSample.e1y =  dax*wc;
        mVertexSample.e1z = (B[0]*day-B[1]*dax)*wc;
        mVertexSample.e1w = TURN(dax,day,dx,dy,mMiter);
    }

    bool first = true;
    do {
        //qDebug("ITER: B= %f %f - C= %f %f", B[0],B[1],C[0],C[1]);
        if (D) {
            //qDebug("ITER: D = %f %f",D[0],D[1]);
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
            mVertexSample.e3x = -ddy*wc;
            mVertexSample.e3y =  ddx*wc;
            mVertexSample.e3z = (C[0]*ddy-C[1]*ddx)*wc;
            mVertexSample.e3w = TURN(dx,dy,ddx,ddy,mMiter);
        }
        //emit vertex
        double offset_x, offset_y, k;
        if (first) {
            //extent status does not matter here, because we
            //at start of the line
            EXTENT(offset_x, offset_y, dax, day, dx, dy);
            k = reallyClosed ? 0.0 : mPenW;
            mVertexSample.x = B[0] - mPenW*offset_y - dx*k;
            mVertexSample.y = B[1] + mPenW*offset_x - dy*k;
            EMIT_VERTEX(mVertexSample, indexCount++);
            //qDebug("VERTEXS %f %f", mVertexSample.x, mVertexSample.y);
            mVertexSample.x = B[0] + mPenW*offset_y - dx*k;
            mVertexSample.y = B[1] - mPenW*offset_x - dy*k;
            EMIT_VERTEX(mVertexSample, indexCount++);
            //qDebug("VERTEXS %f %f", mVertexSample.x, mVertexSample.y);
            first = false;
        }
        bool bad = EXTENT(offset_x, offset_y, dx, dy, ddx, ddy);
        //qDebug("EXTENT %f %f -- %f %f -- %f %f", offset_x, offset_y, dx, dy, ddx, ddy);
        k = D ? 0.0 : mPenW;
        //qDebug("K = %f", k);
        mVertexSample.x = C[0] - mPenW*offset_y + ddx*k;
        mVertexSample.y = C[1] + mPenW*offset_x + ddy*k;
        EMIT_VERTEX(mVertexSample, indexCount++);
        //qDebug("VERTEX  %f %f", mVertexSample.x, mVertexSample.y);
        mVertexSample.x = C[0] + mPenW*offset_y + ddx*k;
        mVertexSample.y = C[1] - mPenW*offset_x + ddy*k;
        EMIT_VERTEX(mVertexSample, indexCount++);
        //qDebug("parms = %f %f %f; %f %f %f; %f %f %f;",
        //       mVertexSample.e1x/wc,mVertexSample.e1y/wc,mVertexSample.e1z/wc,
        //       mVertexSample.e2x/wc,mVertexSample.e2y/wc,mVertexSample.e2z/wc,
        //       mVertexSample.e3x/wc,mVertexSample.e3y/wc,mVertexSample.e3z/wc);
        //qDebug("VERTEX  %f %f", mVertexSample.x, mVertexSample.y);
        if (bad) {
            auto scopy = mVertexSample;
            //set line ending feature
            scopy.e3x = dx*wc;
            scopy.e3y = dy*wc;
            scopy.e3z = (-C[0]*dx-C[1]*dy)*wc;
            scopy.e3w = 0.0;
            //add line cap
            scopy.x = C[0] - mPenW*(offset_y - offset_x);
            scopy.y = C[1] + mPenW*(offset_x + offset_y);
            EMIT_VERTEX(scopy, indexCount++);
            //qDebug("VERTEXB %f %f", scopy.x, scopy.y);
            scopy.x = C[0] + mPenW*(offset_y + offset_x);
            scopy.y = C[1] - mPenW*(offset_x - offset_y);
            EMIT_VERTEX(scopy, indexCount++);
            //qDebug("VERTEXB %f %f", scopy.x, scopy.y);
            //restart line sequence
            int indices[3] = { PRIM_RESTART, indexCount-3, indexCount-4 };
            mIndices.append((const char*)indices, sizeof(indices));
        }

        //if we are rasterizing a line, just stop here
        if (!D || C == src) break;

        //iterate to next vertex
        A = B;
        B = C;
        C = D;
        dx = ddx;
        dy = ddy;
        dl = ddl;
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
            mVertexSample.e3x = ddx*wc;
            mVertexSample.e3y = ddy*wc;
            mVertexSample.e3z = (-C[0]*ddx-C[1]*ddy)*wc;
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
void PainterSource::fillPath() {
    if (!mPathClosed) return;
    if (mPathSize < 3) return;
    assert(mPathExt >= mPathSize);
    assert(uint(mTessBuffer.size()) >= mPathExt*3*sizeof(double));

    mVertexSample.z = mDepthIndex;
    mDepthIndex += DEPTH_INCREMENT;
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
        qDebug("vertex in %f %f", dst[i*3], dst[i*3+1]);
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
        mVertices.resize(mTessFirstVertex);
    }
}
//====utility functions=============================
void PainterSource::insertRCS(int prev) {
    int ni = mIndices.size();

    int extra_restart = PRIM_RESTART;
    mIndices.append((const char*)&extra_restart, sizeof(extra_restart));
    mUpdated = true;
    /*do {//check if RCS can be combined
        if (mRCS.isEmpty())
            break;
        auto last = mRCS.last();
        if (last.mTex1 != mTex1 || last.mTex2 != mTex2)
            break;

        mRCS.last().mNumIndices = 1+ni/4-mRCS.last().mFirstIndex;

        return;
    } while(0);*/
    RenderCommand rc;
    rc.mFirstIndex = prev/4;
    rc.mNumIndices = (ni - prev)/4;
    rc.mTex1 = mTex1;
    rc.mTex2 = mTex2;
    mRCS.append(rc);
}
//====tesselator=================================
void PainterSource::tessBegin(GLenum type, void *self) {
    PainterSource * p = (PainterSource*)self;
    if (p->mTessBroken) return;
    p->mTessType        = type;
    p->mTessVertexBase  = p->mVertices.size()/
                          sizeof(VertexData);
    p->mTessVertexCount = 0;
    qDebug("tessBegin %s", type == GL_TRIANGLE_FAN ? "fan" :
                           type == GL_TRIANGLE_STRIP ? "strip" :
                                                       "other");
}
void PainterSource::tessVertex(void *vertex_data, void *self) {
    PainterSource * p = (PainterSource*)self;
    if (p->mTessBroken) return;

    double * v = (double*)vertex_data;
    p->mVertexSample.x = v[0];
    p->mVertexSample.y = v[1];
    p->mVertices.append((const char*)&p->mVertexSample,
                        sizeof(VertexData));
    ++p->mTessVertexCount;

    if (p->mTessBuffer.data() < (char*)vertex_data)
        qDebug("vertex below");
    if (p->mTessBuffer.data()+p->mTessBuffer.size() >= (char*)vertex_data)
        qDebug("vertex above");

    qDebug("tessVertex %f %f %f", v[0], v[1], v[2]);
}
void PainterSource::tessEnd(void *self) {
    PainterSource * p = (PainterSource*)self;
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
    qDebug("tessEnd");

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
void PainterSource::tessError(GLenum *__errno_location(), void *self) {
    Q_UNUSED(__errno_location);
    PainterSource * p = (PainterSource*)self;
    p->mTessBroken = true;
    p->mVertices.resize(p->mTessFirstVertex);
    p->mIndices.resize(p->mTessFirstIndex);
    qDebug("tessEnd");

}
void PainterSource::tessCombineData(GLdouble coords[],
                              void *vertex_data[],
                              GLfloat weight[],
                              void **outData,
                              void *self)
{
    Q_UNUSED(vertex_data);
    Q_UNUSED(weight);
    PainterSource * p = (PainterSource*)self;
    double * cd = (double*)p->mTessBuffer.data();
    int n = p->mPathExt;
    assert((3*n+3)*sizeof(double) < uint(p->mTessBuffer.size()));
    cd[n+0] = coords[0];
    cd[n+1] = coords[1];
    cd[n+2] = coords[2];
    *outData = cd+n;
    p->mPathExt++;
    qDebug("tessCombine!");
}
void PainterSource::setLinearGradient(GLuint gtexture,
                                     float startX,
                                     float startY,
                                     float stopX,
                                     float stopY)
{
    mTex1 = gtexture;
    mTex2 = 0;
    mVertexSample.fill_r = 0.0;
    mVertexSample.fill_g = 0.0;
    mVertexSample.fill_b = 0.0;

    float dx = stopX - startX;
    float dy = stopY - startY;
    float len= sqrt(dx*dx+dy*dy);
    if (len > 0.000001f) {
        len = 1.0f/len;
        dx *= len; dy *= len;

        mVertexSample.fill_xx= dx;
        mVertexSample.fill_xy= dy;
        mVertexSample.fill_x0= -dx*startX-dy*startY;//to pin start to 0
        mVertexSample.fill_yx= 0.0;
        mVertexSample.fill_yy= 0.0;
        mVertexSample.fill_y0= 0.0;
    } else {
        mVertexSample.fill_xx= 0.0;
        mVertexSample.fill_xy= 0.0;
        mVertexSample.fill_x0= 1.0;
        mVertexSample.fill_yx= 0.0;
        mVertexSample.fill_yy= 0.0;
        mVertexSample.fill_y0= 0.0;
    }
    mVertexSample.fill_G1 = 0.0;
}
void PainterSource::setRadialGradient(GLuint gtexture,
                                     float centerX,
                                     float centerY,
                                     float radius)
{
    mTex1 = gtexture;
    mTex2 = 0;
    mVertexSample.fill_r = 0.0;
    mVertexSample.fill_g = 0.0;
    mVertexSample.fill_b = 0.0;

    if (radius > 0.000001f) {
        float k = 1.0f/radius;

        mVertexSample.fill_xx= k;
        mVertexSample.fill_xy= 0.0;
        mVertexSample.fill_x0= -k*centerX;//to pin start to 0
        mVertexSample.fill_yx= 0.0;
        mVertexSample.fill_yy= k;
        mVertexSample.fill_y0= -k*centerY;
    } else {
        mVertexSample.fill_xx= 0.0;
        mVertexSample.fill_xy= 0.0;
        mVertexSample.fill_x0= 1.0;
        mVertexSample.fill_yx= 0.0;
        mVertexSample.fill_yy= 0.0;
        mVertexSample.fill_y0= 0.0;
    }
    mVertexSample.fill_G1 = 1.0;
}
void PainterSource::setTexture(GLuint texture, float origX, float origY, float scale) {
    mTex1 = 0;
    mTex2 = texture;
    mVertexSample.fill_r = 0.0;
    mVertexSample.fill_g = 0.0;
    mVertexSample.fill_b = 0.0;
    mVertexSample.fill_xx= scale;
    mVertexSample.fill_xy= 0.0;
    mVertexSample.fill_yx= 0.0;
    mVertexSample.fill_yy= scale;
    mVertexSample.fill_x0= -origX;
    mVertexSample.fill_y0= -origY;
    mVertexSample.fill_G1 = 0.0;
}
void PainterSource::setA(float a)
{
    mVertexSample.fill_a = qBound(0.0f, a, 1.0f);
}
void PainterSource::setRGB(float r, float g, float b)
{
    mVertexSample.fill_r = qBound(0.0f, r, 1.0f);
    mVertexSample.fill_g = qBound(0.0f, g, 1.0f);
    mVertexSample.fill_b = qBound(0.0f, b, 1.0f);
    mVertexSample.fill_G1 = 0.0;
}
void PainterSource::setRGBA(float r, float g, float b, float a)
{
    mVertexSample.fill_r = qBound(0.0f, r, 1.0f);
    mVertexSample.fill_g = qBound(0.0f, g, 1.0f);
    mVertexSample.fill_b = qBound(0.0f, b, 1.0f);
    mVertexSample.fill_a = qBound(0.0f, a, 1.0f);
    mVertexSample.fill_G1 = 0.0;
}
void PainterSource::setPenSize(float pen, float miter, bool rounded)
{
    mPenW  = qMax(1.0f, pen*0.5f);
    mMiter = qBound(0.0f, miter, 1.0f);
    mRounded = rounded;
}

static QString __vertex() {
    return QString::fromLatin1(__shader_vertex);
}

static QString __fragment() {
    return QString::fromLatin1(__shader_fragment);
}
