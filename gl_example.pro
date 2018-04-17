#-------------------------------------------------
#
# Project created by QtCreator 2018-04-01T19:01:48
#
#-------------------------------------------------

QT       += core gui

greaterThan(QT_MAJOR_VERSION, 4): QT += widgets

TARGET = gl_example
TEMPLATE = app

# The following define makes your compiler emit warnings if you use
# any feature of Qt which has been marked as deprecated (the exact warnings
# depend on your compiler). Please consult the documentation of the
# deprecated API in order to know how to port your code away from it.
DEFINES += QT_DEPRECATED_WARNINGS

# You can also make your code fail to compile if you use deprecated APIs.
# In order to do so, uncomment the following line.
# You can also select to disable deprecated APIs only up to a certain version of Qt.
#DEFINES += QT_DISABLE_DEPRECATED_BEFORE=0x060000    # disables all the APIs deprecated before Qt 6.0.0


SOURCES += \
    main.cpp \
    painter_base.cpp \
    painter_destination.cpp \
    painter_source.cpp \
    testwindow.cpp \
    painter_example.cpp \
    painter_destination_qt.cpp

HEADERS += \
    painter_base.h \
    painter_destination.h \
    painter_source.h \
    testwindow.h \
    painter_example.h \
    painter_destination_qt.h

win32: LIBS += -lGLU32
else: LIBS += -lGLU
LIBS += $$PWD/libtess2/Build/libtess2.a

RESOURCES += \
    resources.qrc

DISTFILES += \
    notmyfault.txt
