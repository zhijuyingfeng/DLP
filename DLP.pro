TEMPLATE = app
CONFIG += console c++11
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += \
        bigint.cpp \
        main.cpp \
        mpn.cpp

HEADERS += \
    bigint.h \
    mpn.h

QMAKE_CXXFLAGS += -O3
