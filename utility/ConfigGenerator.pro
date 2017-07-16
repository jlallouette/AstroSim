TEMPLATE = app
TARGET = ConfigGenerator
QT += core gui widgets
HEADERS += ListParamGroup.h \
    ListValueParam.h \
    ParamGroup.h \
    SingleValueParam.h \
    SingleParam.h \
    ParamDescription.h \
    ConfigGenerator.h
SOURCES += ListParamGroup.cpp \
    ListValueParam.cpp \
    ParamGroup.cpp \
    SingleValueParam.cpp \
    SingleParam.cpp \
    ParamDescription.cpp \
    main.cpp \
    ConfigGenerator.cpp
FORMS += ConfigGenerator.ui
RESOURCES += 
