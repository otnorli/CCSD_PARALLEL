TEMPLATE = app
CONFIG += console
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp \
    fill_alpha.cpp \
    matrix_size_setter.cpp \
    hartree_fock_solver.cpp

HEADERS += \
    fill_alpha.h \
    matrix_size_setter.h \
    hartree_fock_solver.h

