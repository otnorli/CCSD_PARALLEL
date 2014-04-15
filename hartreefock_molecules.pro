TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    hartree_fock_solver.cpp \
    hartree_integrals.cpp \
    matrix_size_setter.cpp \
    fill_alpha.cpp \
    initializer.cpp \
    ccsd_memory_optimized.cpp

LIBS += -llapack -lblas -lmpi

HEADERS += \
    hartree_fock_solver.h \
    hartree_integrals.h \
    matrix_size_setter.h \
    fill_alpha.h \
    initializer.h \
    ccsd_memory_optimized.h

# MPI Settings
QMAKE_CXX = mpicxx
QMAKE_CXX_RELEASE = $$QMAKE_CXX
QMAKE_CXX_DEBUG = $$QMAKE_CXX
QMAKE_LINK = $$QMAKE_CXX
QMAKE_CC = mpicc

QMAKE_CFLAGS += $$system(mpicc --showme:compile)
QMAKE_LFLAGS += $$system(mpicxx --showme:link)
QMAKE_CXXFLAGS += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
QMAKE_CXXFLAGS_RELEASE += $$system(mpicxx --showme:compile) -DMPICH_IGNORE_CXX_SEEK
