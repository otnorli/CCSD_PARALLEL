TEMPLATE = app
CONFIG += console
CONFIG -= qt

SOURCES += main.cpp \
    hartree_fock_solver.cpp \
    hartree_integrals.cpp \
    matrix_size_setter.cpp \
    fill_alpha.cpp \
    coupled_cluster_integrals.cpp \
    relax_the_structure.cpp \
    diis.cpp \
    initializer.cpp \
    coupled_cluster_v2.cpp \
    ccsd_memory_optimized.cpp \
    ccsd_v2_optimized.cpp

LIBS += -llapack -lblas

HEADERS += \
    hartree_fock_solver.h \
    hartree_integrals.h \
    matrix_size_setter.h \
    fill_alpha.h \
    coupled_cluster_integrals.h \
    relax_the_structure.h \
    diis.h \
    initializer.h \
    coupled_cluster_v2.h \
    ccsd_memory_optimized.h \
    ccsd_v2_optimized.h
