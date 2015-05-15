noinst_LTLIBRARIES = libdna_ode_core.la
libdna_ode_core_la_SOURCES = config_read.c config_read.h \
                             langevin.c langevin.h \
                             dna_ode_math.c dna_ode_math.h \
                             ode_body_functions.c ode_body_functions.h \
                             potentials.c potentials.h \
                             external_potentials.c external_potentials.h \
                             objects_creation.c objects_creation.h \
                             random.c random.h \
                             mechanical_models.c mechanical_models.h \
                             dna.c dna.h \
                             bead.c bead.h \
                             writhe.c writhe.h \
                             trajectory_tools.c trajectory_tools.h \
                             utils.c utils.h\
                             curves.c curves.h \
                             dna_ode_core.h

libdna_ode_core_la_CPPFLAGS = $(ODE_CFLAGS) $(CONFIG_CFLAGS)
AM_CFLAGS = -I$(top_srcdir)/dna_ode/src
AM_LDFLAGS = $(GSL_LIBS) $(ODE_LIBS) $(CONFIG_LIBS)
