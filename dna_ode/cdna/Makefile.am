noinst_LTLIBRARIES = libcdna.la
libcdna_la_SOURCES = cdna.c cdna.h \
			     cdna_init_functions.c \
			     cdna_objects_creation.c \
			     cdna_fin_functions.c \
			     cdna_dynamics.c \
			     cdna_output_functions.c \
			     cdna_print_functions.c\
			     cdna_joints_solution.c

libcdna_la_CPPFLAGS = $(ODE_CFLAGS) -I$(top_srcdir)/dna_ode/core -I$(top_srcdir)/dna_ode/src -I$(top_srcdir)/dna_ode
libcdna_la_LIBADD = $(top_builddir)/dna_ode/core/libdna_ode_core.la -L$(top_builddir)/dna_ode/core

AM_LDFLAGS = $(GSL_LIBS) $(ODE_LIBS) $(CONFIG_LIBS)
