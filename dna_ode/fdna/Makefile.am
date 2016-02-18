noinst_LTLIBRARIES = libfdna.la
libfdna_la_SOURCES = fdna.c fdna.h \
		       fdna_init_functions.c \
		       fdna_objects_creation.c \
		       fdna_fin_functions.c \
		       fdna_dynamics.c \
		       fdna_output_functions.c \
		       fdna_joints_solution.c \
		       fdna_print_functions.c
		       
libfdna_la_CFLAGS = $(ODE_CFLAGS) -I$(top_srcdir)/dna_ode/core -I$(top_srcdir)/dna_ode/src -I$(top_srcdir)/dna_ode
libfdna_la_LIBADD = $(top_builddir)/dna_ode/core/libdna_ode_core.la -L$(top_builddir)/dna_ode/core

AM_LDFLAGS = $(GSL_LIBS) $(ODE_LIBS) $(CONFIG_LIBS)
