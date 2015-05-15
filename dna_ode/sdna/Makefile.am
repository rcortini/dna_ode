noinst_LTLIBRARIES = libsdna.la
libsdna_la_SOURCES = sdna.c sdna.h\
			   sdna_init_functions.c\
			   sdna_objects_creation.c\
			   sdna_joints_solution.c\
			   sdna_fin_functions.c\
			   sdna_dynamics.c\
			   sdna_output_functions.c\
			   sdna_print_functions.c



libsdna_la_CPPFLAGS = $(ODE_CFLAGS) -I$(top_srcdir)/dna_ode/core -I$(top_srcdir)/dna_ode/src -I$(top_srcdir)/dna_ode
libsdna_la_LIBADD = $(top_builddir)/dna_ode/core/libdna_ode_core.la -L$(top_builddir)/dna_ode/core

AM_LDFLAGS = $(GSL_LIBS) $(ODE_LIBS) $(CONFIG_LIBS)
