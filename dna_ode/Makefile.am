SUBDIRS = core sdna fdna

DNA_ODE=$(top_builddir)/dna_ode

bin_PROGRAMS = dna_ode
dna_ode_SOURCES = dna_ode.c dna_ode.h\
		  dna_ode-run.c dna_ode-run.h\
		  dna_ode-getp.c\
		  src/simcontext.c src/simcontext.h\
		  src/simparameters.c src/simparameters.h\
		  src/output_functions.c src/output_functions.h\
		  src/print_functions.c src/print_functions.h\
		  src/collision_callback.c src/collision_callback.h\
		  src/fin_functions.c src/fin_functions.h\
		  dna_ode-analyze.c dna_ode-analyze.h\
		  src/ta_info.c\
		  src/ta_angle_hist.c\
		  src/ta_distance_hist.c\
		  src/ta_distance_map.c\
		  src/ta_body_trace.c\
		  src/ta_bead_force.c\
		  src/ta_last_segment_position.c\
		  src/ta_frame_dump.c


dna_ode_CPPFLAGS = $(GSL_CFLAGS) $(ODE_CFLAGS) -I$(top_srcdir)/dna_ode/src -I$(top_srcdir)/dna_ode/core
dna_ode_LDADD = $(top_builddir)/dna_ode/core/libdna_ode_core.la -L$(top_builddir)/dna_ode/core\
		$(DNA_ODE)/sdna/libsdna.la

AM_LDFLAGS = $(GSL_LIBS) $(ODE_LIBS) $(CONFIG_LIBS)
AM_CFLAGS = -I$(top_srcdir)/dna_ode/src -Isdna
