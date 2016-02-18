conf_files = fdna.conf sdna.conf

conf_file_sources = dna_parameters.conf \
		    environment_parameters.conf \
		    ode_parameters.conf \
		    mechanical_models.conf \
		    thermostat_parameters.conf \
		    forces_parameters.conf \
		    external_forces_parameters.conf \
		    bead_parameters.conf \
		    surface_parameters.conf \
		    time_parameters.conf \
		    file_names.conf \
		    fdna.conf.in \
		    sdna.conf.in

SUFFIXES = .conf.in .conf

$(conf_files) : $(conf_file_sources)

.conf.in.conf : $(conf_file_sources)
	@echo "Building $@"
	@cat $< |\
	sed -e '/@ODE_PARAMETERS@/{r ode_parameters.conf' -e 'd}' |\
	sed -e '/@ENVIRONMENT_PARAMETERS@/{r environment_parameters.conf' -e 'd}' |\
	sed -e '/@BEAD_PARAMETERS@/{r bead_parameters.conf' -e 'd}' |\
	sed -e '/@SURFACE_PARAMETERS@/{r surface_parameters.conf' -e 'd}' |\
	sed -e '/@DNA_PARAMETERS@/{r dna_parameters.conf' -e 'd}' |\
	sed -e '/@MECHANICAL_MODEL@/{r mechanical_models.conf' -e 'd}' |\
	sed -e '/@THERMOSTAT_PARAMETERS@/{r thermostat_parameters.conf' -e 'd}' |\
	sed -e '/@FORCES_PARAMETERS@/{r forces_parameters.conf' -e 'd}' |\
	sed -e '/@EXTERNAL_FORCES_PARAMETERS@/{r external_forces_parameters.conf' -e 'd}' |\
	sed -e '/@FILE_NAMES@/{r file_names.conf' -e 'd}' |\
	sed -e '/@TIME_PARAMETERS@/{r time_parameters.conf' -e 'd}' |\
	sed -e s,@SYSTEM_ID@,$*,g |\
	tee > $@
	@cp $@ ..

all : $(conf_files) all-am

CLEANFILES = $(conf_files)

EXTRA_DIST = $(conf_file_sources)
