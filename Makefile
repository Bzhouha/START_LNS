
FFLAG= -cpp $(SLEPC_FC_INCLUDES) $(PETSC_FC_INCLUDES)
LIBS = ${SLEPC_LIB} ${PETSC_LIB}

EXE = START_LNS
FC = 
IDIR = 
CFLAGS = -O2 -g -J$(OBJS_DIR) $(IDIR) -fallow-argument-mismatch
LFLAGS =  

OBJS_DIR = obj/
EXE_DIR = bin/

SRC_DIR_F90d1 = ./src/third_party/BeFoR64/

SRC_DIR_f90d2 = ./src/third_party/cfgio/

SRC_DIR_f90d3 = ./src/

SRC_DIR_F90d4 = ./src/third_party/PENF/

SRC_DIR_F90d5 = ./src/third_party/stringifor/

VPATH = $(SRC_DIR_F90d1):$(OBJS_DIR):$(SRC_DIR_f90d2):$(OBJS_DIR):$(SRC_DIR_f90d3):$(OBJS_DIR):$(SRC_DIR_F90d4):$(OBJS_DIR):$(SRC_DIR_F90d5):$(OBJS_DIR)
AOBJS = $(addprefix $(OBJS_DIR), $(OBJS_F90d1) $(OBJS_f90d2) $(OBJS_f90d3) $(OBJS_F90d4) $(OBJS_F90d5))

SRCS_F90d1 = \
befor64_pack_data_m.F90 \
befor64.F90 

SRCS_f90d2 = \
cfgio_mod.f90 \
string_conv_mod.f90

SRCS_f90d3 = \
mod_difference.f90 \
mod_forming.f90 \
mod_mftools.f90 \
mod_solving.f90 \
mod_loading.f90 \
cfgio_adapter.f90 \
mod_output.f90 \
mod_metrics.f90 \
mod_parameters.f90 \
main.f90 \
mod_cubes.f90 \
mod_flowtype.f90 \
petsc_viewer.f90 \
loaders.f90 \
mod_points.f90 

SRCS_F90d4 = \
penf.F90 \
penf_b_size.F90 \
penf_global_parameters_variables.F90 \
penf_stringify.F90 

SRCS_F90d5 = \
stringifor.F90 \
stringifor_string_t.F90 

OBJS_F90d1 = \
befor64_pack_data_m.o \
befor64.o 

OBJS_f90d2 = \
cfgio_mod.o \
string_conv_mod.o 

OBJS_f90d3 = \
mod_difference.o \
mod_forming.o \
mod_mftools.o \
mod_solving.o \
mod_loading.o \
cfgio_adapter.o \
mod_output.o \
mod_metrics.o \
mod_parameters.o \
main.o \
mod_cubes.o \
mod_flowtype.o \
petsc_viewer.o \
loaders.o \
mod_points.o 

OBJS_F90d4 = \
penf.o \
penf_b_size.o \
penf_global_parameters_variables.o \
penf_stringify.o 

OBJS_F90d5 = \
stringifor.o \
stringifor_string_t.o 

all : $(EXE)

$(EXE) : $(OBJS_F90d1) $(OBJS_f90d2) $(OBJS_f90d3) $(OBJS_F90d4) $(OBJS_F90d5)
	@mkdir -p $(EXE_DIR)
	$(FC) ${FFLAG} -o $(EXE_DIR)$(EXE) $(AOBJS) $(LFLAGS) $(LIBS)

$(OBJS_F90d1):
	@mkdir -p $(OBJS_DIR)
	$(FC) ${FFLAG} $(CFLAGS) -c $(SRC_DIR_F90d1)$(@:.o=.F90) -o $(OBJS_DIR)$@
$(OBJS_f90d2):
	@mkdir -p $(OBJS_DIR)
	$(FC) ${FFLAG} $(CFLAGS) -c $(SRC_DIR_f90d2)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_f90d3):
	@mkdir -p $(OBJS_DIR)
	$(FC) ${FFLAG} $(CFLAGS) -c $(SRC_DIR_f90d3)$(@:.o=.f90) -o $(OBJS_DIR)$@

$(OBJS_F90d4):
	@mkdir -p $(OBJS_DIR)
	$(FC) ${FFLAG} $(CFLAGS) -c $(SRC_DIR_F90d4)$(@:.o=.F90) -o $(OBJS_DIR)$@

$(OBJS_F90d5):
	@mkdir -p $(OBJS_DIR)
	$(FC) ${FFLAG} $(CFLAGS) -c $(SRC_DIR_F90d5)$(@:.o=.F90) -o $(OBJS_DIR)$@

clean ::
	@rm -f $(OBJS_DIR)*.o
	@rm -f $(OBJS_DIR)*.mod
	@echo " -               Clean               -"
	@echo " =      Objs have been deleted.      ="
	@echo " =   Mods have been sweeping away.   ="

install:
	@mkdir -p ~/bin
	@cp bin/START_LNS ~/bin/
	@echo " -              Install              -"
	@echo " =     hlns have been installed.     ="

befor64_pack_data_m.o: \
	befor64_pack_data_m.F90 \
	penf.o
cfgio_mod.o: \
	cfgio_mod.f90 \
	string_conv_mod.o \
	penf.o
string_conv_mod.o: \
	string_conv_mod.f90
module_gas.o: \
	module_gas.f90 \
	mod_parameters.o \
	penf.o
penf.o: \
	penf.F90 \
	penf_b_size.o \
	penf_global_parameters_variables.o \
	penf_stringify.o
penf_b_size.o: \
	penf_b_size.F90 \
	penf_global_parameters_variables.o
mod_difference.o: \
	mod_difference.f90 \
	mod_parameters.o \
	penf.o
penf_global_parameters_variables.o: \
	penf_global_parameters_variables.F90
mod_forming.o: \
	mod_forming.f90 \
	mod_flowtype.o \
	mod_parameters.o \
	mod_cubes.o \
	mod_difference.o \
	mod_mftools.o \
	penf.o 
penf_stringify.o: \
	penf_stringify.F90 \
	penf_b_size.o \
	penf_global_parameters_variables.o
mod_mftools.o: \
	mod_mftools.f90 \
	mod_parameters.o \
	penf.o
stringifor.o: \
	stringifor.F90 \
	penf.o \
	stringifor_string_t.o
mod_solving.o: \
	mod_solving.f90 \
	mod_points.o \
	mod_parameters.o \
	mod_cubes.o \
	mod_forming.o \
	mod_metrics.o 
stringifor_string_t.o: \
	stringifor_string_t.F90 \
	befor64.o \
	penf.o
mod_loading.o: \
	mod_loading.f90 \
	cfgio_adapter.o \
	loaders.o \
	petsc_viewer.o \
	mod_parameters.o
cfgio_adapter.o: \
	cfgio_adapter.f90 \
	cfgio_mod.o \
	mod_parameters.o
mod_output.o: \
	mod_output.f90 \
	mod_parameters.o \
	penf.o
mod_metrics.o: \
	mod_metrics.f90 \
	mod_parameters.o \
	mod_difference.o \
	penf.o
petsc_viewer.o: \
	petsc_viewer.f90 \
	mod_parameters.o
mod_parameters.o: \
	mod_parameters.f90 \
	mod_flowtype.o \
	penf.o
loaders.o: \
	loaders.f90 \
	mod_parameters.o
main.o: \
	main.f90 \
	mod_loading.o \
	mod_solving.o \
	mod_output.o
mod_cubes.o: \
	mod_cubes.f90 \
	mod_flowtype.o \
	mod_parameters.o \
	penf.o
mod_flowtype.o: \
	mod_flowtype.f90 \
	penf.o
mod_points.o: \
	mod_points.f90 \
	mod_parameters.o \
	mod_flowtype.o \
	penf.o
befor64.o: \
	befor64.F90 \
	befor64_pack_data_m.o \
	penf.o

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
