
FFLAG= -cpp -fallow-argument-mismatch $(SLEPC_FC_INCLUDES) $(PETSC_FC_INCLUDES)
LIBS = ${SLEPC_LIB} ${PETSC_LIB}

EXE = START_LNS
FC = 
IDIR = 
CFLAGS = -g -I$(OBJS_DIR) $(IDIR)
LFLAGS =  

VPATH = $(SRC_DIR_F90d1):$(OBJS_DIR):$(SRC_DIR_f90d2):$(OBJS_DIR):$(SRC_DIR_f90d3):$(OBJS_DIR):$(SRC_DIR_F90d4):$(OBJS_DIR):$(SRC_DIR_F90d5):$(OBJS_DIR)
AOBJS = $(addprefix $(OBJS_DIR), $(OBJS_F90d1) $(OBJS_f90d2) $(OBJS_f90d3) $(OBJS_F90d4) $(OBJS_F90d5))

SRCS_F90d1 = \
befor64_pack_data_m.F90 \
befor64.F90 

SRCS_f90d2 = \
cfgio_mod.f90 \
string_conv_mod.f90

SRCS_f90d3 = \
module_bf_point_org.f90 \
mod_difference.f90 \
mod_forming.f90 \
mod_mf_tools.f90 \
mod_solving.f90 \
petsc_loader.f90 \
cfgio_adapter.f90 \
petsc_output.f90 \
compute_metrics.f90 \
global_parameters.f90 \
main.f90 \
matrix_used_as_cofficient.f90 \
module_baseflow_org.f90 \
petsc_viewer.f90 \
loaders.f90 \
module_bf_points.f90 

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
module_bf_point_org.o \
mod_difference.o \
mod_forming.o \
mod_mf_tools.o \
mod_solving.o \
petsc_loader.o \
cfgio_adapter.o \
petsc_output.o \
compute_metrics.o \
global_parameters.o \
main.o \
matrix_used_as_cofficient.o \
module_baseflow_org.o \
petsc_viewer.o \
loaders.o \
module_bf_points.o 

OBJS_F90d4 = \
penf.o \
penf_b_size.o \
penf_global_parameters_variables.o \
penf_stringify.o 

OBJS_F90d5 = \
stringifor.o \
stringifor_string_t.o 


SRC_DIR_F90d1 = ./src/third_party/BeFoR64/

SRC_DIR_f90d2 = ./src/third_party/cfgio/

SRC_DIR_f90d3 = ./src/

SRC_DIR_F90d4 = ./src/third_party/PENF/

SRC_DIR_F90d5 = ./src/third_party/stringifor/


OBJS_DIR = ./obj/
EXE_DIR = bin/

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
	@rm -f $(OBJS_DIR)*.*
	@rm -f *.mod 
	@rm -f $(EXE_DIR)$(EXE)
	@echo " -               Clean               -"
	@echo " =      Objs have been deleted.      ="
	@echo " =   Mods have been sweeping away.   ="

install:
	mkdir -p ~/bin
	cp bin/START_LNS ~/bin/

befor64_pack_data_m.o: \
	befor64_pack_data_m.F90 \
	penf.o
cfgio_mod.o: \
	cfgio_mod.f90 \
	string_conv_mod.o \
	penf.o
module_bf_point_org.o: \
	module_bf_point_org.f90 \
	module_baseflow_org.o
string_conv_mod.o: \
	string_conv_mod.f90
module_gas.o: \
	module_gas.f90 \
	global_parameters.o \
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
	global_parameters.o \
	penf.o
penf_global_parameters_variables.o: \
	penf_global_parameters_variables.F90
mod_forming.o: \
	mod_forming.f90 \
	module_bf_point_org.o \
	global_parameters.o \
	matrix_used_as_cofficient.o \
	mod_difference.o \
	mod_mf_tools.o \
	penf.o 
penf_stringify.o: \
	penf_stringify.F90 \
	penf_b_size.o \
	penf_global_parameters_variables.o
mod_mf_tools.o: \
	mod_mf_tools.f90 \
	global_parameters.o \
	penf.o
stringifor.o: \
	stringifor.F90 \
	penf.o \
	stringifor_string_t.o
mod_solving.o: \
	mod_solving.f90 \
	module_bf_points.o \
	global_parameters.o \
	matrix_used_as_cofficient.o \
	mod_forming.o \
	compute_metrics.o 
stringifor_string_t.o: \
	stringifor_string_t.F90 \
	befor64.o \
	penf.o
petsc_loader.o: \
	petsc_loader.f90 \
	cfgio_adapter.o \
	loaders.o \
	petsc_viewer.o \
	global_parameters.o
cfgio_adapter.o: \
	cfgio_adapter.f90 \
	cfgio_mod.o \
	global_parameters.o
petsc_output.o: \
	petsc_output.f90 \
	global_parameters.o \
	penf.o
compute_metrics.o: \
	compute_metrics.f90 \
	global_parameters.o \
	mod_difference.o \
	penf.o
petsc_viewer.o: \
	petsc_viewer.f90 \
	global_parameters.o
global_parameters.o: \
	global_parameters.f90 \
	module_bf_point_org.o \
	penf.o
loaders.o: \
	loaders.f90 \
	global_parameters.o
main.o: \
	main.f90 \
	petsc_loader.o \
	mod_solving.o \
	petsc_output.o
matrix_used_as_cofficient.o: \
	matrix_used_as_cofficient.f90 \
	module_bf_point_org.o \
	global_parameters.o \
	penf.o
module_baseflow_org.o: \
	module_baseflow_org.f90 \
	penf.o
module_bf_points.o: \
	module_bf_points.f90 \
	module_bf_point_org.o \
	global_parameters.o \
	module_baseflow_org.o \
	penf.o
befor64.o: \
	befor64.F90 \
	befor64_pack_data_m.o \
	penf.o

include ${SLEPC_DIR}/lib/slepc/conf/slepc_common
