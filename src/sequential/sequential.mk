# Authors: Jonathan, Adrián, Adrià

# Specific compiler
FC := gfortran

# Directory for this module
SEQ_DIR := src/sequential
COM_DIR := src/common

# Build dir
BUILD_DIR := build

# Output dirs for sequential
SEQ_OBJ_DIR := $(BUILD_DIR)/obj/sequential
MOD_OUT_DIR := $(BUILD_DIR)/modules

# Compiler flags specific to sequential
SEQ_FLAGS := -O3 \
             -I$(SEQ_DIR)/initialization \
             -I$(SEQ_DIR)/energy \
             -I$(SEQ_DIR)/MC_update \
             -I$(COM_DIR)/initialization \
             -I$(COM_DIR)/energy \
             -I$(COM_DIR)/MC_update \
             -J$(MOD_OUT_DIR)

# Executable name
SEQ_EXEC := $(BUILD_DIR)/mainglobal_seq.exe

# Source files
SEQ_SRC := $(COM_DIR)/initialization/constants.f90 \
           $(SEQ_DIR)/initialization/system.f90 \
           $(SEQ_DIR)/initialization/io_module.f90 \
           $(COM_DIR)/initialization/init_config.f90 \
           $(COM_DIR)/initialization/centerPolymer.f90 \
           $(SEQ_DIR)/energy/nonbonded.f90 \
           $(COM_DIR)/energy/bonded.f90 \
           $(COM_DIR)/energy/energy.f90 \
           $(COM_DIR)/MC_update/MC_proposal.f90 \
           $(SEQ_DIR)/MC_update/MC_loop.f90 \
           $(SEQ_DIR)/mainglobal.f90

# Object files (in build directory)
SEQ_OBJ := $(patsubst $(SEQ_DIR)/%.f90, $(SEQ_OBJ_DIR)/%.o, \
           $(patsubst $(COM_DIR)/%.f90, $(SEQ_OBJ_DIR)/%.o, $(SEQ_SRC)))

# Rules specific to the sequential module
$(SEQ_EXEC): $(SEQ_OBJ)
	@mkdir -p $(dir $@)
	$(FC) -o $@ $(SEQ_OBJ)

# Pattern rule for compiling Fortran objects in this module
$(SEQ_OBJ_DIR)/%.o: $(SEQ_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_OUT_DIR)
	$(FC) $(SEQ_FLAGS) -c $< -o $@

$(SEQ_OBJ_DIR)/%.o: $(COM_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_OUT_DIR)
	$(FC) $(SEQ_FLAGS) -c $< -o $@

# Local clean target
clean-seq:
	rm -rf $(SEQ_OBJ_DIR) $(SEQ_EXEC) $(MOD_OUT_DIR)/*.mod