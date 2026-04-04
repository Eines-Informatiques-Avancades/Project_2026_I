# =========================
# Makefile MPI – Parallel Build
# =========================

# Compiler
FC := mpifort

# Project roots
PAR_DIR := src/parallel
COM_DIR := src/common

# Build directories
BUILD_DIR := build
PAR_OBJ_DIR := $(BUILD_DIR)/obj/parallel
MOD_OUT_DIR := $(BUILD_DIR)/modules

# Compiler flags
PAR_FLAGS := -O3 \
             -I$(PAR_DIR)/initialization \
             -I$(PAR_DIR)/energy \
             -I$(PAR_DIR)/MC_update \
             -I$(COM_DIR)/initialization \
             -I$(COM_DIR)/energy \
             -I$(COM_DIR)/MC_update \
             -J$(MOD_OUT_DIR)

# Executable
PAR_EXEC := $(BUILD_DIR)/mainglobal_mpi.exe

# Source files
PAR_SRC := $(COM_DIR)/initialization/constants.f90 \
           $(PAR_DIR)/initialization/system.f90 \
           $(PAR_DIR)/initialization/io_module.f90 \
           $(COM_DIR)/initialization/init_config.f90 \
           $(COM_DIR)/initialization/centerPolymer.f90 \
           $(PAR_DIR)/energy/nonbonded.f90 \
           $(COM_DIR)/energy/bonded.f90 \
           $(COM_DIR)/energy/energy.f90 \
           $(COM_DIR)/MC_update/MC_proposal.f90 \
           $(PAR_DIR)/MC_update/MC_loop.f90 \
           $(PAR_DIR)/mainglobal.f90

# Object files (mirroring folder structure taking PAR_DIR and COM_DIR sources into PAR_OBJ_DIR)
PAR_OBJ := $(patsubst $(PAR_DIR)/%.f90, $(PAR_OBJ_DIR)/%.o, \
           $(patsubst $(COM_DIR)/%.f90, $(PAR_OBJ_DIR)/%.o, $(PAR_SRC)))

# =========================
# Build executable
# =========================
$(PAR_EXEC): $(PAR_OBJ)
	@mkdir -p $(dir $@)
	$(FC) -o $@ $(PAR_OBJ)

# =========================
# Pattern rule for object files
# =========================
$(PAR_OBJ_DIR)/%.o: $(PAR_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_OUT_DIR)
	$(FC) $(PAR_FLAGS) -c $< -o $@

$(PAR_OBJ_DIR)/%.o: $(COM_DIR)/%.f90
	@mkdir -p $(dir $@) $(MOD_OUT_DIR)
	$(FC) $(PAR_FLAGS) -c $< -o $@

# =========================
# Clean build files
# =========================
clean-par:
	rm -rf $(PAR_OBJ_DIR) $(PAR_EXEC) $(MOD_OUT_DIR)/*.mod