include Makefile.vars

progs_cpp:=$(filter %.cpp,$(shell ls src/progs))

progs:= $(patsubst %.cpp, run_%, $(progs_cpp))

.PHONY: all clean clean_gcov print
.PRECIOUS: %.o

ifeq ($(BLAS_LIB), mkl)
BLAS_LIB:=-lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -lm -ldl -liomp5 
BLAS_DEF:=-D_MKL
BLAS_LDIR:=-L$(BLAS_DIR)/lib/intel64 $(BLAS_EXTRA_LIB)
endif
ifeq ($(BLAS_LIB), openblas)
BLAS_LIB:=-lopenblas
BLAS_DEF:=-D_OPENBLAS
BLAS_LDIR:=-L$(BLAS_DIR)/lib
endif
ifeq ($(BLAS_PAR), blas)
BLAS_PAR:=-D_SELF_BLAS_PAR
endif
ifeq ($(BLAS_PAR), nopt)
BLAS_PAR:=-D_NOPT_BLAS_PAR
endif

ifeq ($(USE_GRPP), yes)
GRPP_DEF:=-D_USE_GRPP
GRPP_LIB:=-llibgrpp
GRPP_H:=-I$(GRPP_H_DIR)
GRPP_LDIR:=-L$(GRPP_L_DIR)
endif
ifeq ($(USE_GRPP), no)
GRPP_DEF:=
GRPP_H_DIR:=
GRPP_L_DIR:=
GRPP_LIB:=
endif

ifeq ($(USE_XDR),yes)
XDR_DEF:=-D_XDR_FULL
XDR_LIB:=-ltirpc
XDR_H:=-I$(XDR_DIR)
endif

ifeq ($(USE_XDR),read_only)
XDR_DEF:=-D_XDR_READ
XDR_LIB:=-ltirpc
XDR_H:=-I$(XDR_DIR)
endif



LIBS:=-lgfortran -lint2 -lpthread $(XDR_LIB) $(BLAS_LIB) $(GRPP_LIB)

DEFINITIONS:=-DADD_ -DLIBINT -DLBFGS_RETURN_STEP -DLBFGS_CALC_HESS -D_RI -D_LIBINT_INITIAL_SHELLS -DLIBINT2_DISABLE_BOOST_CONTAINER_SMALL_VECTOR $(BLAS_DEF) $(BLAS_PAR) -DNOPT_LIB="\"$(NOPT_lib_DIR)\"" -DL_MAX=$(L_MAX) -DRI_L_MAX=$(RI_L_MAX) -DGTO_MAX=$(GTO_MAX) $(extra_def) $(GRPP_DEF) $(XDR_DEF)

INCLUDE_DIRS:=-Iinclude -I$(BLAS_DIR)/include -I$(LIBINT_H_DIR) -I$(EIGEN_DIR) $(XDR_H) $(GRPP_H)
ifeq ($(MANUAL_BOOST_PATH),yes)
INCLUDE_DIRS:= $(INCLUDE_DIRS) -I$(BOOST_DIR)
endif


LIB_DIRS:=$(BLAS_LDIR) -L$(LIBINT_L_DIR) $(GRPP_LDIR)



all:$(progs)

print:
	@echo $(progs_cpp)
	@echo $(progs)

run_%:src/progs/%.o src/molecule.o src/molecule2.o src/chem_data.o src/MO2.o src/matr.o src/timer.o src/etc.o src/doCI_matr.o src/libint_link.o src/SCF.o src/converger_2_1.o src/l-bfgs_2_1.o src/libint_functions.o src/ecp.o src/mol_link.o src/CI.o src/inp_out.o src/keywords.o src/doCI_data.o src/inp_par_read.o src/from_hash.o src/XMCQDPT.o src/res_fit.o src/binary_files.o src/xmc_read.o src/aldet.o src/RI.o src/nopa_pt.o src/CAS.o src/davidson.o src/basis_lib_read.o src/common_vars.o src/CDAS_PT.o src/CDAS_PT_rel.o src/U_CDAS_PT.o src/PT_tensors_IPEA.o src/PT_tensors_EE.o src/PT_tensors_EE_rel.o src/UPT_tensors_EE.o src/geom.o src/symmetry.o src/gv_solver.o src/trcamm.o src/z_matrix.o src/QM_calc.o src/pseudo_potential.o src/complex_diag.o src/CIS.o src/aldet_rel.o src/jacobi.o src/grabbers.o src/casci_solver.o src/aldet_casci_wrap.o

	$(CXX) $^ $(OPT_LEVEL) -fopenmp $(LIB_DIRS) $(LIBS)  -o $@

%.o:%.cpp
	$(CXX) -o $@ -c $< $(OPT_LEVEL) -fopenmp $(DEFINITIONS) $(INCLUDE_DIRS) -fmax-errors=5

%.o:%.c
	$(CC) -o $@ -c $< $(OPT_LEVEL) -fopenmp $(DEFINITIONS) $(INCLUDE_DIRS) -fmax-errors=5

clean:
	rm $(progs) src/*.o src/progs/*.o
clean_gcov:
	rm src/*.gcda src/progs/*.gcda src/*.gcno src/progs/*.gcno
