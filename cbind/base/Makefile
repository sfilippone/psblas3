TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)lib
INCLUDEDIR=$(TOP)/include
HERE=..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(INCLUDEDIR) 
CINCLUDES=-I$(INCLUDEDIR)

FOBJS= psb_objhandle_mod.o psb_base_cbind_mod.o psb_cpenv_mod.o \
      psb_base_tools_cbind_mod.o psb_base_string_cbind_mod.o \
      psb_base_psblas_cbind_mod.o psb_d_tools_cbind_mod.o  \
      psb_d_serial_cbind_mod.o psb_d_psblas_cbind_mod.o

COBJS= psb_c_base.o  psb_c_dbase.o 
CMOD=psb_base_cbind.h psb_c_base.h psb_c_dbase.h 
OBJS=$(FOBJS) $(COBJS)

LIBMOD=psb_base_cbind_mod$(.mod) psb_cpenv_mod$(.mod)  psb_objhandle_mod$(.mod)\
       psb_base_tools_cbind_mod$(.mod) psb_base_string_cbind_mod$(.mod) psb_base_psblas_cbind_mod$(.mod)\
       psb_d_tools_cbind_mod$(.mod) psb_d_serial_cbind_mod$(.mod) psb_d_psblas_cbind_mod$(.mod)
LOCAL_MODS=$(LIBMOD) 
LIBNAME=$(CBINDLIBNAME)	


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(LIBMOD) $(CMOD) $(HERE)

$(COBJS): $(CMOD)

psb_base_cbind_mod.o: psb_cpenv_mod.o psb_objhandle_mod.o psb_base_tools_cbind_mod.o \
	 psb_base_string_cbind_mod.o psb_base_psblas_cbind_mod.o
psb_base_tools_cbind_mod.o: psb_objhandle_mod.o  psb_base_string_cbind_mod.o
psb_d_tools_cbind_mod.o psb_d_serial_cbind_mod.o: psb_base_tools_cbind_mod.o

psb_base_psblas_cbind_mod.o: psb_d_psblas_cbind_mod.o
psb_d_psblas_cbind_mod.o: psb_objhandle_mod.o  psb_base_string_cbind_mod.o 
psb_cpenv_mod.o: psb_base_string_cbind_mod.o psb_objhandle_mod.o 

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
