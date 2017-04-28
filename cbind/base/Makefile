TOP=../..
include $(TOP)/Make.inc
LIBDIR=$(TOP)lib
INCLUDEDIR=$(TOP)/include
HERE=..

FINCLUDES=$(FMFLAG). $(FMFLAG)$(HERE) $(FMFLAG)$(INCLUDEDIR) 
CINCLUDES=-I. -I$(HERE) -I$(INCLUDEDIR)

FOBJS= psb_objhandle_mod.o psb_base_cbind_mod.o psb_cpenv_mod.o \
      psb_base_tools_cbind_mod.o psb_base_string_cbind_mod.o \
      psb_base_psblas_cbind_mod.o \
      psb_s_tools_cbind_mod.o  psb_s_serial_cbind_mod.o psb_s_psblas_cbind_mod.o \
      psb_d_tools_cbind_mod.o  psb_d_serial_cbind_mod.o psb_d_psblas_cbind_mod.o \
      psb_c_tools_cbind_mod.o  psb_c_serial_cbind_mod.o psb_c_psblas_cbind_mod.o \
      psb_z_tools_cbind_mod.o  psb_z_serial_cbind_mod.o psb_z_psblas_cbind_mod.o \
      psb_s_comm_cbind_mod.o psb_d_comm_cbind_mod.o \
      psb_c_comm_cbind_mod.o psb_z_comm_cbind_mod.o 

COBJS= psb_c_base.o psb_c_sbase.o psb_c_dbase.o psb_c_cbase.o psb_c_zbase.o \
	psb_c_dcomm.o
CMOD=psb_base_cbind.h psb_c_base.h psb_c_sbase.h  psb_c_dbase.h  psb_c_cbase.h  psb_c_zbase.h \
	psb_c_dcomm.h
OBJS=$(FOBJS) $(COBJS)

LIBMOD=psb_base_cbind_mod$(.mod) psb_cpenv_mod$(.mod)  psb_objhandle_mod$(.mod)\
       psb_base_tools_cbind_mod$(.mod) psb_base_string_cbind_mod$(.mod) psb_base_psblas_cbind_mod$(.mod)\
       psb_s_tools_cbind_mod$(.mod) psb_s_serial_cbind_mod$(.mod) psb_s_psblas_cbind_mod$(.mod) \
       psb_d_tools_cbind_mod$(.mod) psb_d_serial_cbind_mod$(.mod) psb_d_psblas_cbind_mod$(.mod) \
       psb_c_tools_cbind_mod$(.mod) psb_c_serial_cbind_mod$(.mod) psb_c_psblas_cbind_mod$(.mod) \
       psb_z_tools_cbind_mod$(.mod) psb_z_serial_cbind_mod$(.mod) psb_z_psblas_cbind_mod$(.mod) \
       psb_s_comm_cbind_mod$(.mod) psb_d_comm_cbind_mod$(.mod) \
       psb_c_comm_cbind_mod$(.mod) psb_z_comm_cbind_mod$(.mod) 

LOCAL_MODS=$(LIBMOD) 
LIBNAME=$(CBINDLIBNAME)	


lib:  $(OBJS) $(CMOD)
	$(AR) $(HERE)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(HERE)/$(LIBNAME)
	/bin/cp -p $(LIBMOD) $(CMOD) $(HERE)

$(COBJS): $(CMOD)

psb_base_cbind_mod.o: psb_cpenv_mod.o psb_objhandle_mod.o psb_base_tools_cbind_mod.o \
	 psb_base_string_cbind_mod.o psb_base_psblas_cbind_mod.o \
      psb_s_tools_cbind_mod.o  psb_s_serial_cbind_mod.o psb_s_psblas_cbind_mod.o \
      psb_d_tools_cbind_mod.o  psb_d_serial_cbind_mod.o psb_d_psblas_cbind_mod.o \
      psb_c_tools_cbind_mod.o  psb_c_serial_cbind_mod.o psb_c_psblas_cbind_mod.o \
      psb_z_tools_cbind_mod.o  psb_z_serial_cbind_mod.o psb_z_psblas_cbind_mod.o  \
      psb_s_comm_cbind_mod.o psb_d_comm_cbind_mod.o \
      psb_c_comm_cbind_mod.o psb_z_comm_cbind_mod.o 

psb_base_tools_cbind_mod.o: psb_objhandle_mod.o  psb_base_string_cbind_mod.o

psb_s_tools_cbind_mod.o psb_s_serial_cbind_mod.o \
psb_d_tools_cbind_mod.o psb_d_serial_cbind_mod.o \
psb_c_tools_cbind_mod.o psb_c_serial_cbind_mod.o \
psb_z_tools_cbind_mod.o psb_z_serial_cbind_mod.o \
psb_s_psblas_cbind_mod.o psb_d_psblas_cbind_mod.o \
psb_c_psblas_cbind_mod.o psb_z_psblas_cbind_mod.o \
psb_d_comm_cbind_mod.o : psb_base_tools_cbind_mod.o psb_objhandle_mod.o  psb_base_string_cbind_mod.o

psb_base_psblas_cbind_mod.o: psb_s_psblas_cbind_mod.o psb_d_psblas_cbind_mod.o psb_c_psblas_cbind_mod.o psb_z_psblas_cbind_mod.o

psb_cpenv_mod.o: psb_base_string_cbind_mod.o psb_objhandle_mod.o 

veryclean: clean
	/bin/rm -f $(HERE)/$(LIBNAME)

clean:
	/bin/rm -f $(OBJS) $(LOCAL_MODS)

veryclean: clean
