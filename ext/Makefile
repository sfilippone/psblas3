include ../Make.inc
#
# Libraries used
#
LIBDIR=../lib
INCDIR=../include
MODDIR=../modules
#
# Compilers and such
#
#CCOPT= -g
FINCLUDES=$(FMFLAG). $(FMFLAG)$(INCDIR) $(FMFLAG)$(MODDIR) $(FIFLAG). 
CINCLUDES=
LIBNAME=libpsb_ext.a


FOBJS=  psb_d_ell_mat_mod.o psb_d_hll_mat_mod.o \
	psb_s_hll_mat_mod.o psb_s_ell_mat_mod.o \
	psb_c_hll_mat_mod.o psb_c_ell_mat_mod.o \
	psb_z_hll_mat_mod.o psb_z_ell_mat_mod.o \
	psb_d_dia_mat_mod.o psb_d_hdia_mat_mod.o \
	psb_s_dia_mat_mod.o psb_s_hdia_mat_mod.o \
	psb_c_dia_mat_mod.o psb_c_hdia_mat_mod.o \
	psb_z_dia_mat_mod.o psb_z_hdia_mat_mod.o \
	psb_s_dns_mat_mod.o psb_d_dns_mat_mod.o \
	psb_c_dns_mat_mod.o psb_z_dns_mat_mod.o \
	psi_ext_util_mod.o psi_i_ext_util_mod.o \
	psi_s_ext_util_mod.o psi_c_ext_util_mod.o \
	psi_d_ext_util_mod.o psi_z_ext_util_mod.o \
	psb_ext_mod.o

COBJS= 

OBJS=$(COBJS) $(FOBJS)

lib: objs ilib 
	ar cur $(LIBNAME) $(OBJS)
	/bin/cp -p $(LIBNAME) $(LIBDIR)

objs: $(OBJS) iobjs
	/bin/cp -p *$(.mod) $(MODDIR)



psb_ext_mod.o: psb_s_dia_mat_mod.o   psb_d_dia_mat_mod.o  \
	psb_c_dia_mat_mod.o   psb_z_dia_mat_mod.o  \
        psb_d_ell_mat_mod.o psb_d_hll_mat_mod.o  \
	psb_s_hll_mat_mod.o psb_s_ell_mat_mod.o \
	psb_c_hll_mat_mod.o psb_c_ell_mat_mod.o \
	psb_z_hll_mat_mod.o psb_z_ell_mat_mod.o \
	psb_s_hdia_mat_mod.o psb_d_hdia_mat_mod.o \
	psb_c_hdia_mat_mod.o psb_z_hdia_mat_mod.o \
	psb_s_dns_mat_mod.o psb_d_dns_mat_mod.o \
	psb_c_dns_mat_mod.o psb_z_dns_mat_mod.o 

#	psb_d_rsb_mat_mod.o psb_d_hdia_mat_mod.o 
psi_ext_util_mod.o: psi_i_ext_util_mod.o \
	psi_s_ext_util_mod.o psi_c_ext_util_mod.o \
	psi_d_ext_util_mod.o psi_z_ext_util_mod.o 

psb_s_dia_mat_mod.o psb_c_dia_mat_mod.o psb_d_dia_mat_mod.o psb_z_dia_mat_mod.o: psi_ext_util_mod.o
psb_s_hdia_mat_mod.o psb_c_hdia_mat_mod.o psb_d_hdia_mat_mod.o psb_z_hdia_mat_mod.o: psi_ext_util_mod.o
psb_s_hll_mat_mod.o psb_c_hll_mat_mod.o psb_d_hll_mat_mod.o psb_z_hll_mat_mod.o: psi_ext_util_mod.o

ilib: objs
	$(MAKE) -C impl lib LIBNAME=$(LIBNAME)

iobjs: $(OBJS)
	$(MAKE) -C impl objs 

clean: cclean iclean
	/bin/rm -f  $(FOBJS) *$(.mod) *.a

cclean: 
	/bin/rm -f  $(COBJS) 
iclean:
	$(MAKE) -C impl clean

veryclean: clean 
	/bin/rm -f $(HERE)/$(LIBNAME) $(LIBMOD) *$(.mod)




