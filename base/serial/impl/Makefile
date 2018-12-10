include ../../../Make.inc

#
# The object files
#
BOBJS=psb_base_mat_impl.o \
  psb_s_base_mat_impl.o psb_d_base_mat_impl.o psb_c_base_mat_impl.o psb_z_base_mat_impl.o
#\
  psb_s_lbase_mat_impl.o psb_d_lbase_mat_impl.o psb_c_lbase_mat_impl.o psb_z_lbase_mat_impl.o 
SOBJS=psb_s_csr_impl.o psb_s_coo_impl.o psb_s_csc_impl.o psb_s_mat_impl.o
#\
 psb_s_lcoo_impl.o psb_s_lcsr_impl.o 
DOBJS=psb_d_csr_impl.o psb_d_coo_impl.o psb_d_csc_impl.o psb_d_mat_impl.o
#\
 psb_d_lcoo_impl.o  psb_d_lcsr_impl.o 
COBJS=psb_c_csr_impl.o psb_c_coo_impl.o psb_c_csc_impl.o psb_c_mat_impl.o
#\
 psb_c_lcoo_impl.o  psb_c_lcsr_impl.o 
ZOBJS=psb_z_csr_impl.o psb_z_coo_impl.o psb_z_csc_impl.o psb_z_mat_impl.o
#\
 psb_z_lcoo_impl.o  psb_z_lcsr_impl.o 

OBJS=$(BOBJS) $(SOBJS)  $(DOBJS)  $(COBJS)  $(ZOBJS) 

#
# Where the library should go, and how it is called. 
# Note that we are regenerating most of libsparker.a on the fly. 
LIBDIR=../..
INCDIR=../..
MODDIR=../../modules
FINCLUDES=$(FMFLAG). $(FMFLAG)$(MODDIR) $(FMFLAG)$(INCDIR) 
LIBFILE=$(LIBDIR)/$(LIBNAME)

#
# No change should be needed below 
#


default: lib

lib: $(OBJS)
	$(AR) $(LIBDIR)/$(LIBNAME) $(OBJS)
	$(RANLIB) $(LIBDIR)/$(LIBNAME)

# A bit excessive, but safe
$(OBJS): $(MODDIR)/psb_base_mod.o

clean: cleanobjs

veryclean: cleanobjs

cleanobjs:
	/bin/rm -f $(OBJS) 

