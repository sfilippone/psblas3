/* This header file replaces every call to a BLACS routine by C interface
   with the same call performed by Fortran interface */

#ifndef CTOF_BLACS
#define CTOF_BLACS
#endif

/* Variables necessary for invocations where 
   constant arguments are used */

static int i1, i2, i3, i4, i5, i6, i7;

/* Support routines:
   Initialization */

#define Cblacs_pinfo(mypnum, nprocs) \
        blacs_pinfo_(mypnum, nprocs)
#define Cblacs_setup(mypnum, nprocs) \
        blacs_setup_(mypnum, nprocs)
#define Cblacs_get(icontxt, what, val) \
        {i1 = icontxt; i2 = what; \
        blacs_get_(&i1, &i2,val);}
#define Cblacs_set(icontxt, what, val) \
        {i1 = icontxt; i2 = what; \
        blacs_set_(&i1, &i2, &val);}
#define Cblacs_gridinit(icontxt, order, nprow, npcol) \
        {i1 = nprow; i2 = npcol; \
        blacs_gridinit_(icontxt, order, &i1, &i2);}
#define Cblacs_gridmap(icontxt, pmap, ldpmap, nprow, npcol) \
        {i1 = ldpmap; i2 = nprow; i3 = npcol; \
        blacs_gridmap_(icontxt, pmap, &i1, &i2, &i3);}

/* Support routines:
   Destruction */

#define Cblacs_freebuff(icontxt, wait) \
        {i1 = icontxt; i2 = wait; \
        blacs_freebuff_(&i1, &i2);}
#define Cblacs_gridexit(icontxt) \
        {i1 = icontxt; \
        blacs_gridexit_(&i1);}
#define Cblacs_abort(icontxt, errornum) \
        {i1 = icontxt; i2 = errornum; \
        blacs_abort_(&i1, &i2);}
#define Cblacs_exit(doneflag) \
        {i1 = doneflag; \
        blacs_exit_(&i1);}

/* Support routines:
   Informational and Miscellaneous */

#define Cblacs_gridinfo(icontxt,nprow,npcol,myprow,mypcol) \
        {i1 = icontxt; \
        blacs_gridinfo_(&i1, nprow, npcol, myprow, mypcol);}
#define Cblacs_pnum(icontxt, prow, pcol) \
        {i1 = icontxt; i2 = prow; i3 = pcol; \
        blacs_pnum_(&i1, &i2, &i3);}
#define Cblacs_pcoord(icontxt, pnum, prow, pcol) \
        {i1 = icontxt; i2 = pnum; \
        blacs_pcoord_(&i1, &i2, prow, pcol);}
#define Cblacs_barrier(icontxt, scope) \
        {i1 = icontxt; \
        blacs_barrier_(&i1, scope);}

/* Support routines:
   Unofficial */ 

#define Csetpvmtids(ntasks, tids) \
        {i1 = ntasks; \
        setpvmtids_(&i1, tids);}
#define Cdcputime() \
        dcputime_()
#define Cdwalltime() \
        dwalltime_()
#define Cksendid(icontxt, rdest, cdest) \
        {i1 = icontxt; i2 = rdest; i3 = cdest; \
        ksendid_(&i1, &i2, &i3);}
#define Ckrecvid(icontxt, rsrc, csrc) \
        {i1 = icontxt; i2 = rsrc; i3 = csrc; \
        krecvid_(&i1, &i2, &i3);}
#define Ckbsid(icontxt, scope) \
        {i1 = icontxt; \
        kbsid_(&i1, scope);}
#define Ckbrid(icontxt, scope, rsrc, csrc) \
        {i1 = icontxt; i2 = rsrc; i3 = csrc; \
        kbrid_(&i1, scope, &i2, &i3);}

/* Point to Point :
   Integer */

#define Cigesd2d(icontxt, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        igesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cigerv2d(icontxt, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        igerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Citrsd2d(icontxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        itrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Citrrv2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        itrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Point to Point :
   Single precision real */

#define Csgesd2d(icontxt, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        sgesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Csgerv2d(icontxt, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        sgerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cstrsd2d(icontxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        strsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Cstrrv2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        strsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Point to Point :
   Double precision real */

#define Cdgesd2d(icontxt, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        dgesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdgerv2d(icontxt, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        dgerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdtrsd2d(icontxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        dtrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdtrrv2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        dtrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Point to Point :
   Single precision complex */

#define Ccgesd2d(icontxt, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        cgesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Ccgerv2d(icontxt, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        cgerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cctrsd2d(icontxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        ctrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Cctrrv2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        ctrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Point to Point :
   Double precision complex */

#define Czgesd2d(icontxt, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        zgesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Czgerv2d(icontxt, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        zgerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cztrsd2d(icontxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        ztrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Cztrrv2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        ztrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Integer */

#define Cigebs2d(icontxt, scope, top, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        igebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Cigebr2d(icontxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        igebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Citrbs2d(icontxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        itrbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Citrbr2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        igebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Single precision real */

#define Csgebs2d(icontxt, scope, top, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        sgebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Csgebr2d(icontxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        sgebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cstrbs2d(icontxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        strbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Cstrbr2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        sgebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Double precision real */

#define Cdgebs2d(icontxt, scope, top, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        dgebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Cdgebr2d(icontxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        dgebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdtrbs2d(icontxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        dtrbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Cdtrbr2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        dgebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Single precision complex */

#define Ccgebs2d(icontxt, scope, top, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        cgebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Ccgebr2d(icontxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        cgebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cctrbs2d(icontxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        ctrbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Cctrbr2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        cgebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Double precision complex */

#define Czgebs2d(icontxt, scope, top, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        zgebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Czgebr2d(icontxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        zgebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cztrbs2d(icontxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; \
        ztrbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Cztrbr2d(icontxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        zgebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Combines:
   Integer */

#define Cigsum2d(icontxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        igsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cigamx2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        igamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Cigamn2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        igamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}

/* Combines:
   Single precision real */

#define Csgsum2d(icontxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        sgsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Csgamx2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        sgamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Csgamn2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        sgamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}

/* Combines:
   Double precision real */

#define Cdgsum2d(icontxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        dgsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdgamx2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        dgamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Cdgamn2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        dgamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}

/* Combines:
   Single precision complex */

#define Ccgsum2d(icontxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        cgsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Ccgamx2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        cgamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Ccgamn2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        cgamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}

/* Combines:
   Double precision complex */

#define Czgsum2d(icontxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        zgsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Czgamx2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        zgamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Czgamn2d(icontxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = icontxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        zgamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}



