/*
 *             Parallel Sparse BLAS  v2.0
 *   (C) Copyright 2006 Salvatore Filippone    University of Rome Tor Vergata
 *                      Alfredo Buttari        University of Rome Tor Vergata
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions
 * are met:
 *   1. Redistributions of source code must retain the above copyright
 *      notice, this list of conditions and the following disclaimer.
 *   2. Redistributions in binary form must reproduce the above copyright
 *      notice, this list of conditions, and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *   3. The name of the PSBLAS group or the names of its contributors may
 *      not be used to endorse or promote products derived from this
 *      software without specific written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 * TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 * PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE PSBLAS GROUP OR ITS CONTRIBUTORS
 * BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 */
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
#define Cblacs_get(ictxt, what, val) \
        {i1 = ictxt; i2 = what; \
        blacs_get_(&i1, &i2,val);}
#define Cblacs_set(ictxt, what, val) \
        {i1 = ictxt; i2 = what; \
        blacs_set_(&i1, &i2, &val);}
#define Cblacs_gridinit(ictxt, order, nprow, npcol) \
        {i1 = nprow; i2 = npcol; \
        blacs_gridinit_(ictxt, order, &i1, &i2);}
#define Cblacs_gridmap(ictxt, pmap, ldpmap, nprow, npcol) \
        {i1 = ldpmap; i2 = nprow; i3 = npcol; \
        blacs_gridmap_(ictxt, pmap, &i1, &i2, &i3);}

/* Support routines:
   Destruction */

#define Cblacs_freebuff(ictxt, wait) \
        {i1 = ictxt; i2 = wait; \
        blacs_freebuff_(&i1, &i2);}
#define Cblacs_gridexit(ictxt) \
        {i1 = ictxt; \
        blacs_gridexit_(&i1);}
#define Cblacs_abort(ictxt, errornum) \
        {i1 = ictxt; i2 = errornum; \
        blacs_abort_(&i1, &i2);}
#define Cblacs_exit(doneflag) \
        {i1 = doneflag; \
        blacs_exit_(&i1);}

/* Support routines:
   Informational and Miscellaneous */

#define Cblacs_gridinfo(ictxt,nprow,npcol,myprow,mypcol) \
        {i1 = ictxt; \
        blacs_gridinfo_(&i1, nprow, npcol, myprow, mypcol);}
#define Cblacs_pnum(ictxt, prow, pcol) \
        {i1 = ictxt; i2 = prow; i3 = pcol; \
        blacs_pnum_(&i1, &i2, &i3);}
#define Cblacs_pcoord(ictxt, pnum, prow, pcol) \
        {i1 = ictxt; i2 = pnum; \
        blacs_pcoord_(&i1, &i2, prow, pcol);}
#define Cblacs_barrier(ictxt, scope) \
        {i1 = ictxt; \
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
#define Cksendid(ictxt, rdest, cdest) \
        {i1 = ictxt; i2 = rdest; i3 = cdest; \
        ksendid_(&i1, &i2, &i3);}
#define Ckrecvid(ictxt, rsrc, csrc) \
        {i1 = ictxt; i2 = rsrc; i3 = csrc; \
        krecvid_(&i1, &i2, &i3);}
#define Ckbsid(ictxt, scope) \
        {i1 = ictxt; \
        kbsid_(&i1, scope);}
#define Ckbrid(ictxt, scope, rsrc, csrc) \
        {i1 = ictxt; i2 = rsrc; i3 = csrc; \
        kbrid_(&i1, scope, &i2, &i3);}

/* Point to Point :
   Integer */

#define Cigesd2d(ictxt, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        igesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cigerv2d(ictxt, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        igerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Citrsd2d(ictxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        itrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Citrrv2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        itrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Point to Point :
   Single precision real */

#define Csgesd2d(ictxt, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        sgesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Csgerv2d(ictxt, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        sgerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cstrsd2d(ictxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        strsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Cstrrv2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        strsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Point to Point :
   Double precision real */

#define Cdgesd2d(ictxt, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        dgesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdgerv2d(ictxt, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        dgerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdtrsd2d(ictxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        dtrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdtrrv2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        dtrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Point to Point :
   Single precision complex */

#define Ccgesd2d(ictxt, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        cgesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Ccgerv2d(ictxt, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        cgerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cctrsd2d(ictxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        ctrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Cctrrv2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        ctrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Point to Point :
   Double precision complex */

#define Czgesd2d(ictxt, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        zgesd2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Czgerv2d(ictxt, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        zgerv2d_(&i1, &i2, &i3, A, &i4, &i5, &i6);}
#define Cztrsd2d(ictxt, uplo, diag, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        ztrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}
#define Cztrrv2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        ztrsd2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Integer */

#define Cigebs2d(ictxt, scope, top, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        igebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Cigebr2d(ictxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        igebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Citrbs2d(ictxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        itrbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Citrbr2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        igebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Single precision real */

#define Csgebs2d(ictxt, scope, top, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        sgebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Csgebr2d(ictxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        sgebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cstrbs2d(ictxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        strbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Cstrbr2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        sgebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Double precision real */

#define Cdgebs2d(ictxt, scope, top, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        dgebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Cdgebr2d(ictxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        dgebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdtrbs2d(ictxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        dtrbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Cdtrbr2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        dgebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Single precision complex */

#define Ccgebs2d(ictxt, scope, top, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        cgebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Ccgebr2d(ictxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        cgebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cctrbs2d(ictxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        ctrbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Cctrbr2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        cgebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Broadcasts :
   Double precision complex */

#define Czgebs2d(ictxt, scope, top, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        zgebs2d_(&i1, scope, top, &i2, &i3, A, &i4);}
#define Czgebr2d(ictxt, scope, top, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        zgebr2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cztrbs2d(ictxt, scope, top, uplo, diag, m, n, A, lda) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; \
        ztrbs2d_(&i1, scope, top, uplo, diag, &i2, &i3, A, &i4);}
#define Cztrbr2d(ictxt, uplo, diag, m, n, A, lda, rsrc, csrc) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rsrc; i6 = csrc; \
        zgebr2d_(&i1, uplo, diag, &i2, &i3, A, &i4, &i5, &i6);}

/* Combines:
   Integer */

#define Cigsum2d(ictxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        igsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cigamx2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        igamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Cigamn2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        igamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}

/* Combines:
   Single precision real */

#define Csgsum2d(ictxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        sgsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Csgamx2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        sgamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Csgamn2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        sgamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}

/* Combines:
   Double precision real */

#define Cdgsum2d(ictxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        dgsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Cdgamx2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        dgamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Cdgamn2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        dgamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}

/* Combines:
   Single precision complex */

#define Ccgsum2d(ictxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        cgsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Ccgamx2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        cgamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Ccgamn2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        cgamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}

/* Combines:
   Double precision complex */

#define Czgsum2d(ictxt, scope, top, m, n, A, lda, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = rdest; i6 = cdest; \
        zgsum2d_(&i1, scope, top, &i2, &i3, A, &i4, &i5, &i6);}
#define Czgamx2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        zgamx2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}
#define Czgamn2d(ictxt, scope, top, m, n, A, lda, RA, CA, RCflag, rdest, cdest) \
        {i1 = ictxt; i2 = m; i3 = n; i4 = lda; i5 = RCflag; i6 = rdest; i7 = cdest; \
        zgamn2d_(&i1, scope, top, &i2, &i3, A, &i4, RA, CA, &i5, &i6, &i7);}



