/*
 *              Sp3MM_for_AlgebraicMultiGrid
 *    (C) Copyright 2021-2022
 *        Andrea Di Iorio      
 * 
 *  Redistribution and use in source and binary forms, with or without
 *  modification, are permitted provided that the following conditions
 *  are met:
 *    1. Redistributions of source code must retain the above copyright
 *       notice, this list of conditions and the following disclaimer.
 *    2. Redistributions in binary form must reproduce the above copyright
 *       notice, this list of conditions, and the following disclaimer in the
 *       documentation and/or other materials provided with the distribution.
 *    3. The name of the Sp3MM_for_AlgebraicMultiGrid or the names of its contributors may
 *       not be used to endorse or promote products derived from this
 *       software without specific written permission.
 * 
 *  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 *  ``AS IS'' AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED
 *  TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
 *  PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE Sp3MM_for_AlgebraicMultiGrid GROUP OR ITS CONTRIBUTORS
 *  BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 *  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 *  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 *  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 *  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 *  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 *  POSSIBILITY OF SUCH DAMAGE.
 */ 

#ifndef OMPCHUNKSDIVIDE
#define OMPCHUNKSDIVIDE

/*
 * chunks distribution for runtime schedule setted to dynamic, TODO guided
 * generic interface with elements, srcMatrix, configuration
 * TODO the srcMatrix is leaved for advanced chunk distribution...
 * configuration is expected to have a valid number of threadNum setted
 */
#include "config.h"
#include "omp.h"
//distribution of @rows|blocks of @matrix, exploiting @config
typedef void (CHUNKS_DISTR )           (ulong,spmat*,CONFIG*);
typedef void (*CHUNKS_DISTR_INTERF )   (ulong,spmat*,CONFIG*);

//NOOP chunks division for manual configuration via export OMP_SCHEDULE  
inline void chunksNOOP(ulong r,spmat* mat,CONFIG* cfg){ return; }
#include "ompGetICV.h"
//fair division of @r elements from matrix @mat with threads num in @cfg
inline void chunksFair(ulong r,spmat* mat,CONFIG* cfg){
	assert(cfg->threadNum > 0); //configured target thread num
	omp_sched_t k,kind; int chunk_size,chunk_size_new=0,monotonic;
	omp_get_schedule(&kind,&chunk_size);
	k = kind;
	#if _OPENMP >= 201811   //OMP_SCHED_SCHEDULE modifier from 5.0
	monotonic = omp_sched_monotonic & kind;
	if(monotonic)   k = kind-omp_sched_monotonic;
	#endif
	(void) monotonic;   //else no unused warning
	switch(k){
		case omp_sched_static :
			DEBUG   printf("static it's already fair\n");
			return;
		case omp_sched_dynamic: 
			chunk_size_new = MAX(1,r/cfg->threadNum);
			break; 
		//case omp_sched_guided :
		//case omp_sched_auto   :
		//default:
	}
	if(chunk_size == chunk_size_new)	return;
	omp_set_schedule(kind,chunk_size_new);
	VERBOSE printf("chunksFair:\tchunk adapted to %d\n",chunk_size_new);
	DEBUG   ompGetRuntimeSchedule(NULL);
}
//fair division of @r elements from matrix @mat of threads in @cfg
//subdividing the fair share with a factor of @FAIR_CHUNKS_FOLDING
inline void chunksFairFolded(ulong r,spmat* mat,CONFIG* cfg){
	assert(cfg->threadNum > 0); //configured target thread num
	omp_sched_t k,kind; int chunk_size,chunk_size_new=0,monotonic;
	omp_get_schedule(&kind,&chunk_size);
	k = kind;
	#if _OPENMP >= 201811   //OMP_SCHED_SCHEDULE modifier from 5.0
	monotonic = omp_sched_monotonic & kind;
	if(monotonic)   k = kind-omp_sched_monotonic;
	#endif
	(void) monotonic;   //else no unused warning
	switch(k){
		case omp_sched_static :
			DEBUG   printf("static it's already fair\n");
			return;
		case omp_sched_dynamic: 
			chunk_size_new = MAX(1,r/(cfg->threadNum*FAIR_CHUNKS_FOLDING));
			break; 
		//case omp_sched_guided :
		//case omp_sched_auto   :
		//default:
	}
	if(chunk_size == chunk_size_new)	return;
	omp_set_schedule(kind,chunk_size_new);
	DEBUG{  //check with ICV get
		printf("chunksFairFolded:\tchunk adapted to %d\n",chunk_size_new);
		ompGetRuntimeSchedule(NULL);
	}
}
#endif
