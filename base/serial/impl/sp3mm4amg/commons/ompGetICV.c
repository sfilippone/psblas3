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
#include <stdlib.h>
#include <stdio.h>
#include <omp.h>
#include <config.h>

char* SCHEDULES[]={"OMP_SCHED_STATIC","OMP_SCHED_DYNAMIC","OMP_SCHED_GUIDED","OMP_SCHED_AUTO"};
void ompGetRuntimeSchedule(int* kind_chunk_monotonic){
	/*
	 * export OMP_SCHEDULE="[modifier:]kind[, chunk]"
	 * typedef enum omp_sched_t {
	 * //schedule kinds
	 * omp_sched_static     = 0x1,
	 * omp_sched_dynamic    = 0x2,
	 * omp_sched_guided     = 0x3,
	 * omp_sched_auto       = 0x4,
	 * //schedule modifier  //TODO API>=5.0
	 * omp_sched_monotonic  = 0x80000000u   //TODO OMP >=5
	 * } omp_sched_t;
	 */
	omp_sched_t k,kind; int chunk_size,monotonic=0;
	omp_get_schedule(&kind,&chunk_size);
	k=kind; //[monotonic OFF]
	#if _OPENMP >= 201811   //OMP_SCHED_SCHEDULE modifier from 5.0
	monotonic = omp_sched_monotonic & kind;
	if(monotonic)   k = kind - omp_sched_monotonic;
	#endif
	printf("omp sched gather:\tkind:%s\tomp chunkSize:%d\tmonotonic:%s\tfairChunkFolding:%d\n",
	  SCHEDULES[k-1],chunk_size,monotonic?"Y":"N",FAIR_CHUNKS_FOLDING);
	if (kind_chunk_monotonic){
		kind_chunk_monotonic[0] = k;
		kind_chunk_monotonic[1] = chunk_size;
		kind_chunk_monotonic[2] = monotonic;
	}
}

float ompVersionMacroMap(){
	switch ( _OPENMP ){
		case 200505:	return 2.5; 
		case 200805:	return 3.0;
		case 201107:	return 3.1;
		case 201307:	return 4.0;
		case 201511:	return 4.5;
		case 201811:	return 5.0;
		case 202011:	return 5.1;
	}
}
//WRAPPER to print all ICV vars
#ifdef OMP_GET_ICV_MAIN
int main(){
#else
void ompGetAllICV(){
#endif
	printf("export OMP_DISPLAY_ENV=VERBOSE for full ICV details\n"); 
	printf("omp MAX THREADS USABLE\t%d\n",omp_get_max_threads());
	ompGetRuntimeSchedule(NULL);
	printf("omp API version:\t %1.1f\n",ompVersionMacroMap()); 
}
