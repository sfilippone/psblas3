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
#ifndef MACROS
#define MACROS

#include <stdio.h>
#include <stdlib.h>

///aux macro-functions
#define	ABS(a)				        ((a) > 0   ? (a) : -(a))
#define	MIN(a,b)			        ((a) < (b) ? (a) : (b))
#define MAX(a,b)			        ((a) > (b) ? (a) : (b))
#define AVG(a,b)                    ((a)/2 + (b)/2 + ((a)%2+(b)%2)/2)
#define SWAP(a,b)                   (a)=(a)^(b);(b)=(b)^(a);(a)=(a)^(b)
#define IN_RANGE(i,start,end)		( (start) <= (i) && (i) <= (end) )

//ceil(x/y) with integers
#define INT_DIV_CEIL(x,y)		    ( ( (x) - 1) / (y) + 1 )
//2D ROW MAJOR indexing wrap compute
#define IDX2D(i,j,nCols)            ((j) + (i)*(nCols))

///distribuite reminder @rem in group givin an extra +1 to the first @rem
#define UNIF_REMINDER_DISTRI(i,div,rem) \
    ( (div) + ( (i) < (rem) ? 1 : 0 ) )
#define UNIF_REMINDER_DISTRI_STARTIDX(i,div,rem) \
    ( (i) * (div) + MIN( (i),(rem) ) )
#define UNIF_REMINDER_DISTRI_ENDIDX(i,div,rem) \
    ( (i+1) * (div) + MIN( (i),(rem) ) )
//shorter name alias
#define unifRemShare(i,div,rem)			UNIF_REMINDER_DISTRI(i,div,rem)
#define unifRemShareStart(i,div,rem)	UNIF_REMINDER_DISTRI_STARTIDX(i,div,rem)
#define unifRemShareEnd(i,div,rem)		UNIF_REMINDER_DISTRI_ENDIDX( (i) , (div) , (rem) )
#define unifRemShareBlock(i,div,rem)	unifRemShareStart(i,div,rem), unifRemShare(i,div,rem)

#define STATIC_ARR_ELEMENTS_N(arr)  (sizeof( (arr) ) / (sizeof(*(arr))))  
////STRING UTILS
#define _STR(s) #s
#define STR(s) _STR(s)
///CONCATENATE
//Concatenate preprocessor tokens A and B WITHOUT   expanding macro definitions
#define _CAT(a,b)		a ## b
#define _CAT3(a,b,c)	a ## b ## c
#define _CAT4(a,b,c,d)	a ## b ## c ## d
//Concatenate preprocessor tokens A and B           EXPANDING macro definitions
#define CAT(a,b)		_CAT(a,b)
#define CAT3(a,b,c)		_CAT3(a,b,c)
#define CAT4(a,b,c,d)	_CAT4(a,b,c,d)

#define _UNDEF			_
////PRINTS
#define CHIGHLIGHT                  "\33[1m\33[92m"
#define CCC                         CHIGHLIGHT
#define CHIGHLIGHTERR               "\33[31m\33[1m\33[44m"
#define CCCERR                      CHIGHLIGHTERR
#define CEND                        "\33[0m"
#define hprintsf(str,...)           printf( CHIGHLIGHT str CEND,__VA_ARGS__ ) 
#define hprintf(str)                printf( CHIGHLIGHT str CEND) 
#define ERRPRINTS(str,...)          fprintf( stderr, CHIGHLIGHTERR str CEND,__VA_ARGS__ )
#define ERRPRINT(str)               fprintf( stderr, CHIGHLIGHTERR str CEND )

#include <assert.h> 

///aux types
typedef unsigned char  uchar;
typedef unsigned short ushort;
typedef unsigned int   uint;
typedef unsigned long  ulong;
typedef char bool;
#define TRUE    (!0)
#define FALSE   0
#define true    TRUE
#define false   FALSE
#define T   TRUE
#define F   FALSE
//smart decimal type custom precision def build macro _DECIMAL_TRGT_PREC 
#ifndef _DECIMAL_TRGT_PREC
//dflt floating point precision & formatting chars
	#define _DECIMAL_TRGT_PREC	double
	#define _DECIMAL_TRGT_PREC_PR 	"%lf"
#else 
    //TODO SELECT WITH RESPECT TO THE EXPORTED TARGET DECIMAL TYPE
	#define _DECIMAL_TRGT_PREC_PR 	"%f"
#endif
typedef _DECIMAL_TRGT_PREC	decimal;


///EXTRA INCLUDE    --- cuda 
///assertionn are disabled at compile time by defining the NDEBUG preprocessor macro before including assert.h	s
//#ifdef ASSERT 	#include <assert.h> #endif

#endif 	//MACROS
