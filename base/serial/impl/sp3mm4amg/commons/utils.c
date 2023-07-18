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

//various aux functions
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
//#include <stddef.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <errno.h>
#include <limits.h>

#include "macros.h"
#include "sparseMatrix.h"
#include "utils.h"


int urndFd; //will point to urandom device file

///IO
//UNBUFFERED IO
int write_wrap(int fd,void* src,size_t count){
	ssize_t wr;
	size_t written=0;
	while (written < count){
		wr=write(fd,src+written,count-written);
		if (wr<0){
			perror("write");
			return wr;
		}
		written+=wr;
	}
	return 0;
}

int read_wrap(int fd,void* dst,size_t count){
	ssize_t rd;
	size_t readed=0;
	while (readed < count){
		rd=read(fd,dst+readed,count-readed);
		if (rd<0){
			perror("read");
			return rd;
		}
		readed+=rd;
	}
	return 0;
}

int readALL_wrap(int fd,char** dst,size_t* count){
	ssize_t rd=!0; //to allow count > fsize
	size_t readed=0;
	char allocated=0;   //flag if required *dst allocation
	if (!(*dst)){  //allocate dst buffer of same size of file
		off_t seekCurr=lseek(fd,0,SEEK_CUR);
		off_t fsize=lseek(fd,0,SEEK_END);
		if( seekCurr==-1 || fsize==-1 || lseek(fd,seekCurr,SEEK_SET)==-1){
			perror("lseek");
			return EXIT_FAILURE;
		}
		*count=fsize;
		if (!(*dst=malloc(fsize))){
			fprintf(stderr,"malloc read_wrap file size buf error\n");
			return EXIT_FAILURE;
		}
		allocated=!0;
	}
	//read loop
	while (readed < *count && rd > 0){
		rd=read(fd,(*dst)+readed,*count-readed);
		if (rd<0){
			perror("read");
			if (allocated) free(*dst);
			return rd;
		}
		readed+=rd;
	}
	if (readed < *count) (*dst)[readed]='\0';	//TODO NEEDED?
	return EXIT_SUCCESS;
}
int init_urndfd(){ // wrap init urndFd
	if((urndFd=open(DRNG_DEVFILE,O_RDONLY))<0){
			perror("open DRNG_DEVFILE");
			return EXIT_FAILURE;
	}
	return EXIT_SUCCESS;
}
int createNewFile(char* const outFpath){
	int mode=S_IRWXU;
	errno = 0;
	int outFd=open(outFpath, O_WRONLY | O_CREAT | O_TRUNC, mode);
	if (errno==EEXIST)      outFd=open(outFpath, O_WRONLY | O_TRUNC, mode);
	if (outFd<0)            perror("open outFd failed ");
	return outFd;
}
///BUFFERED IO 
//TODO series of 0 returns.... fix
int fread_wrap(FILE* fp,void* dst,size_t count){
	int rd=0,toRead;
	size_t readed=0;
	while (readed < count){
		toRead = count - readed;
		rd=fread(dst+readed,1,toRead,fp);
		if (rd != toRead){ //TODO SHORT ITEM COUNT
			if (feof(fp))	return EOF;
			else if (ferror(fp)){
				ERRPRINT("fread_wrap errd\n");
				return -2;
			}
			//else ERRPRINT("W** SHORT ITEM COUNT RETURN MEANS!?"); //TODO OO
		}
		readed+=rd;
	}
	return rd;
}
///STRUCTURED DATA IO
int writeDoubleVector(char* fpath,double* v,idx_t size){
	int fd,out=EXIT_FAILURE;
	if ( (fd = createNewFile(fpath) ) < 0 ) goto _end;
	write_wrap(fd,v,size * sizeof(*v)); //raw write of vector
	out = EXIT_SUCCESS;
	DEBUG printf("written double vector into %s RAW of %d elements\n",fpath,size);

	_end:
	if (close(fd) == EOF)  perror("close errd\n");
	return out;
}

///STRUCTURED DATA IO -- BUFFERED: FSCANF - FPRINTF
int writeDoubleVectorAsStr(char* fpath,double* v,idx_t size){
	int out=EXIT_FAILURE;
	FILE* fp = fopen(fpath,"w");
	if (!fp){
		perror("fopen vector file write");
		return EXIT_FAILURE;
	}
	for (idx_t i=0; i<size; i++){
		if (fprintf(fp,DOUBLE_STR_FORMAT,v[i]) < 0){
			ERRPRINT("fprintf to out vector file errd\n");
			goto _end;
		} 
	}
	out = EXIT_SUCCESS;
	DEBUG printf("written vector into %s of %d elements\n",fpath,size);

	_end:
	if (fclose(fp) == EOF)  perror("fclose errd\n");
	return out;

}

double* readDoubleVector(char* fpath,idx_t* size){
	int rd,hleft;
	double *out,*tmp;
	idx_t i=0,vectorSize=RNDVECTORSIZE;
	if (*size)   vectorSize = *size;
	FILE* fp = fopen(fpath,"r");
	if (!fp){
		perror("fopen vector file");
		return NULL;
	}
	if (!(out = malloc(vectorSize * sizeof(*out)))){ 
		ERRPRINT("vector read malloc fail for file\n");
		goto _err;
	}
	while (1){
		if (i >= vectorSize ){ //realloc the array
			vectorSize *= VECTOR_STEP_REALLOC;
			if (!(tmp=realloc(out,vectorSize*sizeof(*out)))){
				ERRPRINTS("realloc errd to ~~ %d MB\n",vectorSize >> 20);
				goto _err;
			}
			out = tmp;
			DEBUG   printf("reallocd to ~~ %d MB\n",vectorSize >> 20);
		}
		idx_t toRead = MIN(VECTOR_READ_BLOCK,(vectorSize-i));
		if((rd = fread_wrap(fp,out + i,sizeof(*out)*toRead)) == -2){
			ERRPRINT("fread_wrap errd\n");
			goto _err;
		}
		if(feof(fp))	//TODO rd == EOF not work always != W*F???
			break;
		i += rd/sizeof(out);
		if( (hleft = rd % sizeof(out)) ){
			DEBUG hprintf("half left double read... rounding down fptr\n");
			if(fseek(fp,-hleft,SEEK_CUR)){
				perror("fseek in readDoubleVector");
				goto _err;
			}
		}
	}
	//REALLOC THE ARRAY TO THE FINAL SIZE
	assert( i > 0 );
	if (!(tmp = realloc(out,*size*sizeof(*out)))){
		ERRPRINT("realloc errd\n");
		goto _err;
	}
	out = tmp;
	DEBUG printf("readed vector from %s of %d elements\n",fpath,*size);
	goto _free;

	_err:
	free(out);
	out = NULL;
	_free:
	if (fclose(fp) == EOF)  perror("fclose errd\n");
	return out;
}

double* readDoubleVectorStr(char* fpath,idx_t* size){
	int fscanfOut;
	double *out,*tmp;
	idx_t i=0,vectorSize=RNDVECTORSIZE;
	if (*size)   vectorSize = *size;
	FILE* fp = fopen(fpath,"r");
	if (!fp){
		perror("fopen vector file");
		return NULL;
	}
	if (!(out = malloc(vectorSize * sizeof(*out)))){ 
		ERRPRINT("vector read malloc fail for file\n");
		goto _err;
	}
	while (1){
		if (i >= vectorSize ){ //realloc the array
			vectorSize *= VECTOR_STEP_REALLOC;
			if (!(tmp=realloc(out,vectorSize*sizeof(*out)))){
				ERRPRINTS("realloc errd to ~~ %d MB\n",vectorSize >> 20);
				goto _err;
			}
			out = tmp;
			DEBUG   printf("reallocd to ~~ %d MB\n",vectorSize >> 20);
		}
		fscanfOut = fscanf(fp,DOUBLE_STR_FORMAT,out + i++ );
		if ( fscanfOut == EOF && ferror(fp)){
			perror("invalid fscanf");
			goto _err;
		}
		if ( fscanfOut != 1 || fscanfOut == EOF )   break;  //end of vector
	}
	//REALLOC THE ARRAY TO THE FINAL SIZE
	assert( i > 0 );
	*size = --i;
	if (!(tmp = realloc(out,*size*sizeof(*out)))){
		ERRPRINT("realloc errd\n");
		goto _err;
	}
	out = tmp;
	DEBUG printf("readed vector from %s of %d elements\n",fpath,*size);
	goto _free;

	_err:
	free(out);
	out = NULL;
	_free:
	if (fclose(fp) == EOF)  perror("fclose errd\n");
	return out;
}

///CFG-AUX
int getConfig(CONFIG* conf){
	int changes=EXIT_FAILURE;
	char *varVal,*ptr;
	idx_t val;
	if ((varVal = getenv(GRID_ROWS))){
		val=strtoul(varVal,&ptr,10);
		if (ptr==varVal || val>= INT_MAX){
			perror("strtol errd");
		} else {
			conf->gridRows = val;
		}
		changes = EXIT_SUCCESS;
	}
	if ((varVal = getenv(GRID_COLS))){
		val=strtoul(varVal,&ptr,10);
		if (ptr==varVal || val>= INT_MAX){
			perror("strtol errd");
		} else {
			conf->gridCols = val;
		}
		changes = EXIT_SUCCESS;
	}
	return changes;
}

/////LIB-SORTING -- WRAPPERS
//comparing functions
static inline int cmp_idx_t(const void* a, const void*b){
	idx_t aa=*((idx_t*) a), bb = *((idx_t*) b);
	return aa==bb?0:aa>bb?1:-1;
}
static inline int cmpidx_t(const void* a, const void*b){
	idx_t aa=*((idx_t*) a), bb = *((idx_t*) b);
	return aa==bb?0:aa>bb?1:-1;
}
static inline int cmpuint(const void* a, const void*b){
	uint aa=*((uint*) a), bb = *((uint*) b);
	return aa==bb?0:aa>bb?1:-1;
}
static inline int cmpRbNode(const void* a, const void* b){
	rbNode *aa=(rbNode*) a, *bb = (rbNode*) b;
	return cmp_idx_t(&aa->key,&bb->key);
}
//sorting functions 
void sort_idx_t(idx_t* arr, idx_t len){
	qsort(arr,len,sizeof(*arr),cmp_idx_t);
}
void sortidx_t(idx_t* arr, idx_t len){
	qsort(arr,len,sizeof(*arr),cmpidx_t);
}
void sortuint(uint* arr, uint len){
	qsort(arr,len,sizeof(*arr),cmpuint);
}
void sortRbNode(rbNode* arr, idx_t len){
	qsort(arr,len,sizeof(*arr),cmpRbNode);
}

////Tests AUX
inline void assertArrNoRepetitions(idx_t* arrSorted, idx_t arrLen){
	if (arrLen > 0 )	return;
	for (idx_t i=1,last=arrSorted[0]; i<arrLen; last = arrSorted[i++]) 
		assert( arrSorted[i] != last );
}
///MATH UTILS

static inline int rndDouble_sinAll(double* d){
	if(read_wrap(urndFd,(void*) d,sizeof(*d))){
		ERRPRINT("read_wrap failed to read rnd double\n");
		return EXIT_FAILURE;
	}
	*d = sin(*d) * MAXRND;
	return EXIT_SUCCESS;
}
long _rndHold;  //permanent storage of rnd longs
static inline int rndDouble_sinDecimal(double* d){
	if(read_wrap(urndFd,(void*) &_rndHold,sizeof(_rndHold))){
		ERRPRINT("read_wrap failed to read holder for rnd double\n");
		return EXIT_FAILURE;
	}
	*d = (_rndHold % (idx_t) ceil(MAXRND)) + sin(_rndHold);
	return EXIT_SUCCESS;
}
   
void statsAvgVar(double* values,uint numVals, double* out){
	double sum=0,sumSquare=0;
	for (uint i=0;  i<numVals;  i++){
		sum += values[i];
		sumSquare += values[i]*values[i];
	}
	out[0]  =   sum/numVals;							//AVG
	out[1]  =   sumSquare/numVals - (out[0] * out[0]);  //VAR
}

/// MATRIX - VECTOR COMPUTE UTILS
int fillRndVector(idx_t size, double* v){
	for( idx_t x=0; x<size; ++x ){
		if(rndDouble_sinAll( v+x )) return EXIT_FAILURE;
		#ifdef RNDVECTMIN
		v[x] += RNDVECTMIN;
		#endif
	}
	return EXIT_SUCCESS;
}

//convention true result in @a, toCheck in @b
int doubleVectorsDiff(double* a, double* b, idx_t n,double* diffMax){
	int out = EXIT_SUCCESS;
	double diff,diffAbs,_diffMax=0;
	if (diffMax)	*diffMax = 0;
	else diffMax = &_diffMax;
	for (idx_t i=0; i<n; i++){
		diff	= a[i] - b[i];
		diffAbs = ABS( diff );
		if( diffAbs > DOUBLE_DIFF_THREASH ){
			out = EXIT_FAILURE;
			ERRPRINTS("DOUBLE VECTORS DIFF: DOUBLE_DIFF_THREASH=%lf\t<\t"
				  "|%13lg| = %lf %% of @a[%d]\n",
			  	  DOUBLE_DIFF_THREASH,diff,100*diffAbs/ABS(a[i]),i);
			#ifdef DOUBLE_VECT_DIFF_EARLY_EXIT
			/*#pragma message("DOUBLE_VECT_DIFF_EARLY_EXIT: only first diff double reported")*/
			return EXIT_FAILURE;
			#endif
		}
		if ( ABS(*diffMax) < diffAbs )   *diffMax = diff;
	}
	DEBUG{
		printf("\nchecked diff %s"CEND" between 2 double vector of %d elements"
			   "\tmax diff: %le %s threash: %le\n", !out?CCC"OK":CCCERR"ERR!",
				n,*diffMax,*diffMax<DOUBLE_DIFF_THREASH?"<":">=",
			   DOUBLE_DIFF_THREASH);
		if (!*diffMax){ //self diff check uselss TODO REMOVE
			if (!memcpy(a,b,n*sizeof(*a)))
				printf("EXACT MATCHING AMONG THE 2 DOUBLE VECTORS\n!");
		}
	}
	return out;
}

void printMatrix(double* mat,idx_t m,idx_t n,char justNZMarkers){
	printf("printing matrix: %d x %d\n",m,n);
	idx_t i,j;
	for (i=0;i<m;i++){
		for (j=0;j<n;j++){
			if (justNZMarkers)  printf("%s",mat[IDX2D(i,j,n)]?".":" ");
			else				printf("%1.1lf ",mat[IDX2D(i,j,n)]);
		}
		printf("\n");
	}
	printf("\n");
}


void printVector(double* v,idx_t size){
	for( idx_t i=0;i<size;i++)   printf("%1.1lf ",v[i]);
	printf("\n");
}

////VAR -- MISC
/* search for pattern in @str from NULL terminated @patterns array
 * return pointer of pattern in @str if any, otherwise NULL
 */
static inline char* searchPatternsInStr(char** patterns,char* str){
	char* out = NULL;
	for (char** p=patterns;  !out && *p;	p++)	out = strstr(str,*p);
	return out;
}
/* search for pattern in @str from NULL terminated @patterns array
 * return pointer of pattern in @str if any, otherwise NULL
 */
static inline char* searchPatternInStrs(char* pattern,char** strings){
	char* out = NULL;
	for (char** str=strings;  !out && *str;	str++)	out = strstr(*str,pattern);
	return out;
}

inline int appendArr(idx_t val,APPENDARRAY* list){
	return 0;   //TODO
}

///DECOMPRESSION IN TMPFS  ; TODO vector not puttable here... figure out nicer way..
char* COMPRESS_EXTENSIONS[] = { ".gz", ".xz", ".bz2", ".zip", NULL};
#define STD_DECOMPR_FLAGS " -d -c "
char* DECOMPRESS_CMDS[] = {"gzip" STD_DECOMPR_FLAGS, "xz" STD_DECOMPR_FLAGS, "bzip2" STD_DECOMPR_FLAGS,
			   "unzip -c", NULL }; //zip is too old ;)

int extractInTmpFS(char* path, char* tmpFsDecompressPath){
	char* ext = searchPatternsInStr(COMPRESS_EXTENSIONS,path);
	if (!ext)   return -1;	//NOT SUPPORTED COMPRESS EXTENSION -> TRY REGULAR PATH
	//search first 2 char after dot of ext in DECOMPRESS_CMDS to get the decompress cmd
	if (strlen(ext) < 3 ){
		ERRPRINTS("NOT SUPPORTED DECOMPRESSION:\textension %s too short to be matched",ext);
		return -1;
	}
	char extStart[3];
	extStart[0] = ext[1];
	extStart[1] = ext[2];
	extStart[2] = '\0';
	char* decompressCmdBase = searchPatternInStrs(extStart,DECOMPRESS_CMDS);
	if (!decompressCmdBase){
		ERRPRINTS("NOT SUPPORTED DECOMPRESS FOR %s\n",ext);
		return -1;
	}
	uint cmdLen = strlen(decompressCmdBase) + strlen(path) + 4 + strlen(tmpFsDecompressPath);
	char* decompressCmd = alloca(cmdLen+1);
	if (snprintf(decompressCmd,cmdLen+1,"%s %s > %s",
		  decompressCmdBase,path,tmpFsDecompressPath) < 0){
		ERRPRINT("extractInTmpFS, snprintf errd\n");
	}
	VERBOSE printf("decompressing %s --> %s\ncmd:\t%s\n",path,tmpFsDecompressPath,decompressCmd);
	return system(decompressCmd);
}
