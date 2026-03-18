#pragma once

/*
 * spGPU - Sparse matrices on GPU library.
 * 
 * Copyright (C) 2010 - 2012 
 *     Davide Barbieri - University of Rome Tor Vergata
 *
 */
 
#include "stdio.h"
#include "stdlib.h"

#ifndef _WIN32
#include <execinfo.h>

inline
void printTrace (void)
{
	void *array[32];
	size_t size;
	char **strings;
	size_t i;
     
	size = backtrace (array, 32);
	strings = backtrace_symbols (array, size);
     
	printf ("---- Obtained %zd stack frames.\n", size);
     
	for (i = 0; i < size; i++)
		printf ("%s\n", strings[i]);
		    
	free (strings);
}
#endif

inline void __assert(int e, const char* w)
{
	if (!e)
	{
		printf("%s\n",w);
		
#ifndef _WIN32
		printTrace();
#endif
		
		exit(0);
	}
}

