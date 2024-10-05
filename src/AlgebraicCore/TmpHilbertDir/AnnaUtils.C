//   Copyright (c)  2006  Anna Bigatti,  John Abbott

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include <stdlib.h>
#include "AnnaUtils.h"

// This file deliberately left empty!


/*  ==========  memory management  ==========  */

// // #ifdef COCOA_WITH_RUM
 
// void rum_init_all(void);
// /* initialize all the rum-stacks 0 slots */
 
// void rum_init(int rum, int size);
// /* initialize the rum-stack rum for size slots */
 
// void * rum_malloc(int size);
 
// void rum_free(int size, void * p);

// /*********************  rum.c   **********************/

 
// #define RUM_N_STACKS 100
 
// typedef struct rum_aux {
//   void ** slots;
//   int top;
//   int len;
// } rum_stack;
 
// rum_stack rums[RUM_N_STACKS];
 
// #define RUM_MAX_ALLOCATION 100000
 
// short rum_status[RUM_MAX_ALLOCATION]; 
// /* very large number: if the program attempt to allocate more then
//    RUM_MAX_ALLOCATION bytes at a time, you cannot use RUM */
 
// #define rum_NORMAL 0
// #define rum_EMPTY 1
// #define rum_FULL 2
// #define rum_FULL_AND_EMPTY 3
 
// void rum_empty(int rum);
 
// /* TODO: examinate the possibility to use shifts */
 
// #define rum_is_empty(rum) (rum_status[rum] & rum_EMPTY)
// #define rum_is_full(rum) (rum_status[rum] & rum_FULL)
 
// void rum_init_all(void){
//   register int i;
//   for (i=0;i<RUM_N_STACKS;++i) {
//     rums[i].len=0;
//     rums[i].top=-1; }
//   for (i=0;i<RUM_MAX_ALLOCATION;++i)
//     rum_status[i]=rum_FULL_AND_EMPTY;
// }
 
// void rum_init(int rum, int size) {
//   if (rum<RUM_N_STACKS) {
//     int i;
//     if (rums[rum].len!=0) rum_empty(rum);
//     rums[rum].slots=(void **)malloc(size*sizeof(void *));
//     rums[rum].len=size; 
//     rums[rum].top=-1;
//     rum_status[rum]=rum_EMPTY;
// /*
//     for (i=1;i<=size/2;++i) rum_free(rum,(void *)malloc(rum));
// */
// }
// }
 
// void * rum_malloc(int size)
// {
//   if (rum_is_empty(size)) return (void *)malloc(size);
//   if (rums[size].top==0) {
//     register int i;
//     for (i=0;i<30;++i) rum_free(size,(void *)malloc(size));
//     /* rum_status[size]=rum_EMPTY; */
//     /* printf("Room # %d was empty.\n",size);*/ }
//   else
//     rum_status[size]=rum_NORMAL;
//   return rums[size].slots[rums[size].top--];
// }
 
 
// void rum_free(int size, void * p)
// {
//   if (rum_is_full(size)) free(p);
//   else {
//     if (rums[size].top==rums[size].len-2) {
//       rum_status[size]=rum_FULL;
//       /* printf("Room # %d is full.\n",size);*/ }
//     else
//       rum_status[size]=rum_NORMAL;
//     rums[size].slots[++rums[size].top]=p;
//   }
// }
 
 
// /* auxiliary functions */
 
// void rum_empty(int rum)
// {
//   register int i;
//   for (i=0;i<=rums[rum].top;++i) free(rums[rum].slots[i]);
//   free(rums[rum].slots);
//   rum_status[rum]=rum_FULL_AND_EMPTY;
// }
  
/*********************************************/
// #endif

// void *myrealloc(void *oldptr, size_t old_size, size_t new_size,char * m)
// {
//   void *ret;
  
//   ret = realloc (oldptr, new_size);
//   if (ret == 0)
//     {
//       /* perror ("cannot allocate in libmp"); */
//       printf("* myrealloc *");
//       abort ();
//     }

//   return ret;
// }



// void *mycalloc(int len, size_t size, char * m) {
//   void *p = calloc(len, size);
//   if (p == NULL) {
//   		printf("* mycalloc *");
// 		abort();
//   }


//   return p;
// }

