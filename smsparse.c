//*********************************************************************************************************************************************
//  smsparse.c
//  Routine developed to multiply sparse matrix in order to produce the matrix F=A'A+miu1*L'L+miu2I, so that Fm=b, where b=A'*d;
//  Here, to improve the manipulation with sparse matrix I use the package CSPARSE - A Concise Sparse Matrix Package in C developed by Timothy Davis
//  CSPARSE is a C library which implements a number of direct methods for sparse linear systems, by Timothy Davis.
//  http://www.suitesparse.com
//  Refer to "Direct Methods for Sparse Linear Systems," Timothy A. Davis,
//  SIAM, Philadelphia, 2006.  No detailed user guide is included in this package; the user guide is the book itself.
//  Created by Carlos Alberto Chaves on 09/05/2015.
//  University of Sao Paulo - University of Michigan
//  carlos.chaves@iag.usp.br; cchaves@umich.edu; calbertochaves@gmail.com (main)
//  Version 1.0
//*********************************************************************************************************************************************
//Sorry, but I usually don't comment my codes. I don't have much patience neither time to do it!
//I usually don't share my codes, but if I did it, please don't share it with anyone else without my permission.
//I don't share it because I don't think so I'm a good programmer. Thus, I believe that other people can do a better job than me,
//even though all the concepts behind of this routine are right! Obviously, it is very difficult to test the routine in all possible
//situations, so if you find some bug, please email me and let me know the problem.

# include <stdlib.h>
# include <limits.h>
# include <math.h>
# include <stdio.h>
# include "csparse.h"

int main(int numargs, char *arg[]){
  
    float
    miu1,
    miu2;
    cs *A1;
    cs *A2;
    cs *C1;
    cs *C2;
    cs *Eye;
    int i;
    int m;
    cs *T1;
    cs *T2;

    
    FILE *inputfile1;
    FILE *inputfile2;

    inputfile1=fopen(arg[1],"r");
    inputfile2=fopen(arg[2],"r");
    miu1=atof(arg[3]);
    miu2=atof(arg[4]);
  //Load the triplet matrix T from standard input.
    T1=cs_load(inputfile1);
    T2=cs_load(inputfile2);
    
    
    A1=cs_triplet(T1);
    A2=cs_triplet(T2);

    cs_spfree(T1);
    cs_spfree(T2);
    
    /*
     Compute C2=A1+A2*miu1.
     */
    
    C1=cs_add (A1,A2,1.0,miu1);
//  M = number of rows of C1.
    m = C1->m;
    
    /*
    Create triplet identity matrix.
    */
    T1 =cs_spalloc(m,m,m,1,1);
    
    for (i=0; i < m; i++){
        cs_entry(T1,i,i,1);
    }
    /*
     Eye = speye ( m )
     */
    Eye = cs_triplet(T1);
    cs_spfree(T1);
    
    /*
     Compute C2 = C1 + Eye * miu2.
     */
    C2=cs_add(C1,Eye,1.0,miu2);

    cs_print (C2,0);
    
    cs_spfree(A1);
    cs_spfree(A2);
    cs_spfree(Eye);
    cs_spfree(C2);
    cs_spfree(C1);
/*
  Terminate.
*/
    
  return(0);
}
