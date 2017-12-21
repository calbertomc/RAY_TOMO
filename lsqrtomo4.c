// *********************************************************************************************************************************
//  lsqrtomo.c
//  This routine is a modified verison from the original routine named lsqrSOL.m. lsqrSOL was developed by Chris Paige and Michael Saunders
//  The LSQR is a method to solve sparse linear equations and sparse least-square problems. You can obtain the matlab version of
//  this routine at http://www.stanford.edu/group/SOL/software/lsqr/.
//  I'm using it to solve the large sparse matrix due to the seismic tomography method.
//  This routine was adapted to C language by Carlos Alberto Chaves on 05/28/2014. This version is an improved one because I only use the elements of
//  the sparse matrix which are different of zero. See all documentation to understanding the procedure.
//  University of Sao Paulo - University of Michigan
//  carlos.chaves@iag.usp.br; cchaves@umich.edu; calbertochaves@gmail.com (main)
//  Version 1.0 - 08/02/2014
// **********************************************************************************************************************************
//Sorry, but I usually don't comment my codes. I don't have much patience neither time to do it!
//I usually don't share my codes, but if I did it, please don't share it with anyone else without my permission.
//I don't share it because I don't think so I'm a good programmer. Thus, I believe that other people can do a better job than me,
//even though all the concepts behind of this routine are right!


/*LSQR solves Ax = b or min ||b - Ax||_2 if damp = 0,*/
/*or   min ||(b) - (  A   )x||   otherwise.*/
/*          ||(0)   (damp*I) ||_2*/
/*A is an m by n matrix (ideally sparse),*/
/*or a function handle such that*/
/*   y = A(x,1) returns y = A*x   (where x will be an n-vector);*/
/*   y = A(x,2) returns y = A'*x  (where x will be an m-vector).*/

/*-----------------------------------------------------------------------*/
/* LSQR uses an iterative (conjugate-gradient-like) method.*/
/* For further information, see */
/* 1. C. C. Paige and M. A. Saunders (1982a).*/
/*    LSQR: An algorithm for sparse linear equations and sparse least squares,*/
/*    ACM TOMS 8(1), 43-71.*/
/* 2. C. C. Paige and M. A. Saunders (1982b).*/
/*    Algorithm 583.  LSQR: Sparse linear equations and least squares problems,*/
/*    ACM TOMS 8(2), 195-209.*/
/* 3. M. A. Saunders (1995).  Solution of sparse rectangular systems using*/
/*    LSQR and CRAIG, BIT 35, 588-604.*/

/* Input parameters:*/
/* m   , n     are the dimensions of A.*/
/* atol, btol  are stopping tolerances.  If both are 1.0e-9 (say),*/
/*             the final residual norm should be accurate to about 9 digits.*/
/*             (The final x will usually have fewer correct digits,*/
/*             depending on cond(A) and the size of damp.)*/
/* conlim      is also a stopping tolerance.  lsqr terminates if an estimate*/
/*             of cond(A) exceeds conlim.  For compatible systems Ax = b,*/
/*             conlim could be as large as 1.0e+12 (say).  For least-squares*/
/*             problems, conlim should be less than 1.0e+8.*/
/*             Maximum precision can be obtained by setting*/
/*             atol = btol = conlim = zero, but the number of iterations*/
/*             may then be excessive.*/
/* itnlim      is an explicit limit on iterations (for safety).*/
/* show = 1    gives an iteration log,*/
/* show = 0    suppresses output.*/

/* Output parameters:*/
/* x           is the final solution.*/
/* istop       gives the reason for termination.*/
/* istop       = 1 means x is an approximate solution to Ax = b.*/
/*             = 2 means x approximately solves the least-squares problem.*/
/* r1norm      = norm(r), where r = b - Ax.*/
/* r2norm      = sqrt( norm(r)^2  +  damp^2 * norm(x)^2 )*/
/*             = r1norm if damp = 0.*/
/* Anorm       = estimate of Frobenius norm of Abar = [  A   ].*/
/*                                                    [damp*I]*/
/* Acond       = estimate of cond(Abar).*/
/* Arnorm      = estimate of norm(A'*r - damp^2*x).*/
/* xnorm       = norm(x).*/
/* var         (if present) estimates all diagonals of (A'A)^{-1} (if damp=0)*/
/*             or (A'A + damp^2*I)^{-1} if damp > 0.*/
/*             This is well defined if A has full column rank or damp > 0.*/
/*             More precisely, var = diag(Dk*Dk'), where Dk is the n*k*/
/*             matrix of search directions after k iterations.  Theoretically*/
/*             Dk satisfies Dk'(A'A + damp^2*I)Dk = I for any A or damp.*/

/*        1990: Derived from Fortran 77 version of LSQR.*/
/* 22 May 1992: bbnorm was used incorrectly.  Replaced by Anorm.*/
/* 26 Oct 1992: More input and output parameters added.*/
/* 01 Sep 1994: Print log reformatted.*/
/* 14 Jun 1997: show  added to allow printing or not.*/
/* 30 Jun 1997: var   added as an optional output parameter.*/
/* 07 Aug 2002: Output parameter rnorm replaced by r1norm and r2norm.*/
/* 03 May 2007: Allow A to be a matrix or a function handle.*/
/* 04 Sep 2011: Description of y = A(x,1) and y = A(x,2) corrected.*/
/* 04 Sep 2011: I would like to allow an input x0.*/
/*              If damp = 0 and x0 is nonzero, we could compute*/
/*              r0 = b - A*x0, solve min ||r0 - A*dx||, and return*/
/*              x = x0 + dx.  The current updating of "xnorm" would*/
/*              give norm(dx), which we don't really need.  Instead*/
/*              we would compute xnorm = norm(x0+dx) directly.*/

/*              If damp is nonzero,  we would have to solve the bigger system*/
/*                 min ||(   r0   ) - (  A   )dx||*/
/*                     ||(-damp*x0)   (damp*I)  ||_2*/
/*              with no benefit from the special structure.*/
/*              Forget x0 for now and leave it to the user.*/

/*              Michael Saunders, Systems Optimization Laboratory,*/
/*              Dept of MS&E, Stanford University.*/
/*-----------------------------------------------------------------------*/


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#define TAM 100

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
/*Calculating the norm -L2 of a vector*/
double normvec(int n, double *vt){
	int i;
	double res=0.;
	for(i=1; i<=n; i++){
		res=res+pow(vt[i],2);
	}
	res=sqrt(res);
return res;
}

//////////
double *resthypo(double mddepth, double mdlat, double mdlon, int mhypo, double *x){
    int
    i,
    k1=1,
    k2=2,
    k3=3,
    k4=4;
    
    for(i=1;i<=(mhypo/4);i++){
        
        if(fabs(x[k1])>0){x[k1]=0.0;}
        if(fabs(x[k2])>mddepth){x[k2]=(fabs(x[k2])/x[k2])*mddepth;}
        if(fabs(x[k3])>mddepth){x[k3]=(fabs(x[k3])/x[k3])*mdlat;}
        if(fabs(x[k4])>mddepth){x[k4]=(fabs(x[k4])/x[k4])*mdlon;}
        
        k1=k1+4;
        k2=k2+4;
        k3=k3+4;
        k4=k4+4;
    }
    
    
    return x;
}


/////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
char *printsol(int istop,char *string){
    
    switch(istop){
        case 1:
            strcpy(string,"Ax - b is small enough, given atol, btol");
            break;
        case 2:
            strcpy(string,"The least-squares solution is good enough, given atol");
            break;
        case 3:
            strcpy(string,"The estimate of cond(Abar) has exceeded conlim");
            break;
        case 4:
            strcpy(string,"Ax - b is small enough for this machine");
            break;
        case 5:
            strcpy(string,"The least-squares solution is good enough for this machine");
            break;
        case 6:
            strcpy(string,"Cond(Abar) seems to be too large for this machine");
            break;
        case 7:
            strcpy(string,"The iteration limit has been reached");
            break;
    }
    return string;
    
}
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

int main(){
    
    int
    mevt,
    *Acol,
    *Arow,
    *Gcol,
    *Grow,
    *Gtcol,
    *Gtrow,
    *countG,
    *countGt,
    *countA,
    yhypo,
    ch,
    i,
    j,
    k,
    f,
    n, // parameters dimension
    m, //data dimension
    p=1,
    nlinesA=1,
    nlinesG=1,
    itnlim,
    itn=0,
    stop=0,
    istop=0;
    double
    mddepth,
    mdlat,
    mdlon,
    sum1=0,
    alpha=0,
    beta,
    theta,
    rho,
    psi,
    phi,
    tau,
    delta,
    gamma,
    rhobar,
    rhobar1,
    phibar,
    gambar,
    Anorm,
    Arnorm,
    dknorm,
    ddnorm=0.0,
    xnorm=0.0,
    xxnorm=0.0,
    bnorm,
    rnorm,
    r1norm,
    r2norm,
    rhs,
    zbar,
    Acond=0.0,
    *var,
    *A,
    *G,
    *Gt,
    *d,
    *u,
    *v,
    *x,
    *w,
    *dk,
    vec[10],
    damp,
    dampsq,
    atol,
    btol,
    ctol=0,
    rtol,
    conlim,
    r1sq,
    res1,
    res2=0.0,
    test1,
    test2,
    test3,
    cs,
    cs1,
    cs2=-1,
    sn,
    sn1,
    sn2,
    t1,
    t2,
    z;
    char string[TAM];
    
    FILE *inputdhypo;
    FILE *inputmatrix;
    FILE *inputfrechet;
    FILE *inputparameters;
    FILE *inputdata;
    FILE *outputdata;
    FILE *inputvelmodel;
    
    printf("Openning the file iplsqrtomo.txt...\n");
    inputparameters=fopen("iplsqrtomo.txt","r");
    if(inputparameters==NULL){
        printf("File does not exist or has been moved: iplsqrtomo.txt\n");
        printf("Terminating the program...\n");
        printf("Try it again!!!\n");
        fclose(inputparameters);
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    printf("Reading the file iplsqrtomo.txt...\n");
    fscanf(inputparameters,"%lf %lf %lf %d %lf %d %d %d %d %lf %lf %lf",&atol,&btol,&conlim,&itnlim,&damp,&m,&n,&mevt,&yhypo,&mddepth,&mdlat,&mdlon);
    fclose(inputparameters);
    printf("Done!\n");
    dampsq=sqrt(pow(damp,2));
    
    printf("Openning the file matrixF.txt...\n");
    inputmatrix=fopen("matrixF.txt","r");
    if(inputmatrix==NULL){
        printf("File does not exist or has been moved: matrixF.txt\n");
        printf("Terminating the program...\n");
        printf("Try it again!!!\n");
        fclose(inputmatrix);
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    do{
        ch=fgetc(inputmatrix);
        if(ch=='\n'){nlinesA++;}
    }while(ch!=EOF);
    nlinesA=nlinesA-1;
    fclose(inputmatrix);
    
    
    inputmatrix=fopen("matrixF.txt","r");
    printf("Allocating memory...\n");
    A=(double *)malloc((nlinesA+1) * sizeof(double));
    Acol=(int *)malloc((nlinesA+1) * sizeof(int));
    Arow=(int *)malloc((nlinesA+1) * sizeof(int));
    countA=(int *)malloc((n+1) * sizeof(int));
    if(A==NULL || Acol==NULL || Arow==NULL || countA==NULL){
        printf("OUT OF MEMORY\n");
        printf("Terminating the inversion routine...\n");
        fclose(inputmatrix);
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
		exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    printf("Reading the file matrixF.txt...\n");
    for(i=1;i<=nlinesA+1;i++){fscanf(inputmatrix,"%d %d %lf",&Acol[i],&Arow[i],&A[i]);}
    fclose(inputmatrix);
    printf("Done!\n");
// Set up the first vectors u and v for the bidiagonalization.
// These satisfy  beta*u = b,  alfa*v = A'u.
    k=1;
    f=1;
    for(i=k;i<=nlinesA;i++){
        if(Acol[i]==Acol[i+1]){f++;}
        else{
            countA[k]=f;
            k++;
            f=1;
        }
    }
    ///////////////////////////////////////////////////////////////////////////////////////
    printf("Openning the file Frechetvelsparse.st...\n");
    inputfrechet=fopen("Frechetvelsparse.st","r");
	if(inputfrechet==NULL){
        printf("File does not exist or has been moved: Frechetvelsparse.st\n");
        printf("Terminating the program...\n");
        printf("Try it again!!!\n");
        fclose(inputfrechet);
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    do{
        ch=fgetc(inputfrechet);
        if(ch=='\n'){nlinesG++;}
    }while(ch!=EOF);
    nlinesG=nlinesG-1;
    fclose(inputfrechet);
    
    printf("Allocating memory...\n");
	G=(double *)malloc((nlinesG+1) * sizeof(double));
    Gcol=(int *)malloc((nlinesG+1) * sizeof(int));
    Grow=(int *)malloc((nlinesG+1) * sizeof(int));
    Gt=(double *)malloc((nlinesG+1) * sizeof(double));
    Gtcol=(int *)malloc((nlinesG+1) * sizeof(int));
    Gtrow=(int *)malloc((nlinesG+1) * sizeof(int));
    countG=(int *)malloc((m+1) * sizeof(int));
    countGt=(int *)malloc((n+1) * sizeof(int));
    if(G==NULL || Gcol==NULL || Grow==NULL || Gt==NULL || Gtcol==NULL || Gtrow==NULL || countG==NULL || countGt==NULL){
        printf("OUT OF MEMORY\n");
        printf("Terminating the inversion routine...\n");
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        free(G);
        G=NULL;
        free(Gcol);
        Gcol=NULL;
        free(Grow);
        Grow=NULL;
        free(Gtcol);
        Gtcol=NULL;
        free(Gtrow);
        Gtrow=NULL;
        free(countG);
        countG=NULL;
        free(countGt);
        countGt=NULL;
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    printf("Reading the file Frechetvelsparse.st...\n");
    inputfrechet=fopen("Frechetvelsparse.st","r");
    for(i=1;i<=nlinesG;i++){
        fscanf(inputfrechet,"%d %d %lf",&Gcol[i],&Grow[i],&G[i]);
    }
    fclose(inputfrechet);
    
    k=1;
    for(i=k;i<=nlinesG;i++){
        if(Grow[i]==Grow[i+1]){f++;}
        else{
            countG[k]=f;
            k++;
            f=1;
        }
    }
    //Transpose
    for(k=0;k<n;k++){
        for(i=1;i<=nlinesG;i++){
            if(Gcol[i]==k){
                Gt[p]=G[i];
                Gtcol[p]=Grow[i];
                Gtrow[p]=Gcol[i];
                p++;
            }
        }
    }
    k=1;
    f=1;
    for(i=k;i<=nlinesG;i++){
        if(Gtrow[i]==Gtrow[i+1]){f++;}
        else{
            countGt[k]=f;
            k++;
            f=1;
        }
    }
    
    printf("Openning the file obstimec.txt...\n");
    inputdata=fopen("obstimec.txt","r");
    if(inputdata==NULL){
        printf("File does not exist or has been moved: obstimec.txt\n");
        printf("Terminating the program...\n");
        printf("Try it again!!!\n");
        fclose(inputdata);
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        free(G);
        G=NULL;
        free(Gcol);
        Gcol=NULL;
        free(Grow);
        Grow=NULL;
        free(Gtcol);
        Gtcol=NULL;
        free(Gtrow);
        Gtrow=NULL;
        free(countG);
        countG=NULL;
        free(countGt);
        countGt=NULL;
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    printf("Allocating memory...\n");
    u=(double *)malloc((n+1) * sizeof(double));
    d=(double *)malloc((m+1) * sizeof(double));
    if(u==NULL || d==NULL){
        printf("OUT OF MEMORY\n");
        printf("Terminating the inversion routine...\n");
        fclose(inputdata);
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        free(G);
        G=NULL;
        free(Gcol);
        Gcol=NULL;
        free(Grow);
        Grow=NULL;
        free(Gtcol);
        Gtcol=NULL;
        free(Gtrow);
        Gtrow=NULL;
        free(countG);
        countG=NULL;
        free(countGt);
        countGt=NULL;
        free(d);
        d=NULL;
        free(u);
        u=NULL;
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    printf("Reading the file obstimec.txt...\n");
    for(i=1;i<=m;i++){
        fscanf(inputdata,"%lf",&d[i]);
        //printf("d[%d]: %lf\n",i,d[i]);
    }
    fclose(inputdata);
    printf("Done!\n");
    
    printf("Openning the file grid3D.txt...\n");
    inputvelmodel=fopen("grid3D.txt","r");
    if(inputvelmodel==NULL){
        printf("File does not exist or has been moved: grid3D.txt\n");
        printf("Terminating the program...\n");
        printf("Try it again!!!\n");
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        free(G);
        G=NULL;
        free(Gcol);
        Gcol=NULL;
        free(Grow);
        Grow=NULL;
        free(Gtcol);
        Gtcol=NULL;
        free(Gtrow);
        Gtrow=NULL;
        free(countG);
        countG=NULL;
        free(countGt);
        countGt=NULL;
        free(d);
        d=NULL;
        free(u);
        u=NULL;
        fclose(inputvelmodel);
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    
    if(yhypo==1){
        printf("Openning the file dhypo.txt...\n");
        inputdhypo=fopen("dhypo.txt","r");
        if(inputdhypo==NULL){
            printf("File does not exist or has been moved: dhypo.txt\n");
            printf("Terminating the program...\n");
            printf("Try it again!!!\n");
            free(A);
            A=NULL;
            free(Acol);
            Acol=NULL;
            free(Arow);
            Arow=NULL;
            free(countA);
            countA=NULL;
            free(G);
            G=NULL;
            free(Gcol);
            Gcol=NULL;
            free(Grow);
            Grow=NULL;
            free(Gtcol);
            Gtcol=NULL;
            free(Gtrow);
            Gtrow=NULL;
            free(countG);
            countG=NULL;
            free(countGt);
            countGt=NULL;
            free(d);
            d=NULL;
            free(u);
            u=NULL;
            fclose(inputdhypo);
            exit(EXIT_FAILURE);
        }
        printf("Done!\n");
    }
    
    
    printf("Allocating memory...\n");
	x=(double *)malloc((n+1) * sizeof(double));
    if(x==NULL){
        printf("OUT OF MEMORY\n");
        printf("Terminating the inversion routine...\n");
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        free(G);
        G=NULL;
        free(Gcol);
        Gcol=NULL;
        free(Grow);
        Grow=NULL;
        free(Gtcol);
        Gtcol=NULL;
        free(Gtrow);
        Gtrow=NULL;
        free(countG);
        countG=NULL;
        free(countGt);
        countGt=NULL;
        free(d);
        d=NULL;
        free(u);
        u=NULL;
        if(yhypo==1){fclose(inputdhypo);}
        fclose(inputvelmodel);
		exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    printf("Allocating memory...\n");
    v=(double *)malloc((n+1) * sizeof(double));
    if(v==NULL){
        printf("OUT OF MEMORY\n");
        printf("Terminating the inversion routine...\n");
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        free(G);
        G=NULL;
        free(Gcol);
        Gcol=NULL;
        free(Grow);
        Grow=NULL;
        free(Gtcol);
        Gtcol=NULL;
        free(Gtrow);
        Gtrow=NULL;
        free(countG);
        countG=NULL;
        free(countGt);
        countGt=NULL;
        free(d);
        d=NULL;
        free(u);
        u=NULL;
        free(x);
        x=NULL;
        if(yhypo==1){fclose(inputdhypo);}
        fclose(inputvelmodel);
		exit(EXIT_FAILURE);
    }
    printf("Allocating memory...\n");
    w=(double *)malloc((n+1) * sizeof(double));
    if(w==NULL){
        printf("OUT OF MEMORY\n");
        printf("Terminating the inversion routine...\n");
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        free(G);
        G=NULL;
        free(Gcol);
        Gcol=NULL;
        free(Grow);
        Grow=NULL;
        free(Gtcol);
        Gtcol=NULL;
        free(Gtrow);
        Gtrow=NULL;
        free(countG);
        countG=NULL;
        free(countGt);
        countGt=NULL;
        free(d);
        d=NULL;
        free(u);
        u=NULL;
        free(x);
        x=NULL;
        free(v);
        v=NULL;
        if(yhypo==1){fclose(inputdhypo);}
        fclose(inputvelmodel);
		exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    printf("Allocating memory...\n");
    dk=(double *)malloc((n+1) * sizeof(double));
    var=(double *)malloc((n+1) * sizeof(double));
    if(dk==NULL){
        printf("OUT OF MEMORY\n");
        printf("Terminating the inversion routine...\n");
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        free(G);
        G=NULL;
        free(Gcol);
        Gcol=NULL;
        free(Grow);
        Grow=NULL;
        free(Gtcol);
        Gtcol=NULL;
        free(Gtrow);
        Gtrow=NULL;
        free(countG);
        countG=NULL;
        free(countGt);
        countGt=NULL;
        free(d);
        d=NULL;
        free(u);
        u=NULL;
        free(x);
        x=NULL;
        free(v);
        v=NULL;
        free(w);
        w=NULL;
        free(dk);
        dk=NULL;
        if(yhypo==1){fclose(inputdhypo);}
        fclose(inputvelmodel);
		exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    if(yhypo==1){
        printf("Reading the file dhypo.txt...\n");
        for(i=1;i<=mevt;i++){
            //fscanf(inputdhypo,"%lf",&x[i]);
            //printf("x[%d]: %lf\n",i,x[i]);
            x[i]=0.0;
        	v[i]=0.0;
        }
        fclose(inputdhypo);
        printf("Reading the file grid3D.txt...\n");
        for(i=(mevt+1);i<=n;i++){
            fscanf(inputvelmodel,"%lf",&x[i]);
            //printf("x[%d]: %lf\n",i,x[i]);
                v[i]=0.0;
        }
        fclose(inputvelmodel);
    }
    else{
        printf("Reading the file grid3D.txt...\n");
        for(i=1;i<=n;i++){
            fscanf(inputvelmodel,"%lf",&x[i]);
            v[i]=0.0;
        }
        fclose(inputvelmodel);
        
    }
    printf("Done!\n");
    if(conlim > 0){ctol=1.0/conlim;}
    
    sum1=0.0;
    k=1;
    for(i=1;i<=m;i++){
        for(j=1;j<=countG[i];j++){
            sum1=sum1+G[k]*x[Gcol[k]+1];
            k++;
        }
        d[i]=d[i]-sum1;
        sum1=0.0;
    }
    sum1=0.0;
    k=1;
    for(i=1;i<=n;i++){
        for(j=1;j<=countGt[i];j++){
            sum1=sum1+(Gt[k]*d[Gtcol[k]+1]);
            k++;
        }
        u[Gtrow[k-1]+1]=sum1;
        var[i]=0.0;
        sum1=0.0;
    }
    
    beta=normvec(n,u);
	sum1=0.0;
    k=1;
    if(beta > 0){
   		for(i=1;i<=n;i++){u[i]=(1.0/beta)*u[i];}
        for(i=1;i<=n;i++){
            for(j=1;j<=countA[i];j++){
                sum1=sum1+A[k]*u[Arow[k]+1];
                k++;
            }
            v[Acol[k-1]+1]=v[Acol[k-1]+1]+sum1;
            sum1=0.0;
        }
        alpha=normvec(n,v);
	}
	if(alpha > 0){
		for(i=1;i<=n;i++){
   			v[i]=(1.0/alpha)*v[i];
   			w[i]=v[i];
		}
	}
	Arnorm=alpha*beta;
	if (Arnorm == 0){
		printf("The exact solution is  x = 0");
		printf("Terminating the program...\n");
        free(A);
        A=NULL;
        free(Acol);
        Acol=NULL;
        free(Arow);
        Arow=NULL;
        free(countA);
        countA=NULL;
        free(G);
        G=NULL;
        free(Gcol);
        Gcol=NULL;
        free(Grow);
        Grow=NULL;
        free(Gtcol);
        Gtcol=NULL;
        free(Gtrow);
        Gtrow=NULL;
        free(countG);
        countG=NULL;
        free(countGt);
        countGt=NULL;
        free(d);
        d=NULL;
        free(u);
        u=NULL;
        free(x);
        x=NULL;
        free(v);
        v=NULL;
        free(w);
        w=NULL;
        free(dk);
        dk=NULL;
		return 0;
		exit(EXIT_FAILURE); 
	}

	rhobar=alpha;
	phibar=beta;          
	bnorm=beta;
	rnorm=beta;
	r1norm=rnorm;
	r2norm=rnorm;
//------------------------------------------------------------------
//     Main iteration loop.
//------------------------------------------------------------------
	while(stop==0){
  		itn=itn+1;
        printf("iteration: %d\n",itn);
//  Perform the next step of the bidiagonalization to obtain the
//  next beta, u, alfa, v.  These satisfy the relations
//       beta*u  =  A*v  - alfa*u,
//       alfa*v  =  A'*u - beta*v.
        sum1=0.0;
        k=1;
		for(i=1;i<=n;i++){
            for(j=1;j<=countA[i];j++){
                sum1=sum1+A[k]*v[Arow[k]+1];
                k++;
            }
            u[Acol[k-1]+1]=sum1-alpha*u[Acol[k-1]+1];
            sum1=0.0;
        }
		beta=normvec(n,u);
        //printf("beta: %lf\n",beta);
  		if(beta > 0){
			for(i=1;i<=n;i++){
    				u[i]=(1.0/beta)*u[i];
			}
            vec[1]=Anorm;
            vec[2]=alpha;
            vec[3]=beta;
            vec[4]=damp;
            Anorm=normvec(4,vec);
            sum1=0.0;
            k=1;
			for(i=1;i<=n;i++){ 				
                for(j=1;j<=countA[i];j++){
                    sum1=sum1+A[k]*u[Arow[k]+1];
                    k++;
                }
                v[Acol[k-1]+1]=sum1-beta*v[Acol[k-1]+1];
                sum1=0.0;
            }
            alpha=normvec(n,v);
            //printf("alpha: %lf\n",alpha);
            if(alpha > 0){
                for(i=1;i<=n;i++){
                    v[i]=(1.0/alpha)*v[i];
                }
            }
  		}

//  Use a plane rotation to eliminate the damping parameter.
//  This alters the diagonal (rhobar) of the lower-bidiagonal matrix.
        vec[1]=rhobar;
		vec[2]=damp;
		rhobar1=normvec(2,vec);
		cs1=rhobar/rhobar1;
		sn1=damp/rhobar1;
		psi=sn1*phibar;
		phibar=cs1*phibar;

//  Use a plane rotation to eliminate the subdiagonal element (beta)
//  of the lower-bidiagonal matrix, giving an upper-bidiagonal matrix.
		vec[1]=rhobar1;
		vec[2]=beta;
		rho=normvec(2,vec);
		cs=rhobar1/rho;
		sn=beta/rho;
		theta=sn*alpha;
		rhobar=-cs*alpha;
		phi=cs*phibar;
		phibar=sn*phibar;
		tau=sn*phi;

// Update x and w.

		t1=phi/rho;
		t2=-theta/rho;
		for(i=1;i<=n;i++){
			dk[i]=(1.0/rho)*w[i];
			x[i]=x[i]+t1*w[i];
			w[i]=v[i]+t2*w[i];
            var[i]=var[i]+(pow(dk[i],2));
		}
        
        if(yhypo==1){resthypo(mddepth,mdlat,mdlon,mevt,x);}
        
        dknorm=normvec(n,dk);
		ddnorm=ddnorm+pow(dknorm,2);
		
//  Use a plane rotation on the right to eliminate the
//  super-diagonal element (theta) of the upper-bidiagonal matrix.
//  Then use the result to estimate  norm(x).

		delta=sn2*rho;
		gambar=-cs2*rho;
		rhs=phi-delta*z;
		zbar=rhs/gambar;
		xnorm=sqrt(xxnorm+pow(zbar,2));
		vec[1]=gambar;
		vec[2]=theta;		
		gamma=normvec(2,vec);
		cs2=gambar/gamma;
		sn2=theta/gamma;
		z=rhs/gamma;
		xxnorm=xxnorm+pow(z,2);

//  Test for convergence.
//  First, estimate the condition of the matrix  Abar,
//  and the norms of  rbar  and  Abar'rbar.

		Acond=Anorm*sqrt(ddnorm);
		res1=pow(phibar,2);
		res2=res2+pow(psi,2);
		rnorm=sqrt(res1+res2);
		Arnorm=alpha*fabs(tau);

//  07 Aug 2002:
//  Distinguish between
//     r1norm = ||b - Ax|| and
//     r2norm = rnorm in current code
//            = sqrt(r1norm^2 + damp^2*||x||^2).
//     Estimate r1norm from
//     r1norm = sqrt(r2norm^2 - damp^2*||x||^2).
//  Although there is cancellation, it might be accurate enough.

		r1sq=pow(rnorm,2)-dampsq*xxnorm;
		r1norm=sqrt(fabs(r1sq));   
		if(r1sq<0){
			r1norm=-r1norm;
		}
		r2norm=rnorm;

//  Now use these norms to estimate certain other quantities,
//  some of which will be small near a solution.

		test1=rnorm /bnorm;
		test2=Arnorm/(Anorm*rnorm);
		test3=1.0/Acond;
		t1=test1/(1.0+((Anorm*xnorm)/bnorm));
		rtol=btol+((atol*Anorm*xnorm)/bnorm);

//  The following tests guard against extremely small values of
//  atol, btol  or  ctol.  (The user may have set any or all of
//  the parameters  atol, btol, conlim  to 0.)
//  The effect is equivalent to the normal tests using
//  atol = eps,  btol = eps,  conlim = 1/eps.

  		if(itn>=itnlim){
  			istop=7;
  		}
  		if(1+test3<=1){
  			istop=6;
  		}
  		if(1+test2<=1){
            istop=5;
  		}
  		if(1+t1<= 1){
  			istop=4;
  		}

// Allow for tolerances set by the user.

  		if (test3 <= ctol){istop = 3;}
  		if (test2 <= atol){istop = 2;}
  		if (test1 <= rtol){istop = 1;}

        if(istop>0){stop=1;}
        //End of iteration loop.
	}
    printf("The lsqrtomo routine has finished\n");
    printsol(istop,string);
    printf("%s\n",string);
    
    printf("Terminating the program...\n");
    
    outputdata=fopen("paramestim.txt","w");
    for(i=1;i<=n;i++){
        fprintf(outputdata, "%lf\n",x[i]);
    }
    fclose(outputdata);
    
    outputdata=fopen("stdparamestim.txt","w");
    for(i=1;i<=n;i++){
        fprintf(outputdata, "%lf\n",sqrt(var[i]));
    }
    fclose(outputdata);
    
    outputdata=fopen("calctime.txt","w");
    sum1=0.0;
    k=1;
    for(i=1;i<=m;i++){
        for(j=1;j<=countG[i];j++){
            sum1=sum1+G[k]*x[Gcol[k]+1];
            k++;
        }
        fprintf(outputdata, "%lf\n",sum1);
        sum1=0.0;
    }
    fclose(outputdata);
    
    outputdata=fopen("rnorm.txt","w");
    fprintf(outputdata, "%lf\n",r1norm);
    fclose(outputdata);
    outputdata=fopen("xnorm.txt","w");
    fprintf(outputdata, "%lf\n",xnorm);
    fclose(outputdata);
    
    
    free(A);
    A=NULL;
    free(Acol);
    Acol=NULL;
    free(Arow);
    Arow=NULL;
    free(countA);
    countA=NULL;
    free(G);
    G=NULL;
    free(Gcol);
    Gcol=NULL;
    free(Grow);
    Grow=NULL;
    free(Gtcol);
    Gtcol=NULL;
    free(Gtrow);
    Gtrow=NULL;
    free(countG);
    countG=NULL;
    free(countGt);
    countGt=NULL;
    free(d);
    d=NULL;
    free(u);
    u=NULL;
    free(x);
    x=NULL;
    free(v);
    v=NULL;
    free(w);
    w=NULL;
    free(dk);
    dk=NULL;
    free(var);
    var=NULL;
    
	return 0;
}

/*-----------------------------------------------------------------------*/
/* end function lsqrtomo.c*/
/*-----------------------------------------------------------------------*/

