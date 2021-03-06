// *****************************************************************************************************************************
//  genindex1D.c
//  A routine developed to calculate the partial derivative of the traveltime with respect to velocity in order to create the
//  Frechet matrix to be used in the seismic tomography problem of realocating
//  earthquakes and determing the 3-D velocity within the Earth. Indeed, It calculate the partial derivative values and the index of the cell.
//  Output is used by either frechetvel1D.c or frechetvelsparse.c depend on what you want to do!
//  Most of this work was based on works of Inoue et al. [1990]
//  Thurber [1981], Thuber [1983], Nolet [1993] and Pulliam et al. [1993];
//  Created by Carlos Alberto Chaves on 05/18/2014.
//  University of Sao Paulo - University of Michigan
//  carlos.chaves@iag.usp.br; cchaves@umich.edu; calbertochaves@gmail.com (main)
//  Version 2.0

// Updated on 09/01/2015 ---> Now, interp1 is here. Also, to correct a bug regarding inaccuracies of calculating the path length, I use the
// time generated by taup. Sorry, this is the best that I can do this time. I need to finish my ray tracer and I no longer have this problem.
// ******************************************************************************************************************************
//Sorry, but I usually don't comment my codes. I don't have much patience neither time to do it!


#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#define TAM1 500
#define TAM2 30
#define pi 3.141592653589793 //pi number
#define rearth 6371.0 //Earth's radius


double *velblock(int nr,int nt,int np,double *rnew,double *thetanew,double *phinew,double latmin,double lonmin, double dt, double dp,double *coordvnew){
    
    int
    i,
    j,
    k,
    l,
    m,
    n,
    t,
    stop;
    double
    *rmodel,
    *latmodel,
    *lonmodel,
    *vmodel,
    r[3],
    theta[3],
    phi[3];
    char
    filevelmodel[TAM2]={"grid1D.txt"};
    
    FILE *inputvelmodel;
    FILE *inputdepth;
    ///////////////////////////////////////////////////////
    rmodel=(double *)malloc((nr+1) * sizeof(double));
    latmodel=(double *)malloc((nt+1) * sizeof(double));
    lonmodel=(double *)malloc((np+1) * sizeof(double));
    if(rmodel==NULL || latmodel==NULL || lonmodel==NULL){
        printf("OUT OF MEMORY\n");
        exit(EXIT_FAILURE);
    }
    ////////////////////////////////////////////////////////
    vmodel=(double *)malloc((nr+1) * sizeof(double));
    
    if(vmodel==NULL){
        printf("OUT OF MEMORY\n");
        exit(EXIT_FAILURE);
    }
    
    inputvelmodel=fopen(filevelmodel,"r");
    for(l=1;l<=nr-1;l++){
        fscanf(inputvelmodel,"%lf",&vmodel[l]);
    }
    fclose(inputvelmodel);
    ////////////////////////////////////////////////////////
    for(m=1;m<=np;m++){
        lonmodel[m]=lonmin+((double)(m-1)*dp);
    }
    for(m=1;m<=nt;m++){
        latmodel[m]=latmin+((double)(m-1)*dt);
    }
    inputdepth=fopen("ipdepth.txt","r");
    if(inputdepth==NULL){
        printf("File does not exist or has been moved: ipdepth.txt\n");
        fclose(inputdepth);
        exit(EXIT_FAILURE);
    }
    for(m=1;m<=nr;m++){
        fscanf(inputdepth,"%lf",&rmodel[m]);
        rmodel[m]=rearth-rmodel[m];
        
    }
    fclose(inputdepth);
    ////////////////////////////////////////////////////////
    stop=0;
    i=1;
    while(stop==0 && i<=(nr-1)){
        if((rnew[1]>=rmodel[i+1]) && (rnew[1]<=rmodel[i])){
            stop=1;
        }
        i=i+1;
    }
    r[1]=rmodel[i-1];
    r[2]=rmodel[i];
    stop=0;
    j=1;
    while(stop==0 && j<=(nt-1)){
        if((thetanew[1]>=latmodel[j]) && (thetanew[1]<=latmodel[j+1])){
            stop=1;
        }
        j=j+1;
    }
    theta[1]=latmodel[j-1];
    theta[2]=latmodel[j];
    stop=0;
    k=1;
    while(stop==0 && k<=(np-1)){
        if((phinew[1]>=lonmodel[k]) && (phinew[1]<=lonmodel[k+1])){
            stop=1;
        }
        k=k+1;
    }
    phi[1]=lonmodel[k-1];
    phi[2]=lonmodel[k];
    
    coordvnew[1]=vmodel[i-1];
    coordvnew[2]=k-1;
    coordvnew[3]=j-1;
    coordvnew[4]=i-1;
    
    free(rmodel);
    rmodel=NULL;
    free(latmodel);
    latmodel=NULL;
    free(lonmodel);
    lonmodel=NULL;
    free(vmodel);
    vmodel=NULL;
    return coordvnew;
}

int *findcell(int nrp, int nr, int nt, int np, int nor, double *rpi, double *thetapi, double *phipi, double *rmodel, double *latmodel, double *lonmodel,int *coordp){
    
    int
    i,
    j,
    k,
    stop=0,
    ip=2,
    ipi=1,
    t=1;
    double
    **r,
    **theta,
    **phi;
    
    FILE *outfile;
    FILE *outfile2;
    
    r=(double **)malloc(3 * sizeof(double *));
    theta=(double **)malloc(3 * sizeof(double *));
    phi=(double **)malloc(3 * sizeof(double *));
    for(i=0; i<=3;i++){r[i]=(double *)malloc((nrp+1) * sizeof(double));}
    for(i=0; i<=3;i++){theta[i]=(double *)malloc((nrp+1) * sizeof(double));}
    for(i=0; i<=3;i++){phi[i]=(double *)malloc((nrp+1) * sizeof(double));}
    
    
    while(ip<=(nrp-1)){
        i=1;
        stop=0;
        while(stop==0 && i<=(nr-1)){
            if((rpi[ipi]>=rmodel[i+1]) && (rpi[ipi]<=rmodel[i])){
                stop=1;
            }
            i=i+1;
        }
        r[1][t]=rmodel[i-1];
        r[2][t]=rmodel[i];
        stop=0;
        j=1;
        while(stop==0 && j<=(nt-1)){
            if((thetapi[ipi]>=latmodel[j]) && (thetapi[ipi]<=latmodel[j+1])){
                stop=1;
            }
            j=j+1;
        }
        theta[1][t]=latmodel[j-1];
        theta[2][t]=latmodel[j];
        stop=0;
        k=1;
        while(stop==0 && k<=(np-1)){
            if((phipi[ipi]>=lonmodel[k]) && (phipi[ipi]<=lonmodel[k+1])){
                stop=1;
            }
            k=k+1;
        }
        phi[1][t]=lonmodel[k-1];
        phi[2][t]=lonmodel[k];
        while((rpi[ip]>=r[2][t] && rpi[ip]<=r[1][t]) && (thetapi[ip]>=theta[1][t] && thetapi[ip]<=theta[2][t]) && (phipi[ip]>=phi[1][t] && phipi[ip]<=phi[2][t]) && (ip<=(nrp-1))){
            ip=ip+1;
        }
        if(ipi==ip){
            printf("raypath outside the target volume\n");
            printf("Recheck your data or expand the target volume to include this raypath\n");
            outfile=fopen("raypathoutmodel.txt","a+");
            outfile2=fopen("rayout.txt","w");
            fprintf(outfile,"%d\n",nor);
            fprintf(outfile2,"%d\n",nor);
            fclose(outfile);
            fclose(outfile2);
            free(r);
            r=NULL;
            free(theta);
            theta=NULL;
            free(phi);
            phi=NULL;
            exit(EXIT_FAILURE);
            
        }
        else{
            ipi=ip;
        }
        coordp[t]=ip;
        t=t+1;
    }
    
    coordp[0]=t;
    free(r);
    r=NULL;
    free(theta);
    theta=NULL;
    free(phi);
    phi=NULL;
    return coordp;
    
}

/** To interpolate data */
double *interp1(double *x, int sizex, double *y , double *xi, int sizexi, double *yi){
    // indexes
    int idyi=0,
    idx = 0;
    // After computing, to know the current step progress
    double steps;
    
    
    steps=(y[1]-y[0])/(x[1]-x[0]);
    
    // fill left points
    while((idyi < sizexi-1) && (xi[idyi] < x[0])){
        yi[idyi+1]=y[0]-(steps*(x[0]-xi[idyi])) ;
        ++idyi;
    }
    // fill point in the same interval as the original values
    while((idyi<sizexi) && (idx < sizex)){
        while((idx < sizex) && (xi[idyi] >= x[idx])){
            ++idx;
        }
        if(idx != sizex){
            steps=(y[idx]-y[idx-1])/(x[idx]-x[idx-1]);
            yi[idyi+1]=y[idx]-(steps*(x[idx]-xi[idyi])) ;
            ++idyi;
        }
    }
    
    // fill the right points
    if(idyi < sizexi){
        steps=(y[sizex-1]-y[sizex-2])/(x[sizex-1]-x[sizex-2]);
        while(idyi < sizexi){
            yi[idyi+1]=y[sizex-1]+(steps*(xi[idyi]-x[sizex-1])) ;
            ++idyi;
        }
    }
    return yi;
}


int main(){
    int
    sizecoordp,
    ch,
    i,
    j,
    k,
    m,
    n,
    p,
    nr,
    nt,
    np,
    anor, //nor = current ray number
    nrp=0,
    nrpi,
    *index,
    coordp[TAM1];
    double
    dp,
    dt,
    ttaup,
    tcalc=0,
    int1,
    latmin,
    lonmin,
    rnew[3],
    thetanew[3],
    phinew[3],
    vnew[3],
    *coordvnew,
    cospsi,
    deltas,
    dx2,
    *x1,
    *x2,
    *distp,
    *rp,
    *thetap,
    *phip,
    *rpi,
    *thetapi,
    *phipi,
    *rmodel,
    *latmodel,
    *lonmodel;
    char
    int2str[TAM2],
    int2str2[TAM2],
    fileicrayspath[TAM2]={"icraypath"}, //index cell raypath
    filerayspath[TAM2]={"raypath"};
    
    FILE *inputfile;
    FILE *inputraypath;
    FILE *inputparameters;
    FILE *outputindex;
    FILE *inputdepth;
    FILE *outfile;
    FILE *outfile2;
    /**************************************************************************/
    /*Here, the program starts!*/
    /**************************************************************************/
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inputparameters=fopen("ipgenindex1D.txt","r");
    if(inputparameters==NULL){
        printf("File does not exist or has been moved: ipgenindex1D.txt\n");
        fclose(inputparameters);
        exit(EXIT_FAILURE);
    }
    fscanf(inputparameters,"%d %d %d %lf %lf %lf %lf %lf %d %d",&nr,&nt,&np,&dt,&dp,&lonmin,&latmin,&ttaup,&anor,&nrpi);
    fclose(inputparameters);
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inputdepth=fopen("ipdepth.txt","r");
    if(inputdepth==NULL){
        printf("File does not exist or has been moved: ipdepth.txt\n");
        fclose(inputdepth);
        exit(EXIT_FAILURE);
    }
    //////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    rmodel=(double *)malloc((nr+1) * sizeof(double));
    latmodel=(double *)malloc((nt+1) * sizeof(double));
    lonmodel=(double *)malloc((np+1) * sizeof(double));
    if(rmodel==NULL || latmodel==NULL || lonmodel==NULL){
        printf("OUT OF MEMORY\n");
        exit(EXIT_FAILURE);
    }
    ///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    for(i=1;i<=np;i++){
        lonmodel[i]=lonmin+((double)(i-1)*dp);
    }
    for(i=1;i<=nt;i++){
        latmodel[i]=latmin+((double)(i-1)*dt);
    }
    for(i=1;i<=nr;i++){
        fscanf(inputdepth,"%lf",&rmodel[i]);
        rmodel[i]=rearth-rmodel[i];
        //printf("rmodel[%d]: %lf\n",i,rmodel[i]);
        
    }
    fclose(inputdepth);
    /////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    index=(int *)malloc(4 * sizeof(int));
    coordvnew=(double *)malloc(7 * sizeof(double));
    if(index==NULL || coordvnew==NULL){
        printf("OUT OF MEMORY\n");
        exit(EXIT_FAILURE);
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    printf("Calculating index of cells for each raypath...!\n");
    printf("Raypath number %d\n",anor);
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x2=(double *)malloc((nrpi+1) * sizeof(double));
    rpi=(double *)malloc((nrpi+1) * sizeof(double));
    thetapi=(double *)malloc((nrpi+1) * sizeof(double));
    phipi=(double *)malloc((nrpi+1) * sizeof(double));
    if(x2==NULL || rpi==NULL || thetapi==NULL || phipi==NULL){
        printf("OUT OF MEMORY\n");
        exit(EXIT_FAILURE);
    }
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    sprintf(int2str,"%d",anor);
    sprintf(int2str2,"%d",anor);
    strcat(filerayspath,int2str);
    strcat(fileicrayspath,int2str2);
    strcat(filerayspath,".txt");
    strcat(fileicrayspath,".txt");
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inputfile=fopen(filerayspath,"r");
    do{
        ch = fgetc(inputfile);
        if(ch == '\n'){nrp++;}
    }while (ch != EOF);
    nrp=nrp-1;
    fclose(inputfile);//close the file
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    inputfile=fopen(filerayspath,"r");
    distp=(double *)malloc((nrp+1) * sizeof(double));
    rp=(double *)malloc((nrp+1) * sizeof(double));
    thetap=(double *)malloc((nrp+1) * sizeof(double));
    phip=(double *)malloc((nrp+1) * sizeof(double));
    for(i=0;i<nrp;i++){
        fscanf(inputfile,"%lf %lf %lf %lf",&distp[i],&rp[i],&thetap[i],&phip[i]);
        rp[i]=rp[i]/rearth;
        thetap[i]=thetap[i]*(pi/180.0);
        phip[i]=phip[i]*(pi/180.0);
    }
    fclose(inputfile);//close the file
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    x1=(double *)malloc((nrp+1) * sizeof(double));
    for(i=0;i<=nrp;i++){
        x1[i]=1.0+(1.0*i);
    }
    dx2=(double)(nrp-1.0)/(nrpi-1.0);
    for(i=0;i<=nrpi;i++){
        x2[i]=1.0+(dx2*i);
    }
    interp1(x1,nrp,rp,x2,nrpi,rpi);
    interp1(x1,nrp,thetap,x2,nrpi,thetapi);
    interp1(x1,nrp,phip,x2,nrpi,phipi);
    for(i=1;i<=nrpi;i++){
        rpi[i]=rpi[i]*rearth;
        thetapi[i]=thetapi[i]*(180.0/pi);
        phipi[i]=phipi[i]*(180.0/pi);
    }
    
    /////////////////////////////////////////////////////////////////////////////
    findcell(nrpi,nr,nt,np,anor,rpi,thetapi,phipi,rmodel,latmodel,lonmodel,coordp);
    sizecoordp=coordp[0];
    i=1;
    j=1;
    k=coordp[1]-1;
    /////////////////////////////////////////////////////////////////////////////
    while (i<=(sizecoordp-1)){
        rnew[1]=rpi[j];
        thetanew[1]=thetapi[j];
        phinew[1]=phipi[j];
        rnew[2]=rpi[k];
        thetanew[2]=thetapi[k];
        phinew[2]=phipi[k];
        velblock(nr,nt,np,rnew,thetanew,phinew,latmin,lonmin,dt,dp,coordvnew);
        index[1]=(int)coordvnew[2];
        index[2]=(int)coordvnew[3];
        index[3]=(int)coordvnew[4];
        int1=1.0/pow(coordvnew[1],2);
        deltas=0.0;
        if(i<(sizecoordp-1)){
            for(m=j;m<=k-1;m++){
                cospsi=(sin(thetapi[m]*(pi/180.0))*sin(thetapi[m+1]*(pi/180.0)))+(cos(thetapi[m]*(pi/180.0))*cos(thetapi[m+1]*(pi/180.0)))*cos((phipi[m+1]-phipi[m])*(pi/180.0));
                deltas=deltas+sqrt(fabs(pow(rpi[m]/rearth,2)+pow(rpi[m+1]/rearth,2)-(2.0*((rpi[m]*rpi[m+1])/(rearth*rearth))*cospsi)));
            }
        }
        else{
            for(m=j;m<=k-1;m++){
                cospsi=(sin(thetapi[m]*(pi/180.0))*sin(thetapi[m+1]*(pi/180.0)))+(cos(thetapi[m]*(pi/180.0))*cos(thetapi[m+1]*(pi/180.0)))*cos((phipi[m+1]-phipi[m])*(pi/180.0));
                if(m<k-1){
                    
                    deltas=deltas+sqrt(fabs(pow(rpi[m]/rearth,2)+pow(rpi[m+1]/rearth,2)-(2.0*((rpi[m]*rpi[m+1])/(rearth*rearth))*cospsi)));
                }
                else{
                    deltas=deltas+sqrt(fabs(pow(rpi[m]/rearth,2)+pow(rmodel[index[3]]/rearth,2)-(2.0*((rpi[m]*rmodel[index[3]])/(rearth*rearth))*cospsi)));
                }
            }
            
        }
        
        if(deltas!=deltas || int1!=int1){
            printf("The routine has found a Not a Number (NaN) during the calculation\n");
            fclose(outputindex);
            outfile=fopen("raypathoutmodel.txt","a+");
            outfile2=fopen("rayout.txt","w");
            fprintf(outfile,"%d\n",anor);
            fprintf(outfile2,"%d\n",anor);
            fclose(outfile);
            fclose(outfile2);
            printf("It's possible that the index file is wrong.\n");
            printf("Delete it and try to find the problem!\n");
            printf("Check the raypathoutmodel.txt file to see the raypath associated with this problem!\n");
            printf("I have not experienced this problem during my calculation.\n");
            printf("Thus, I'd appreciate if you reported this problem to me.\n");
            printf("Terminating the program...\n");
            free(rp);
            rp=NULL;
            free(thetap);
            thetap=NULL;
            free(phip);
            phip=NULL;
            free(x1);
            x1=NULL;
            free(x2);
            x2=NULL;
            free(rpi);
            rpi=NULL;
            free(thetapi);
            thetapi=NULL;
            free(phipi);
            phipi=NULL;
            free(rmodel);
            rmodel=NULL;
            free(latmodel);
            latmodel=NULL;
            free(lonmodel);
            lonmodel=NULL;
            free(index);
            index=NULL;
            exit(EXIT_FAILURE);
        }
        tcalc=tcalc+((int1*deltas*rearth))*coordvnew[1];
        i=i+1;
        j=k+1;
        k=coordp[i]-1;
    }
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    i=1;
    j=1;
    k=coordp[1]-1;
    outputindex=fopen(fileicrayspath,"w");
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    while (i<=(sizecoordp-1)){
        rnew[1]=rpi[j];
        thetanew[1]=thetapi[j];
        phinew[1]=phipi[j];
        rnew[2]=rpi[k];
        thetanew[2]=thetapi[k];
        phinew[2]=phipi[k];
        velblock(nr,nt,np,rnew,thetanew,phinew,latmin,lonmin,dt,dp,coordvnew);
        index[1]=(int)coordvnew[2];
        index[2]=(int)coordvnew[3];
        index[3]=(int)coordvnew[4];
        int1=1.0/pow(coordvnew[1],2);
        deltas=0.0;
        if(i<(sizecoordp-1)){
            for(m=j;m<=k-1;m++){
                cospsi=(sin(thetapi[m]*(pi/180.0))*sin(thetapi[m+1]*(pi/180.0)))+(cos(thetapi[m]*(pi/180.0))*cos(thetapi[m+1]*(pi/180.0)))*cos((phipi[m+1]-phipi[m])*(pi/180.0));
                deltas=deltas+sqrt(fabs(pow(rpi[m]/rearth,2)+pow(rpi[m+1]/rearth,2)-(2.0*((rpi[m]*rpi[m+1])/(rearth*rearth))*cospsi)));
            }
        }
        else{
            for(m=j;m<=k-1;m++){
                cospsi=(sin(thetapi[m]*(pi/180.0))*sin(thetapi[m+1]*(pi/180.0)))+(cos(thetapi[m]*(pi/180.0))*cos(thetapi[m+1]*(pi/180.0)))*cos((phipi[m+1]-phipi[m])*(pi/180.0));
                if(m<k-1){
                    
                    deltas=deltas+sqrt(fabs(pow(rpi[m]/rearth,2)+pow(rpi[m+1]/rearth,2)-(2.0*((rpi[m]*rpi[m+1])/(rearth*rearth))*cospsi)));
                }
                else{
                    deltas=deltas+sqrt(fabs(pow(rpi[m]/rearth,2)+pow(rmodel[index[3]]/rearth,2)-(2.0*((rpi[m]*rmodel[index[3]])/(rearth*rearth))*cospsi)));
                }
            }
            
        }
        if (ttaup !=0){
            deltas=deltas+(((ttaup-tcalc)/(ttaup))*deltas); //path length correction
        }
        if(deltas!=deltas || int1!=int1){
            printf("The routine has found a Not a Number (NaN) during the calculation\n");
            fclose(outputindex);
            outfile=fopen("raypathoutmodel.txt","a+");
            outfile2=fopen("rayout.txt","w");
            fprintf(outfile,"%d\n",anor);
            fprintf(outfile2,"%d\n",anor);
            fclose(outfile);
            fclose(outfile2);
            printf("It's possible that the index file is wrong.\n");
            printf("Delete it and try to find the problem!\n");
            printf("Check the raypathoutmodel.txt file to see the raypath associated with this problem!\n");
            printf("I have not experienced this problem during my calculation.\n");
            printf("Thus, I'd appreciate if you reported this problem to me.\n");
            printf("Terminating the program...\n");
            free(rp);
            rp=NULL;
            free(thetap);
            thetap=NULL;
            free(phip);
            phip=NULL;
            free(x1);
            x1=NULL;
            free(x2);
            x2=NULL;
            free(rpi);
            rpi=NULL;
            free(thetapi);
            thetapi=NULL;
            free(phipi);
            phipi=NULL;
            free(rmodel);
            rmodel=NULL;
            free(latmodel);
            latmodel=NULL;
            free(lonmodel);
            lonmodel=NULL;
            free(index);
            index=NULL;
            exit(EXIT_FAILURE);
        }
        fprintf(outputindex,"%d %d %d %lf\n",index[1],index[2],index[3],((int1*deltas*rearth))); //Trapezoidal rule
        i=i+1;
        j=k+1;
        k=coordp[i]-1;
    }
    fclose(outputindex);
    free(rp);
    rp=NULL;
    free(thetap);
    thetap=NULL;
    free(phip);
    phip=NULL;
    free(x1);
    x1=NULL;
    free(x2);
    x2=NULL;
    free(rpi);
    rpi=NULL;
    free(thetapi);
    thetapi=NULL;
    free(phipi);
    phipi=NULL;
    free(rmodel);
    rmodel=NULL;
    free(latmodel);
    latmodel=NULL;
    free(lonmodel);
    lonmodel=NULL;
    free(index);
    index=NULL;
    return 0;
    
}


