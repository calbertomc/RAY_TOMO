//*********************************************************************************************************************************************
//  matrixL.c
//  Routine developed
//  Created by Carlos Alberto Chaves on 09/05/2015.
//  University of Sao Paulo - University of Michigan
//  carlos.chaves@iag.usp.br; cchaves@umich.edu; calbertochaves@gmail.com (main)
//  Version 1.0
//*********************************************************************************************************************************************
//Sorry, but I usually don't comment my codes. I don't have much patience neither time to do it!


#include <stdio.h>
#include <stdlib.h>

int main(){
    
    int
    nz,
    ny,
    nx,
    i,
    j,
    m,
    mhypo,
    t;

    FILE *inputparameters;
    FILE *outmatrix;

    printf("Openning the file ipmatrixL.txt...\n");
    inputparameters=fopen("ipmatrixL.txt","r");
    if(inputparameters==NULL){
        printf("File does not exist or has been moved: ipmatrixL.txt\n");
        printf("Terminating the program...\n");
        printf("Try it again!!!\n");
        fclose(inputparameters);
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    printf("Reading the file ipmatrixL.txt...\n");
    fscanf(inputparameters,"%d %d %d %d",&nz,&ny,&nx,&mhypo);
    fclose(inputparameters);
    printf("Done!\n");

    
    
    if(nx>=2 && ny>=2 && nz>=2){
        t=1;
        outmatrix=fopen("matrixL.txt","w");
        printf("Filling matrix\n");
        for(m=1;m<=nz;m++){
            for(i=1;i<=ny;i++){
                for(j=1;j<=nx;j++){
                    if(j==1 && i==1 && m==1){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,3);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+2-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+nx+1-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny)+1)-1,-1);
                        t++;
                    }
                    else if(j==nx && i==1 && m==1){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,3);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+nx-1-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(2*nx)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny)+nx)-1,-1);
                        t++;
                    }
                    else if(j==1 && i==ny && m==1){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,3);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(ny-2)+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(ny-1)+2)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(2*(nx*ny)-nx+1)-1,-1);
                        t++;
                    }
                    else if(j==nx && i==ny && m==1){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,3);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(ny-1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny)-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(2*(nx*ny))-1,-1);
                        t++;
                    }
                    else if(j==1 && i==1 && m==nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,3);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+2)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-2))+1)-1,-1);
                        t++;
                    }
                    else if(j==nx && i==1 && m==nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,3);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+2*nx)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-2))+nx)-1,-1);
                        t++;
                    }
                    else if(j==1 && i==ny && m==nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,3);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(ny-2)+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(ny-1)+2)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-2))+nx*(ny-1)+1)-1,-1);
                        t++;
                    }
                    else if(j==nx && i==ny && m==nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,3);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(ny-1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+(nx*ny)-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1)))-1,-1);
                        t++;
                    }
                    else if(j==1 && i==1 && m>1 && m<nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+2)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+(nx*ny)+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-2))+1)-1,-1);
                        t++;
                    }
                    else if(j==nx && i==1 && m>1 && m<nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+2*nx)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*m)+nx)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-2))+nx)-1,-1);
                        t++;
                    }
                    else if(j==1 && i==ny && m>1 && m<nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*(ny-2)+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*(ny-1)+2)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*m)+(nx*ny)-nx+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-2))+nx*(ny-1)+1)-1,-1);
                        t++;
                    }
                    else if(j==nx && i==ny && m>1 && m<nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*(ny-1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+(nx*ny)-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1)))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m+1)))-1,-1);
                        t++;
                    }
                    else if(j==1 && i>1 && i<ny && m==1){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(i-2)+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(i-1)+2)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*i+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny)+nx*(i-1)+1)-1,-1);
                        t++;
                    }
                    else if(j==nx && i>1 && i<ny && m==1){
                        fprintf(outmatrix,"%d %d %d\n",t,mhypo+t,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(i-1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*i)-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(i+1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny)+nx*i)-1,-1);
                        t++;
                    }
                    else if(i==1 && j>1 && j<nx && m==1){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(j-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(j+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny)+j)-1,-1);
                        t++;
                    }
                    else if(i==ny && j>1 && j<nx && m==1){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(ny-1)+j-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(ny-2)+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(nx*(ny-1)+j+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny)+(nx*(ny-1))+j)-1,-1);
                        t++;
                    }
                    else if(j==1 && i>1 && i<ny && m==nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(i-2)+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(i-1)+2)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*i+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-2))+nx*(i-1)+1)-1,-1);
                        t++;
                    }
                    else if(j==nx && i>1 && i<ny && m==nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(i-1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+(nx*i)-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(i+1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-2))+nx*i)-1,-1);
                        t++;
                    }
                    else if(i==1 && j>1 && j<nx && m==nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+j-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+j+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-2))+j)-1,-1);
                        t++;
                    }
                    else if(i==ny && j>1 && j<nx && m==nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,4);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(ny-1)+j-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(ny-2)+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+nx*(ny-1)+j+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-2))+(nx*(ny-1))+j)-1,-1);
                        t++;
                    }
                    else if(j==1 && i>1 && i<ny && m>1 && m<nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,5);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*(i-2)+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*(i-1)+2)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*i+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-2))+nx*(i-1)+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*m)+nx*(i-1)+1)-1,-1);
                        t++;
                    }
                    else if(j==nx && i>1 && i<ny && m>1 && m<nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,5);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*(i-1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+(nx*i)-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*m))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-2))+nx*i)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*m)+nx*i)-1,-1);
                        t++;
                    }
                    else if(i==1 && j>1 && j<nx && m>1 && m<nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,5);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+j-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+j+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-2))+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*m)+j)-1,-1);
                        t++;
                    }
                    else if(i==ny && j>1 && j<nx && m>1 && m<nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,5);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*(ny-1)+j-1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*(ny-2)+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+nx*(ny-1)+j+1)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+(nx*ny)+(nx*(ny-1))+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-2))+nx*(ny-1)+j)-1,-1);
                        t++;
                    }
                    else if(i>1 && i<ny && j>1 && j<nx && m==1){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,5);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(j+(i-2)*nx)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+(i*nx+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((i-1)*nx+(j-1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((i-1)*nx+(j+1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*m)+(i-1)*nx+j)-1,-1);
                        t++;
                    }
                    else if(i>1 && i<ny && j>1 && j<nx && m==nz){
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,5);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+j+(i-2)*nx)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+i*nx+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+(i-1)*nx+(j-1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(nz-1))+(i-1)*nx+(j+1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-2))+(i-1)*nx+j)-1,-1);
                        t++;
                    }
                    else{
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+t-1,6);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+j+(i-2)*nx)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+i*nx+j)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+(i-1)*nx+(j-1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-1))+(i-1)*nx+(j+1))-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*(m-2))+j+(i-1)*nx)-1,-1);
                        fprintf(outmatrix,"%d %d %d\n",mhypo+t-1,mhypo+((nx*ny*m)+(i-1)*nx+j)-1,-1);
                        t++;
                    }
                }
            }
        }
        fclose(outmatrix);
        printf("Done!\n");
    }
    else{
        printf("This algorithm only works if nx>=2, ny>=2 and nz>=2!\n");
        printf("Please, run it again and modify nx and ny!\n");
    }
    
    return 0;
}

