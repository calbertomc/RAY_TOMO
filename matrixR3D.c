//*********************************************************************************************************************************************
//  matrixR3D.c
//  Routine developed
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

#include <stdio.h>
#include <stdlib.h>

int main(){
    
    int
    nz,
    ny,
    nx,
    np,
    i,
    j,
    k,
    m,
    t,
    v,
    *pos,
    **R;

    FILE *inputparameters;

    printf("Openning the file ipmatrixR3D.txt...\n");
    inputparameters=fopen("ipmatrixR3D.txt","r");
    if(inputparameters==NULL){
        printf("File does not exist or has been moved: ipmatrixR3D.txt\n");
        printf("Terminating the program...\n");
        printf("Try it again!!!\n");
        fclose(inputparameters);
        exit(EXIT_FAILURE);
    }
    printf("Done!\n");
    printf("Reading the file ipmatrixR3D.txt...\n");
    fscanf(inputparameters,"%d %d %d",&nz,&ny,&nx);
    fclose(inputparameters);
    printf("Done!\n");

    
    
    if(nx>=2 && ny>=2 && nz>=2){
        np=nx*ny*nz;
        k=1;
        t=1;
        
        printf("Allocating memory...\n");
        pos=(int *)malloc((5*np) * sizeof(int));
        R=(int **)malloc((np+1) * sizeof(int*));
        printf("Done!\n");
        for(i=1;i<=np;i++){
            R[i]=(int *)malloc((np+1) * sizeof(int));
        }
//        for(i=1;i<=np;i++){
//            for(j=1;j<=np;j++){
//                R[i][j]=0;
//            }
//        }
        printf("Filling matrix\n");
        for(m=1;m<=nz;m++){
            for(i=1;i<=ny;i++){
                for(j=1;j<=nx;j++){
                    if(j==1 && i==1 && m==1){
                        R[t][t]=3;
                        pos[k]=2;
                        k++;
                        pos[k]=nx+1;
                        k++;
                        pos[k]=(nx*ny)+1;
                        k++;
                        t++;
                    }
                    else if(j==nx && i==1 && m==1){
                        R[t][t]=3;
                        pos[k]=nx-1;
                        k++;
                        pos[k]=2*nx;
                        k++;
                        pos[k]=(nx*ny)+nx;
                        k++;
                        t++;
                    }
                    else if(j==1 && i==ny && m==1){
                        R[t][t]=3;
                        pos[k]=nx*(ny-2)+1;
                        k++;
                        pos[k]=nx*(ny-1)+2;
                        k++;
                        pos[k]=2*(nx*ny)-nx+1;
                        k++;
                        t++;
                    }
                    else if(j==nx && i==ny && m==1){
                        R[t][t]=3;
                        pos[k]=nx*(ny-1);
                        k++;
                        pos[k]=(nx*ny)-1;
                        k++;
                        pos[k]=2*(nx*ny);
                        k++;
                        t++;
                    }
                    else if(j==1 && i==1 && m==nz){
                        R[t][t]=3;
                        pos[k]=(nx*ny*(nz-1))+2;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+nx+1;
                        k++;
                        pos[k]=(nx*ny*(nz-2))+1;
                        k++;
                        t++;
                    }
                    else if(j==nx && i==1 && m==nz){
                        R[t][t]=3;
                        pos[k]=(nx*ny*(nz-1))+nx-1;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+2*nx;
                        k++;
                        pos[k]=(nx*ny*(nz-2))+nx;
                        k++;
                        t++;
                    }
                    else if(j==1 && i==ny && m==nz){
                        R[t][t]=3;
                        pos[k]=(nx*ny*(nz-1))+nx*(ny-2)+1;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+nx*(ny-1)+2;
                        k++;
                        pos[k]=(nx*ny*(nz-2))+nx*(ny-1)+1;
                        k++;
                        t++;
                    }
                    else if(j==nx && i==ny && m==nz){
                        R[t][t]=3;
                        pos[k]=(nx*ny*(nz-1))+nx*(ny-1);
                        k++;
                        pos[k]=(nx*ny*(nz-1))+(nx*ny)-1;
                        k++;
                        pos[k]=(nx*ny*(nz-1));
                        k++;
                        t++;
                    }
                    else if(j==1 && i==1 && m>1 && m<nz){
                        R[t][t]=4;
                        pos[k]=(nx*ny*(m-1))+2;
                        k++;
                        pos[k]=(nx*ny*(m-1))+nx+1;
                        k++;
                        pos[k]=(nx*ny*(m-1))+(nx*ny)+1;
                        k++;
                        pos[k]=(nx*ny*(m-2))+1;
                        k++;
                        t++;
                    }
                    else if(j==nx && i==1 && m>1 && m<nz){
                        R[t][t]=4;
                        pos[k]=(nx*ny*(m-1))+nx-1;
                        k++;
                        pos[k]=(nx*ny*(m-1))+2*nx;
                        k++;
                        pos[k]=(nx*ny*m)+nx;
                        k++;
                        pos[k]=(nx*ny*(m-2))+nx;
                        k++;
                        t++;
                    }
                    else if(j==1 && i==ny && m>1 && m<nz){
                        R[t][t]=4;
                        pos[k]=(nx*ny*(m-1))+nx*(ny-2)+1;
                        k++;
                        pos[k]=(nx*ny*(m-1))+nx*(ny-1)+2;
                        k++;
                        pos[k]=(nx*ny*m)+(nx*ny)-nx+1;
                        k++;
                        pos[k]=(nx*ny*(m-2))+nx*(ny-1)+1;
                        k++;
                        t++;
                    }
                    else if(j==nx && i==ny && m>1 && m<nz){
                        R[t][t]=4;
                        pos[k]=(nx*ny*(m-1))+nx*(ny-1);
                        k++;
                        pos[k]=(nx*ny*(m-1))+(nx*ny)-1;
                        k++;
                        pos[k]=(nx*ny*(m-1));
                        k++;
                        pos[k]=(nx*ny*(m+1));
                        k++;
                        t++;
                    }
                    else if(j==1 && i>1 && i<ny && m==1){
                        R[t][t]=4;
                        pos[k]=nx*(i-2)+1;
                        k++;
                        pos[k]=nx*(i-1)+2;
                        k++;
                        pos[k]=nx*i+1;
                        k++;
                        pos[k]=(nx*ny)+nx*(i-1)+1;
                        k++;
                        t++;
                    }
                    else if(j==nx && i>1 && i<ny && m==1){
                        R[t][t]=4;
                        pos[k]=nx*(i-1);
                        k++;
                        pos[k]=(nx*i)-1;
                        k++;
                        pos[k]=nx*(i+1);
                        k++;
                        pos[k]=(nx*ny)+nx*i;
                        k++;
                        t++;
                    }
                    else if(i==1 && j>1 && j<nx && m==1){
                        R[t][t]=4;
                        pos[k]=j-1;
                        k++;
                        pos[k]=nx+j;
                        k++;
                        pos[k]=j+1;
                        k++;
                        pos[k]=(nx*ny)+j;
                        k++;
                        t++;
                    }
                    else if(i==ny && j>1 && j<nx && m==1){
                        R[t][t]=4;
                        pos[k]=nx*(ny-1)+j-1;
                        k++;
                        pos[k]=nx*(ny-2)+j;
                        k++;
                        pos[k]=nx*(ny-1)+j+1;
                        k++;
                        pos[k]=(nx*ny)+(nx*(ny-1))+j;
                        k++;
                        t++;
                    }
                    else if(j==1 && i>1 && i<ny && m==nz){
                        R[t][t]=4;
                        pos[k]=(nx*ny*(nz-1))+nx*(i-2)+1;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+nx*(i-1)+2;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+nx*i+1;
                        k++;
                        pos[k]=(nx*ny*(nz-2))+nx*(i-1)+1;
                        k++;
                        t++;
                    }
                    else if(j==nx && i>1 && i<ny && m==nz){
                        R[t][t]=4;
                        pos[k]=(nx*ny*(nz-1))+nx*(i-1);
                        k++;
                        pos[k]=(nx*ny*(nz-1))+(nx*i)-1;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+nx*(i+1);
                        k++;
                        pos[k]=(nx*ny*(nz-2))+nx*i;
                        k++;
                        t++;
                    }
                    else if(i==1 && j>1 && j<nx && m==nz){
                        R[t][t]=4;
                        pos[k]=(nx*ny*(nz-1))+j-1;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+nx+j;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+j+1;
                        k++;
                        pos[k]=(nx*ny*(nz-2))+j;
                        k++;
                        t++;
                    }
                    else if(i==ny && j>1 && j<nx && m==nz){
                        R[t][t]=4;
                        pos[k]=(nx*ny*(nz-1))+nx*(ny-1)+j-1;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+nx*(ny-2)+j;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+nx*(ny-1)+j+1;
                        k++;
                        pos[k]=(nx*ny*(nz-2))+(nx*(ny-1))+j;
                        k++;
                        t++;
                    }
                    else if(j==1 && i>1 && i<ny && m>1 && m<nz){
                        R[t][t]=5;
                        pos[k]=(nx*ny*(m-1))+nx*(i-2)+1;
                        k++;
                        pos[k]=(nx*ny*(m-1))+nx*(i-1)+2;
                        k++;
                        pos[k]=(nx*ny*(m-1))+nx*i+1;
                        k++;
                        pos[k]=(nx*ny*(m-2))+nx*(i-1)+1;
                        k++;
                        pos[k]=(nx*ny*m)+nx*(i-1)+1;
                        k++;
                        t++;
                    }
                    else if(j==nx && i>1 && i<ny && m>1 && m<nz){
                        R[t][t]=5;
                        pos[k]=(nx*ny*(m-1))+nx*(i-1);
                        k++;
                        pos[k]=(nx*ny*(m-1))+(nx*i)-1;
                        k++;
                        pos[k]=(nx*ny*m);
                        k++;
                        pos[k]=(nx*ny*(m-2))+nx*i;
                        k++;
                        pos[k]=(nx*ny*m)+nx*i;
                        k++;
                        t++;
                    }
                    else if(i==1 && j>1 && j<nx && m>1 && m<nz){
                        R[t][t]=5;
                        pos[k]=(nx*ny*(m-1))+j-1;
                        k++;
                        pos[k]=(nx*ny*(m-1))+nx+j;
                        k++;
                        pos[k]=(nx*ny*(m-1))+j+1;
                        k++;
                        pos[k]=(nx*ny*(m-2))+j;
                        k++;
                        pos[k]=(nx*ny*m)+j;
                        k++;
                        t++;
                    }
                    else if(i==ny && j>1 && j<nx && m>1 && m<nz){
                        R[t][t]=5;
                        pos[k]=(nx*ny*(m-1))+nx*(ny-1)+j-1;
                        k++;
                        pos[k]=(nx*ny*(m-1))+nx*(ny-2)+j;
                        k++;
                        pos[k]=(nx*ny*(m-1))+nx*(ny-1)+j+1;
                        k++;
                        pos[k]=(nx*ny*(m-1))+(nx*ny)+(nx*(ny-1))+j;
                        k++;
                        pos[k]=(nx*ny*(m-2))+nx*(ny-1)+j;
                        k++;
                        t++;
                    }
                    else if(i>1 && i<ny && j>1 && j<nx && m==1){
                        R[t][t]=5;
                        pos[k]=j+(i-2)*nx;
                        k++;
                        pos[k]=i*nx+j;
                        k++;
                        pos[k]=(i-1)*nx+(j-1);
                        k++;
                        pos[k]=(i-1)*nx+(j+1);
                        k++;
                        pos[k]=(nx*ny*m)+(i-1)*nx+j;
                        k++;
                        t++;
                    }
                    else if(i>1 && i<ny && j>1 && j<nx && m==nz){
                        R[t][t]=5;
                        pos[k]=(nx*ny*(nz-1))+j+(i-2)*nx;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+i*nx+j;
                        k++;
                        pos[k]=(nx*ny*(nz-1))+(i-1)*nx+(j-1);
                        k++;
                        pos[k]=(nx*ny*(nz-1))+(i-1)*nx+(j+1);
                        k++;
                        pos[k]=(nx*ny*(m-2))+(i-1)*nx+j;
                        k++;
                        t++;
                    }
                    else{
                        R[t][t]=6;
                        pos[k]=(nx*ny*(m-1))+j+(i-2)*nx;
                        k++;
                        pos[k]=(nx*ny*(m-1))+i*nx+j;
                        k++;
                        pos[k]=(nx*ny*(m-1))+(i-1)*nx+(j-1);
                        k++;
                        pos[k]=(nx*ny*(m-1))+(i-1)*nx+(j+1);
                        k++;
                        pos[k]=(nx*ny*(m-2))+j+(i-1)*nx;
                        k++;
                        pos[k]=(nx*ny*m)+(i-1)*nx+j;
                        k++;
                        t++;
                    }
                }
            }
        }
        k=1;
        for(i=1;i<=np;i++){
            v=R[i][i];
            for(j=1;j<=v;j++){
                R[pos[k]][i]=-1;
                k++;
            }
        }
        //Ws=R'*R;
        //save Wsmatrix.dat Ws -ascii
//        for(i=1;i<=np;i++){
//            for(j=1;j<=np;j++){
//                printf("%d ",R[i][j]);
//            }
//            printf("\n");
//        }
        free(R);
        R==NULL;
        free(pos);
        pos=NULL;
    }
    else{
        printf("This algorithm only works if nx>=2, ny>=2 and nz>=2!\n");
        printf("Please, run it again and modify nx and ny!\n");
    }
    
    return 0;
}
