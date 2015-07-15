#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[]) {
  FILE *fr;
  int Na,Nx,Ny,Nz,N3;
  double *data=NULL;
  double *out=NULL;
  double sum;
  double lx,ly,lz,L3,V;
  int i,j,k,s;
  char *string=NULL;
  size_t len=0;
  ssize_t read;

  string = (char *) malloc(80 * sizeof(char));

  if(argc != 2){
    printf("failed, need 1 arg\n");
    return;
  }
  fr = fopen(argv[1],"r");

  for(i=0;i<2;i++){
    read = getline(&string, &len, fr);
  }

  fscanf(fr,"%d",&Na);
  getline(&string, &len, fr);

  fscanf(fr,"%d",&Nx);
  /* read in voxel volume in x */
  fscanf(fr,"%lf",&lx);
  getline(&string, &len, fr);
  fscanf(fr,"%d",&Ny);
  /* read in voxel volume in y */
  fscanf(fr,"%lf",&ly);
  fscanf(fr,"%lf",&ly);
  getline(&string, &len, fr);
  fscanf(fr,"%d",&Nz);
  /* read in voxel volume in z */
  fscanf(fr,"%lf",&lz);
  fscanf(fr,"%lf",&lz);
  fscanf(fr,"%lf",&lz);
  getline(&string, &len, fr);

  L3= lx * ly * lz;
  N3= Nx * Ny * Nz;
  V = L3 * N3;
  data = (double *) malloc(N3 * sizeof(double));
  out = (double *) malloc(Nx * sizeof(double));

  for(i=0;i<Na;i++) getline(&string,&len,fr);
  
  //printf("%lf %lf %lf\n",lx,ly,lz);
  for(i=0;i<N3;i++){
    if(!feof(fr)){
       fscanf(fr,"%le",&data[i]);
      //printf("%le\n",data[i]);
    }else{
      printf("end of file at i=%d\n",i-1);
      break;
    }
  }


  //printf("%d %d %d\n",Nx, Ny, Nz);

  for(i=0;i<Nx;i++){
    sum=0;
    for(j=0;j<Ny;j++){
      for(k=0;k<Nz;k++){
        s = k + Nz*j + Nz*Ny*i;    
        sum=sum+(data[s]*L3)/V;
      }
    }
    out[k]=sum;
    printf("%e\n",out[k]);
  }

/*
%  for(i=0;i<N3;i++){
%    printf("%e\n",data[i]);
%  }
*/
  fclose(fr);
  free(string);
  free(data);
  free(out);
  data=NULL;
  return 0;
}
