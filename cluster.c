#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "C:\Users\priya\Documents\Important stuff\6th sem\Sem - 6\try\MD\1-rn.c"

#define l0 64
#define L (l0*R0)
#define npart 4096
#define R0 1.05
#define sqrtn pow(npart,0.5)
#define tmax 10000
#define densi pow(sqrtn/(1.05*(sqrtn + 1)),2)
#define temp 0.5
#define v0 0.9
#define box pow(npart/densi,0.5)
#define dt 0.01
#define Rmin 1.0
#define Rmax 1.1


int neigh1[npart],neigh2[npart][4],neigh3[npart][4];
int NNcrack[npart][4];
int lattice[npart];
int ncl,t;
double vir,vir_sum,P,tmp;
double fx[npart],fy[npart];
double T0,rho,V,hL;
double x[npart],y[npart];
double xi[npart],yi[npart];
double xcm,ycm,xcm0,ycm0;
double xm[npart],ym[npart];
double vx[npart],vy[npart];
int find_crack();
int burn_nn();
void starting_position();
void create_neighbour();
void update_nearest_neighbors();
double force();
double integrate();
int burn_nn();
void crackspan();
double print_spanningcrack();
void printlatindex();
void CPC();


int main(void){

	clock_t t1,t2;
	t1 = clock();

	// creating the pointer variable to the files
	FILE *file1, *file2, *file3, *file4, *file5, *fp1, *fp2;

	file1 = fopen("posi.dat", "w");	//opening the files
	file2 = fopen("blender.xyz", "w");
	file3 = fopen("posf.dat", "w");
	file4 = fopen("Energy_1.dat", "w");
	file5 = fopen("withbond.dat", "w");
	fp1=fopen("FILE1.dat","w");
	fp2 = fopen("FILE2.dat","w");


	// variables
	double x[npart], y[npart], vx[npart], vy[npart], fx[npart], fy[npart], vvx[npart], vvy[npart];
	
    // allocate memory for other arrays
	double Ep=0.0,Ec=0.0;
	for(int i=0;i<npart;i++){
	x[i]=0; y[i]=0; vx[i]=0; vy[i]=0; fx[i]=0; fy[i]=0; vvx[i]=0; vvy[i]=0;  // zeroing of vectors
		}



	// initiating position
	starting_position();
	create_neighbour();

	fprintf(file2,"%d \n",npart);
    fprintf(file2,"Particle Center \n");
	for(int j=0;j<npart;j++){
		fprintf(file2,"Default      %lf      %lf      0 \n",xi[j],yi[j]);	// blender
		fprintf(file1,"%lf      %lf\n",xi[j],yi[j]);	// initial position
		}
		
	int crackflag = 0;
	double sumv2,sumvx,sumvy;
	int i,j,ibond=0;
	for(int t=1;t<=tmax;t++){
		force();
		integrate();		
		update_nearest_neighbors();
	if(crackflag==0){
    crackflag=find_crack(t);
    if(crackflag==1){
	int ncl=burn_nn();
	printf("%e %d\n",t*dt,burn_nn());
		//print_state_withbond_Lammpsdata(itime);
	printf("Crack found\n");
	for(int j=0;j<npart;j++){
		fprintf(file3,"%lf      %lf\n",xi[j],yi[j]);	// pos final
	}
	crackspan(t,ncl,NNcrack);
	for(i=0;i<npart;i++)
    for(j=0;j<4;j++)
		if(NNcrack[i][j]>0)
			fprintf(fp2,"%d %d %d\n",i*4+j+1,NNcrack[i][j],i+1);
			ibond++;
	print_spanningcrack(t,ncl,NNcrack);

	printf("%e crack = %d\n",t*dt,find_crack(t));
  fclose(fp1);
  fclose(fp2);
return 0;
      }
    }

		
		vir_sum+=vir;
		Ep = force();
		Ec = integrate();

		
		fprintf(file4,"%d      %lf      %lf \n",t,Ep,Ec);  // energy

	}


	// print_spanningcrack(t,ncl,NNcrack);
	// double temp2 = print_spanningcrack(t,ncl,NNcrack);

 
	t2 = clock();
	printf("%f",(((double)t2 - (double)t1) / 1000000.0F ));
	return 0;
}


//square grid
void starting_position(){
	int i,j,l3;
  	double ix,iy,iz;
  	double sumvx,sumvy,sumvz,sumv2,fs;  

    
    for(i=0;i<npart;i++){
      x[i]=R0*(i/64)+R0/2.;
      y[i]=R0*(i%64)+R0/2.;      
      xi[i]=x[i];    yi[i]=y[i];    //storing x(t=0)
      xcm0+=xi[i];    ycm0+=yi[i];    //cm for t=0;
    }
    
    xcm0/=npart;  ycm0/=npart; //Center of mass for the particles
    
    sumv2=0;sumvx=0;sumvy=0;sumvz=0;
	//uniform velocity distribution
    for (i=0;i<npart;i++){
      vx[i]=(ran2(&iseed)-0.5)*10;
      vy[i]=(ran2(&iseed)-0.5)*10;
    }    
 }


void create_neighbour(){
  int i,j;
  double xr,yr,zr,r,r2;
  for(i=0;i<npart-1;i++){
    for(j=i+1;j<npart;j++){
      xr=xi[i]-xi[j]; 
      yr=yi[i]-yi[j]; 
      
      
      	if (xr>=0.5*box){ xr=xr-box; }
    	if (xr<=-0.5*box){ xr=xr+box; }
        if (yr>=0.5*box){ yr=yr-box; }
        if (yr<=-0.5*box){ yr=yr+box; }

      
      r2=xr*xr+yr*yr;
      r=sqrt(r2);
      if(r<R0+0.01){
	neigh2[i][neigh1[i]]=j;
	neigh2[j][neigh1[j]]=i;

	neigh3[i][neigh1[i]]=1;
	neigh3[j][neigh1[j]]=1;

	neigh1[i]++;
	neigh1[j]++;
      }
      
    }
  }
  
}

void update_nearest_neighbors(){
  
  	int i,j,k;
  	double xr,yr,zr,r,r2;
  
  	for(i=0;i<npart;i++){
    for(k=0;k<neigh1[i];k++){

      if(neigh3[i][k]==1){
		j=neigh2[i][k];
		xr=xi[i]-xi[j]; 
		yr=yi[i]-yi[j]; 
	
		if (xr>box/2)       xr-=box;
		else if (xr<-box/2) xr+=box;
		if (yr>box/2)       yr-=box;
		else if (yr<-box/2) yr+=box;

		r2=xr*xr+yr*yr;
		r=sqrt(r2);
	
		if(r>1.1)neigh3[i][k]=0;
      }
    }
  }

}

// interaction force between particles
double force(){
    double ene=0.0;
	int i,j,k;
	double xr,yr,zr,r,r2,ff = 0,Ep;
    for(i=0;i<npart;i++){
        fx[i]=fy[i]=0.0;
    }
    for(i=0;i<npart;i++){
    	for(k=0;k<neigh1[k];k++){

      	if(neigh3[i][k]>0){

			j=neigh2[i][k];
	
			xr=xi[i]-xi[j]; 
			yr=yi[i]-yi[j]; 
	
		if (xr>box/2)       xr-=box;
		else if (xr<-box/2) xr+=box;
		if (yr>box/2)       yr-=box;
		else if (yr<-box/2) yr+=box;
	
		r2=xr*xr+yr*yr;
		r=sqrt(r2);

		if(r-R0<10e-5){ // scattering: exchnage velocity
	  	tmp=vx[i];vx[i]=vx[j];vx[j]=tmp;
	  	tmp=vy[i];vy[i]=vy[j];vy[j]=tmp;	
	}
	
	if(r>Rmin && r<Rmax){
	  double rR0=r-R0;
	  ff=-(100*rR0-50*rR0*rR0)/r;
	
	  fx[i]=fx[i]+ff*xr;
	  fx[j]=fx[j]-ff*xr;
	  
	  fy[i]=fy[i]+ff*yr;
	  fy[j]=fy[j]-ff*yr;
	  
	  double e=(100*0.5*rR0*rR0-50*rR0*rR0*rR0/3.0);
	  ene=ene+e;
	  vir+=ff*r2;
	}
	
    }
      
    }
  }
    Ep=ene/(double)npart;
    return Ep;
}

// integration of particle positions
double integrate(){
	for(int i=0;i<npart;i++){
		vx[i]+=0.5*dt*fx[i];
    	vy[i]+=0.5*dt*fy[i];
		xi[i]+=vx[i]*dt+0.5*fx[i]*dt*dt;
    	yi[i]+=vy[i]*dt+0.5*fy[i]*dt*dt;
	}

	
	force();

	

	for(int i=0;i<npart;i++){
		vx[i]+=0.5*dt*fx[i];
    	vy[i]+=0.5*dt*fy[i];
	}	

	double sumv2,sumvx,sumvy,Ec;
	for (int i=0;i<npart;i++){	     
      	sumvx+=vx[i];
      	sumvy+=vy[i];
      	sumv2+=(vx[i]*vx[i]+vy[i]*vy[i]);
    }
		Ec = 0.5*(sumv2/npart);
		return Ec;

}

// periodic boundary condition
void CPC(){
	for(int i=0;i<npart;i++){
		if (xi[i]>box/2)       xi[i]-=box;
		else if (xi[i]<-box/2) xi[i]+=box;
		if (yi[i]>box/2)       yi[i]-=box;
		else if (yi[i]<-box/2) yi[i]+=box;
	}
}

int  burn_nn()
{
  static int i,j,k,k1,k2,n,max=0;
  static int iroot,inow,rroot,sum=0;
  static double z;


  int *ip,*ns,*lattice;
  double *size;
  ip = malloc(npart* sizeof(int));
  ns = malloc(npart * sizeof(int));
  lattice = malloc(npart * sizeof(int));
  for(i=0;i<npart;i++)lattice[i]=0;

  size = malloc(npart * sizeof(double));
  for(i=0;i<npart;i++)size[i]=0;

  int ncl=0;
  for(i=0;i<npart;i++){//lattice scan for burning
    if(lattice[i]==0){
      ncl++;
      lattice[i]=ncl;
      n=1;
      ip[n]=i;
      k1=n;
      k2=k1;
    rr:
      for(j=k1;j<=k2;j++){
	iroot=ip[j];

	for(k=0;k<neigh1[iroot];k++){
	  if(neigh3[iroot][k]==1){

	    inow=neigh2[iroot][k];

	    if(lattice[inow]==0){
	      lattice[inow]=ncl;
	      n++;
	      ip[n]=inow;
	    }
	  }
	}
      }
      k1=k2+1;
      k2=n;
      if(k1<=k2)goto rr;
      if(n>max)max=n;
      ns[n]++;
      //printf("n=%d\n",n);
      for(k=1;k<=n;k++)
	ip[k]=0;

    size[ncl]=n;
    }//if(lattice[i]==0)
  }//end lattice scan

	FILE *fpc;
	FILE *fpp;
  char FILE1[50];
    int sum1=0,sum2=0,sum4=0,sum8=0,sum16=0,sum32=0,sum64=0,sum128=0,sum256=0,sum512 = 0;
	int sum_1,sum_2,sum_4,sum_8,sum_16,sum_32,sum_64,sum_128,sum_256,sum_512 = 0;
    int count1=0,count2=0,count4=0,count8=0,count16=0,count32=0,count64=0,count128=0,count256=0,count512 = 0;
	double rad1=0,rad2=0,rad4=0,rad8=0,rad16=0,rad32=0,rad64=0,rad128=0,rad256=0,rad512 = 0;

  fpc=fopen("FILE3.dat","w");
  fpp=fopen("FILE_10.dat","w");
  for(i=1;i<=ncl;i++){
	if(size[i]>0){
		fprintf(fpc,"%f\n",size[i]);
      if((size[i]) ==1)
	  {
		
		sum1 = sum1 + size[i];
		count1++;

	  }
	  else if (size[i] < 4 && size[i] >= 2 )
	  {
		sum2 = sum2 + size[i];
		count2++;
	  }
	  else if (size[i] < 8 && size[i] >= 4 )
	  {
		sum4 = sum4 + size[i];
		count4++;
	  }
	  else if (size[i] < 16 && size[i] >= 8 )
	  {
		sum8 = sum8 + size[i];
		count8++;
	  }
	  else if (size[i] < 32 && size[i] >= 16 )
	  {
		sum16 = sum16 + size[i];
		count16++;
	  }
	  else if (size[i] < 64 && size[i] >= 32 )
	  {
		sum32 = sum32 + size[i];
		count32++;
	  }
	  else if (size[i] < 128 && size[i] >= 64 )
	  {
		sum64 = sum64 + size[i];
		count64++;
	  }
	  else if (size[i] < 256 && size[i] >=128 )
	  {
		sum128 = sum128 + size[i];
		count128++;
	  }
	  else if (size[i] < 512 && size[i] >=128 )
	  {
		sum256 = sum256 + size[i];
		count256++;
	  }
	  
	  
  }}
//   sum_2 = sum2 - sum1;
//   sum_4 = sum4 - sum_2;
//   sum_8 = sum8 - sum_4;
//   sum_16= sum16 - sum_8;
//   sum_32 = sum32 - sum_16;
//   sum_64= sum64 - sum_32;
//   sum_128 = sum128 - sum_64;
//   sum_256 = sum256 - sum_128;

rad1 = pow(count1,0.5)*R0;
rad2 = sqrt(count2*(1.05)*(1.05));
rad4 = sqrt(count4*(1.05)*(1.05));
rad8 = sqrt(count8*(1.05)*(1.05));
rad16 = sqrt(count16*(1.05)*(1.05));
rad32 = sqrt(count32*(1.05)*(1.05));
rad64 = sqrt(count64*(1.05)*(1.05));
rad128 = sqrt(count128*(1.05)*(1.05));
rad256 = sqrt(count256*(1.05)*(1.05));


  fprintf(fpp,"1    %d    	%d 		%f\n",sum1,count1,rad1);
  fprintf(fpp,"2    %d 	  	%d 		%f\n",sum2,count2,rad2);
  fprintf(fpp,"4    %d 		%d 		%f\n",sum4,count4,rad4);
  fprintf(fpp,"8    %d 		%d 		%f\n",sum8,count8,rad8);
  fprintf(fpp,"16   %d 		%d 		%f\n",sum16,count16,rad16);
  fprintf(fpp,"32   %d 		%d 		%f\n",sum32,count32,rad32);
  fprintf(fpp,"64   %d 		%d 		%f\n",sum64,count64,rad64);
  fprintf(fpp,"128  %d 		%d 		%f\n",sum128,count128,rad128);
  fprintf(fpp,"256  %d 		%d 		%f\n",sum256,count256,rad256);


fclose(fpp);
  fclose(fpc);

  return ncl;
}

int  find_crack(int t)
{
  static int i,j,k,k1,k2,n,max=0; 
  static int iroot,inow,rroot,sum=0;
  static double z;
  int l1=64-1,ll1=64*l1,spannedcrack=0;

  int NNcrack[npart][4];
  for(i=0;i<npart;i++)
    for(j=0;j<4;j++)
      NNcrack[i][j]=0;
  

  int xcls[npart],ycls[npart];
  for(i=0;i<npart;i++){xcls[i]=0;ycls[i]=0;}
  

  int *ip,*ns,*lattice,*size;                               
  ip = malloc(npart * sizeof(int)); 
  ns = malloc(npart * sizeof(int));
  lattice = malloc(npart * sizeof(int));
  for(i=0;i<npart;i++)lattice[i]=0;

  size = malloc(npart * sizeof(int));
  for(i=0;i<npart;i++)size[i]=0;

  //Here I do bond percolation analysis of broken bond . If the span
  //of the cluster O(L): system spanning crack found and the function
  //return 1 else 0.
  for(i=0;i<npart;i++)lattice[i]=0;
    //{
    //  lattice[i]=-1;
    //  for(k=0;k<neigh1[i];k++){
    // 	if(neigh3[i][k]==0)lattice[i]=0;
    //  }
    //}
  
  int startcounting;
  int ncl=0;
  for(i=0;i<npart;i++){//lattice scan for burning

    startcounting=0;
    for(k=0;k<neigh1[i];k++)
      if(neigh3[i][k]==0)startcounting=1;
    

    if(lattice[i]==0 &&  startcounting==1){ 
      ncl++;     
      lattice[i]=ncl;
      n=1;       
      ip[n]=i;  
      k1=n;  
      k2=k1;

      xcls[i]=0;ycls[i]=0;

      //printf("ncl=%d\n",ncl);
      //printf("xstart=%d ystart=%d\n",xstart,ystart);
      
    rr:
      for(j=k1;j<=k2;j++){      
	iroot=ip[j];  
	
	for(k=0;k<neigh1[iroot];k++){
	  //for(k=0;k<4;k++){
	  //inow=iroot+iv[k];
	  //if(iroot%l0==l1 && k==0)inow=iroot-l1;
	  //if(iroot/l0==l1 && k==1)inow=iroot-ll1;
	  //if(iroot%l0==0 && k==2)inow=iroot+l1;
	  //if(iroot/l0==0 && k==3)inow=iroot+ll1;
	  
	  inow=neigh2[iroot][k];  
	    
	    if(lattice[inow]==0 && neigh3[iroot][k]==0){
	      lattice[inow]=ncl;
	      n++;
	      ip[n]=inow;

	      if(k==0){
		xcls[inow]=xcls[iroot]+1;
		ycls[inow]=ycls[iroot];
	      }
	      if(k==1){
		xcls[inow]=xcls[iroot];
		ycls[inow]=ycls[iroot]+1;
	      }
	      if(k==2){
		xcls[inow]=xcls[iroot]-1;
		ycls[inow]=ycls[iroot];
	      }
	      if(k==3){
		xcls[inow]=xcls[iroot];
		ycls[inow]=ycls[iroot]-1;
	      }

	    }//if lattice[inow]=0


	    
	    //storing crack information
	    if(lattice[inow]==ncl && neigh3[iroot][k]==0)
	      NNcrack[iroot][k]=ncl;
	    
	    
	}//loop over neighbour
      }
      k1=k2+1;
      k2=n;
      if(k1<=k2)goto rr;
      if(n>max)max=n;
      ns[n]++;
      
      for(k=1;k<=n;k++)
	ip[k]=0;

      size[ncl]=n;

      
      //printf("ncl=%d n=%d\n",ncl,n);
      
      int xmin=0,xmax=0,ymin=0,ymax=0;
      for(k=0;k<npart;k++){
	
	if(lattice[k]==ncl){
	  //printf("x[k]=%d y[k]=%d\n",xcls[k],ycls[k]);
	  if(xcls[k]<xmin)xmin=xcls[k];
	  if(xcls[k]>xmax)xmax=xcls[k];
	  
	  if(ycls[k]<ymin)ymin=ycls[k];
	  if(ycls[k]>ymax)ymax=ycls[k];
	}
      }
      
      
      //printf("xmin=%d ymin=%d\n",xmin,ymin);
      //printf("xmax=%d ymax=%d\n",xmax,ymax);
      //printf("xdiff=%d ydiff=%d\n",xmax-xmin,ymax-ymin); 
      
      
      if(xmax-xmin>=64 || ymax-ymin>=64)spannedcrack=1;
	

	  
	  
	  //printf("ncl=%d\n",ncl);
	  //printlatindex(lattice);
	  
	  
      
      

      
    }//if lattice[]=0
  }//end lattice scan

//   printlatindex(lattice);

  //printf("ncl=%d\n",ncl);

  if(spannedcrack)
    {
		  printlatindex(lattice);
	  crackspan(t,ncl,NNcrack);
	  
      
  }
  return spannedcrack;
}

void printlatindex(int lattice[])
{
	FILE *fpn;
	fpn= fopen("FILE20.dat","w");
  int i,j;
  int count1,count2,count4,count8,count16,count32,count64,count128,count256,count512 = 0;
  for(i=64-1;i>=0;i--){
    for(j=0;j<64;j++)
      fprintf(fpn," %d\n",i,lattice[i*64+j]+1);
  }

  fclose(fpn);
}
void crackspan(int t, int ncl, int NNcrack[][4])
{	
	FILE *fp2;
	fp2 = fopen("FILE2.dat","w");
	int i,j,ibond = 0;
	for(i=0;i<npart;i++)
    for(j=0;j<4;j++)if(NNcrack[i][j]>0)
		fprintf(fp2,"%d %d %d\n",i*4+j+1,NNcrack[i][j],i+1);
		ibond++;
	fclose(fp2);
}
double print_spanningcrack(int t,int ncl,int NNcrack[][4])
{
  FILE *fp1;
  char FILE1[50];
  int i,j,ibond=0;

  for(i=0;i<npart;i++)
    for(j=0;j<4;j++)if(NNcrack[i][j]>0)ibond++;

  
    
  fp1=fopen("FILE.dat","w");
  fprintf(fp1,"%s LAMMPS data file\n","#");
  fprintf(fp1,"%d atoms\n",npart);
  fprintf(fp1,"%d bonds\n",ibond);
  fprintf(fp1,"1 atom types\n");
  fprintf(fp1,"%d bond types\n",ncl+1);
  fprintf(fp1,"0.0 %e xlo xhi\n",L);
  fprintf(fp1,"0.0 %e ylo yhi\n",L);
  fprintf(fp1,"0.0 0.001 zlo zhi\n");

  fprintf(fp1,"\n");
  fprintf(fp1,"Atoms %s bond\n","#");
  fprintf(fp1,"\n");
    for(i=0;i<npart;i++)
    fprintf(fp1,"%d 1 1 %g %g\n",i+1,xi[i],yi[i]);
 
  
  fprintf(fp1,"\n");
  fprintf(fp1,"Bonds\n");
  fprintf(fp1,"\n");
  
 
  for(i=0;i<npart;i++)
    for(j=0;j<4;j++)
      if(NNcrack[i][j]>0)
	fprintf(fp1,"%d %d %d %d\n",i*4+j+1,NNcrack[i][j],i+1,neigh2[i][j]+1);
  fclose(fp1);
}