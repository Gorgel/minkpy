
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include "main.h"
#include "fourier.h"
#include "org.h"

#define pi 3.141592653589793238462643

void minkowski(void*, void*, int, float, float, void*);
void minkowski_calculation(float*,float*,parlist*,int*);

int binko(int,int);
float errorf(float);
float gaussf(float);
float gam(float);
float hermite(int,float);
float omega(int);
int round2(float);
int idxtwo(int,int,int,int);
int idx(int,int,int,int*);
int idx2(int,int,int,int*);

float invariant(float*,int,int,int,float*,float*,parlist*);
void anotherslice(int,float***,float***,float***,float*,parlist*);
void derive(float**,float**,float**,parlist*);
void kollekte(float,float**,float***,float*,float*,float*,float*,parlist*);

void minkowski(void *uv, void *usizev, int numbins, float low, float high, void *vsumv)
{
  parlist par;

  /* cast the void pointers to floats */
  float *u = (float*)uv; /* the data to be analyzed */
  long *uvsize = (long*)usizev; /* size of the 3D uv array */ 
  float *vsum = (float*)vsumv; /* the container of the results */

  int*vsumsize; /* will contain a 1D array with the sizes of the vsum array */
  int nvsum=0; /* will contain total size of the vsum array */

  int i,j,k;

  /* set default values in the parameter structure*/
  par.sigma=2; par.seed=200294;
  par.dim=intvector(1,3); par.dim[1]=0; par.dim[2]=0; par.dim[3]=0; 
  empty(&par.inname); empty(&par.outname);
  empty(&par.colname); empty(&par.grsname);
  par.length=1; par.nongauss=0; par.lo=0; par.hi=0; par.bins=0; 
  par.med=0; par.pixel=1; par.normal=1; par.time=0;
  par.cstyle  =0; par.floattyp=0; par.complement = 0;

  /* Assign array size (dimensions) */
  par.dim[1]=uvsize[0]; par.dim[2]=uvsize[1]; par.dim[3]=uvsize[2];
  par.a=1.0/(float)par.dim[1]; /* grid constant is used everywhere */

  /* Assign binning values to the parameter structure */
  par.bins=numbins;
  par.lo=low;
  par.hi=high;

  vsumsize=intvector(1,3); /* allocate size array */
  vsumsize[1]=2;vsumsize[2]=par.bins+1;vsumsize[3]=5; /* fill size array */
  nvsum=vsumsize[1]*vsumsize[2]*vsumsize[3]; /* total size of the vsum array */
  //vsum = (float*)calloc(nvsum,sizeof(float)); /* allocate vsum array */
  /* initialize result to zero */
  for(i=0;i<vsumsize[1];i++) for(j=0;j<vsumsize[2];j++) 
			       for(k=0;k<vsumsize[3];k++) {
				 vsum[idx2(i,j,k,vsumsize)]=0;
			       }
  /* call Minkowski routine */
  minkowski_calculation(u,vsum,&par,vsumsize);  
 
}

void minkowski_calculation(float*u,float*vsum,parlist*par,int*vsumsize)
{
  int index,i,j,k,x,y,z;
  int *d=par->dim,l=d[1]*d[2]*d[3];
  float a=par->a,b=par->bins/(par->hi-par->lo);

  static float**un,**up,**ud,*ucube[8];
  static float*thr,***v;
  static float mu,sig,tau;

  /* allocate storage space */
  un=matrix(0,d[2]*d[3],0,9);
  up=matrix(0,d[2]*d[3],0,9);
  ud=matrix(0,d[2]*d[3],0,9);
  
  /* define and fill threshold vector */
  thr=vector(0,par->bins);
  for(j=0;j<=par->bins;j++) thr[j]=par->lo+j/b;
  
  /* define data structure for results */
  v   =(float***)malloc(4*sizeof(float**));
  for(i=0;i<3;i++) {
    v[i]=matrix(0,par->bins,0,3);
    for(j=0;j<=par->bins;j++) for(k=0;k<4;k++) {
	v[i][j][k]=0;
      }
  }
  
  /* empty a few things */
  mu=sig=tau=0;

  /* note that we have to start with one, so the routine anotherslice
     can move from zero! */
  
  for(x=1;x<=d[1];x++) {

    /* get another slice from the random field, together with derivatives */
    anotherslice(x,&un,&up,&ud,u,par);
    
    for(y=0,i=0;y<d[2];y++) for(z=0;z<d[3];z++,i++) {
	
	/* corners of cube */
	/* 1 - x-axis, 2 - y-axis, 4 - z-axis, sums - diagonals */
	ucube[0]=un[i];
	ucube[1]=up[i];
	ucube[2]=un[idxtwo(y+1,z  ,d[2],d[3])];
	ucube[3]=up[idxtwo(y+1,z  ,d[2],d[3])];
	ucube[4]=un[idxtwo(y  ,z+1,d[2],d[3])];
	ucube[5]=up[idxtwo(y  ,z+1,d[2],d[3])];
	ucube[6]=un[idxtwo(y+1,z+1,d[2],d[3])];
	ucube[7]=up[idxtwo(y+1,z+1,d[2],d[3])];
	
	/* collect statistics */
	kollekte(1.,ucube,v,&mu,&sig,&tau,thr,par);
      }
  }
  
  /* -------------------------------------------- */
  /* finish calculations of Minkowski functionals */
  /* -------------------------------------------- */

  for(j=0;j<=par->bins;j++) {
    for(k=0;k<4;k++) {
    /* Crofton's formula over lattice planes -- Minkowski functionals */
      v[0][j][k]/=l*pow(a,k);  
    }
    /* Koenderink invariant average -- normalization */
    for(k=0;k<4;k++) v[1][j][k]/=l;
    for(k=1;k<4;k++) v[1][j][k]*=b/(3*omega(k));
  }

  /* collect results */
  for(i=0;i<2;i++) for(j=0;j<=par->bins;j++) for(k=0;k<4;k++) {
	vsum[idx2(i,j,k,vsumsize)]=v[i][j][k];
      }
  /* put thresholds into result array */
  for(i=0;i<2;i++) for(j=0;j<=par->bins;j++) {
      vsum[idx2(i,j,4,vsumsize)]=thr[j];
    }
}

int binko(int n,int k)
{
  if(k==0||k==n||n==0||n==1) return 1;
  if(k==1||k==n-1) return n;
  return binko(n-1,k-1)+binko(n-1,k);
}
float errorf(float x)
{
  if(x> 10) return  1; if(x<-10) return -1;
  if(fabs(x)<1e-6) return x*2/sqrt(pi);
  return qromb(gaussf,0,x);
}
float gaussf(float x)
{
  return 2/sqrt(pi)*exp(-x*x);
}
float gam(float x)
{
  if(x<=0)  return -666;
  if(x==.5) return sqrt(pi);
  if(x==1)  return 1;
  return (x-1)*gam(x-1);
}
float hermite(int n,float x)
{
  if(n<0)  return 0;
  if(n==0) return exp(-.5*x*x)/sqrt(2*pi);
  return x*hermite(n-1,x)-(n-1)*hermite(n-2,x);
}
float omega(int k)
{
  switch(k) {
  case 0 : return 1;
  case 1 : return 2;
  case 2 : return pi;
  case 3 : return pi/.75;
  default : return pow(pi,k/2.)/gam(1+k/2.);
  }
}
int round2(float x)
{
  if(x>=0) return (int)(x+.5); else return (int)(x-.5);
}
int idxtwo(int x,int y,int dx, int dy)
{
  while(x>=dx) x-=dx; while(x<0) x+=dx;
  while(y>=dy) y-=dy; while(y<0) y+=dy;
  return x*dy+y;
}

int idx(int i, int j, int k, int*dim)
{
  while(i>=dim[1]) i-=dim[1]; while(i<0) i+=dim[1];
  while(j>=dim[2]) j-=dim[2]; while(j<0) j+=dim[2];
  while(k>=dim[3]) k-=dim[3]; while(k<0) k+=dim[3];
  return ((i*dim[2]+j)*dim[3]+k);
}

int idx2(int i, int j, int k, int*dim)
{
  return ((i*dim[2]+j)*dim[3]+k);
}

/* calculate the values of a random field uu in a slice um at position
   x also extract the two following slices un and um, and calculate
   first and second derivatives in um and un */

void anotherslice(int x,float***um,float***un,float***up,float*u,parlist*par)
{
  int i,xx,y,z,*d=par->dim,l=d[1]*d[2]*d[3];
  float**h;

  static int num=0,xnul,xone,xtwo;
  static float**unul,**uone,**utwo;

  /* if we are already moving through the field, (*um) contains the
     slice at x-1, (*un) the slice at x, and (*up) the slice at
     x+1. We already have derivatives at x.  If we are only starting,
     the slices are empty and we have to extract the values and
     calculate the derivatives before we can set off. */

  if(x==1) {
    /* save the first batch of slices for later */
    xnul=x-1; while(xnul>=d[1]) xnul-=d[1]; while(xnul<0) xnul+=d[1];
    xone=x;   while(xone>=d[1]) xone-=d[1]; while(xone<0) xone+=d[1];
    xtwo=x+1; while(xtwo>=d[1]) xtwo-=d[1]; while(xtwo<0) xtwo+=d[1];

    if(num==0) {
      unul=matrix(0,d[2]*d[3],0,9);
      uone=matrix(0,d[2]*d[3],0,9);
      utwo=matrix(0,d[2]*d[3],0,9);
    }

    for(y=0,i=0;y<d[2];y++) for(z=0;z<d[3];z++,i++) {
      num+=3;
      (*um)[i][0]=unul[i][0]=u[idx(x-1,y,z,d)];
      (*un)[i][0]=uone[i][0]=u[idx(x  ,y,z,d)];
      (*up)[i][0]=utwo[i][0]=u[idx(x+1,y,z,d)];
    }

    derive(*un,*up,*um,par);

  }

  /* Now we can extract the next slice from the field, move the
     other two downwards, and calculate the missing derivatives */

  h=(*um); (*um)=(*un); (*un)=(*up); (*up)=h;
  xx=x+2; while(xx>=d[1]) xx-=d[1]; while(xx<0) xx+=d[1];
  if(xx==xnul) for(i=0;i<d[2]*d[3];i++) (*up)[i][0]=unul[i][0];
  else if(xx==xone) for(i=0;i<d[2]*d[3];i++) (*up)[i][0]=uone[i][0];
  else if(xx==xtwo) for(i=0;i<d[2]*d[3];i++) (*up)[i][0]=utwo[i][0];
  else 
    for(y=0,i=0;y<d[2];y++) for(z=0;z<d[3];z++,i++) 
      num++, (*up)[i][0]=u[idx(x+2,y,z,d)];

  derive(*un,*up,*um,par);
}

void derive(float**un,float**up,float**um,parlist*par)
{
  int i,y,z,*d=par->dim;
  float a=par->a;
  int inn,inp,inm,ipn,ipp,ipm,imn,imp,imm;
  float unnn,upnn,unpn,unnp,umnn,unmn,unnm;
  float uppn,unpp,upnp,ummn,unmm,umnm,upmn,umpn,unpm,unmp,umnp,upnm, uppp;

  for(y=0,i=0;y<d[2];y++) for(z=0;z<d[3];z++,i++) {
    if(z) {
      inp=idxtwo(y  ,z+1,d[2],d[3]);
      ipn=idxtwo(y+1,z  ,d[2],d[3]);
      ipp=idxtwo(y+1,z+1,d[2],d[3]);
      imn=idxtwo(y-1,z  ,d[2],d[3]);
      imp=idxtwo(y-1,z+1,d[2],d[3]);

      ummn=um[imn][0];
      umnm=umnn; umnn=umnp; umnp=um[inp][0];
      umpn=um[ipn][0];
      unmm=unmn; unmn=unmp; unmp=un[imp][0];
      unnm=unnn; unnn=unnp; unnp=un[inp][0];
      unpm=unpn; unpn=unpp; unpp=un[ipp][0];
      upmn=up[imn][0];
      upnm=upnn; upnn=upnp; upnp=up[inp][0];
      uppn=uppp; uppp=up[ipp][0];
    } else {
      inn=idxtwo(y  ,z  ,d[2],d[3]);
      inp=idxtwo(y  ,z+1,d[2],d[3]);
      inm=idxtwo(y  ,z-1,d[2],d[3]);
      ipn=idxtwo(y+1,z  ,d[2],d[3]);
      ipp=idxtwo(y+1,z+1,d[2],d[3]);
      ipm=idxtwo(y+1,z-1,d[2],d[3]);
      imn=idxtwo(y-1,z  ,d[2],d[3]);
      imp=idxtwo(y-1,z+1,d[2],d[3]);
      imm=idxtwo(y-1,z-1,d[2],d[3]);

      unnn=un[inn][0]; 
      upnn=up[inn][0]; unpn=un[ipn][0]; unnp=un[inp][0];
      umnn=um[inn][0]; unmn=un[imn][0]; unnm=un[inm][0];
      uppn=up[ipn][0]; unpp=un[ipp][0]; upnp=up[inp][0];
      ummn=um[imn][0]; unmm=un[imm][0]; umnm=um[inm][0];
      upmn=up[imn][0]; umpn=um[ipn][0]; unpm=un[ipm][0];
      unmp=un[imp][0]; umnp=um[inp][0]; upnm=up[inm][0];
      uppp=up[ipp][0];
    }
    
    un[i][1]=(upnn-umnn)/(2*a);
    un[i][2]=(unpn-unmn)/(2*a);
    un[i][3]=(unnp-unnm)/(2*a);
    un[i][4]=(upnn-2*unnn+umnn)/(a*a);
    un[i][5]=(unpn-2*unnn+unmn)/(a*a);
    un[i][6]=(unnp-2*unnn+unnm)/(a*a);
    un[i][7]=(uppn-upmn-umpn+ummn)/(4*a*a);
    un[i][8]=(unpp-unpm-unmp+unmm)/(4*a*a);
    un[i][9]=(upnp-upnm-umnp+umnm)/(4*a*a);
  }
}


/* given the values of a random field and its derivatives at the
   corners of a cube (1 - x-axis, 2 - y-axis, 4 - z-axis, sums -
   diagonals), collect the statistics for Minkowski functional
   calculation. if necessary, interpolate new values and go to finer
   cubes */

void kollekte(float weight,float**u,float***v,
	      float*mu,float*sig,float*tau,float*thr,parlist*par)
{
  int ii,i,j,k,l,m,n[4];
  float w[4],x,y,z,xx,yy,zz,p,q,r,**uu,**c,*o,inv[4],t,min,max;
  float a=par->a*pow(weight,1./3.),b=par->bins/(par->hi-par->lo);

  /* Koenderink invariants at origin */
  /* derivatives: 0=u 1=ux 2=uy 3=uz 4=up 5=uq 6=ur 7=uxy 8=uyz 9=uzx */
  o=u[0];
  inv[0]=o[0]; 
  inv[1]=sqrt(o[1]*o[1]+o[2]*o[2]+o[3]*o[3]);
  inv[2]=
    (2*o[1]*o[7]*o[2]-o[4]*(o[2]*o[2]+o[3]*o[3])+
     2*o[2]*o[8]*o[3]-o[5]*(o[3]*o[3]+o[1]*o[1])+
     2*o[3]*o[9]*o[1]-o[6]*(o[1]*o[1]+o[2]*o[2])) 
      / (2*(o[1]*o[1]+o[2]*o[2]+o[3]*o[3]));
  inv[3]=
    (o[1]*o[1]*(o[5]*o[6]-o[8]*o[8])+2*o[1]*o[2]*(o[8]*o[9]-o[7]*o[6])+
     o[2]*o[2]*(o[6]*o[4]-o[9]*o[9])+2*o[2]*o[3]*(o[9]*o[7]-o[8]*o[4])+
     o[3]*o[3]*(o[4]*o[5]-o[7]*o[7])+2*o[3]*o[1]*(o[7]*o[8]-o[9]*o[5]))
      / pow(o[1]*o[1]+o[2]*o[2]+o[3]*o[3],1.5);

  /* minimum and maximum of field at corners */
  for(j=0,min=66e6,max=-66e6;j<8;j++) { 
    if(u[j][0]<min) min=u[j][0]; if(u[j][0]>max) max=u[j][0];
  }
  
  if(weight>pow(.5,3.*par->med)) {

    /* allocate memory */
    uu=matrix(0,7,0,9); c=matrix(0,9,0,7);

    /* calculate coefficients for interpolation */
    for(i=0;i<10;i++) {
      c[i][0]=u[0][i];
      c[i][1]=u[1][i]-u[0][i];
      c[i][2]=u[2][i]-u[0][i];
      c[i][4]=u[4][i]-u[0][i];
      c[i][3]=u[3][i]-u[1][i]-u[2][i]+u[0][i];
      c[i][6]=u[6][i]-u[2][i]-u[4][i]+u[0][i];
      c[i][5]=u[5][i]-u[4][i]-u[1][i]+u[0][i];
      c[i][7]=u[7][i]-u[3][i]-u[6][i]-u[5][i]+u[1][i ]+u[2][i]+u[4][i]-u[0][i];
    }




    for(p=0;p<1;p+=.5) for(q=0;q<1;q+=.5) for(r=0;r<1;r+=.5) {

      /* interpolate values */
      for(x=p,l=0;x<p+1;x+=.5) for(y=q;y<q+1;y+=.5) for(z=r;z<r+1;z+=.5,l++) {
	for(ii=0;ii<10;ii++) {
	  uu[l][ii]=0, m=0;
	  for(i=0,xx=1;i<=1;i++,xx*=x) 
	    for(j=0,yy=1;j<=1;j++,yy*=y) 
	      for(k=0,zz=1;k<=1;k++,zz*=z) 
		uu[l][ii]+=c[ii][m++]*xx*yy*zz;
	}
      }


   

      /* descend to finer cube */
      kollekte(weight/8,uu,v,mu,sig,tau,thr,par);
    }

    /* free memory for cube */
    free_matrix(uu,0,7,0,9); free_matrix(c,0,9,0,7);

  } else {

    /* ------------------------------------ */
    /* collect statistics over random field */
    /* ------------------------------------ */
  
    /* Crofton's formula over lattice planes -- plaquette counts */
    for(i=0;i<4;i++) w[i]=pow(weight,1.-i/3.);
    for(j=0;j<=par->bins;j++) {
      if(thr[j]<min)
	v[0][j][0]+=w[0];
      else if(thr[j]<max) {
	t=thr[j];
	n[3]=u[0][0]>t;
	n[2]=((u[0][0]>t||u[1][0]>t)+
	      (u[0][0]>t||u[2][0]>t)+
	      (u[0][0]>t||u[4][0]>t));
	n[1]=((u[0][0]>t||u[1][0]>t||u[2][0]>t||u[3][0]>t)+
	      (u[0][0]>t||u[2][0]>t||u[4][0]>t||u[6][0]>t)+
	      (u[0][0]>t||u[4][0]>t||u[1][0]>t||u[5][0]>t));
	n[0]=u[0][0]>t||u[1][0]>t||u[2][0]>t||u[3][0]>t||u[4][0]>t||u[5][0]>t||u[6][0]>t||u[7][0]>t;
	v[0][j][0] += w[0]*    n[3];
	v[0][j][1] += w[1]*(-3*n[3]+  n[2])          *2./9.;
	v[0][j][2] += w[2]*( 3*n[3]-2*n[2]+n[1])     *2./9.;
	v[0][j][3] += w[3]*(  -n[3]+  n[2]-n[1]+n[0]);

 
      } 
    }

    /* Koenderink invariant average -- collect invariants in bins */
    for(j=0;inv[0]>=thr[j]&&j<=par->bins;j++) v[1][j][0]+=weight; 
    j=round2(b*(inv[0]-par->lo));
    if(j>=0&&j<=par->bins) for(k=1;k<4;k++) v[1][j][k]+=weight*inv[k];
    
    /* exact Crofton's formula -- forget it, grid is excellent*/
    ;
    
    /* Tomita's formulae -- estimate parameters from field */
    *mu+=weight*inv[0],*sig+=weight*inv[0]*inv[0],*tau+=weight*inv[1]*inv[1];

  }

  return;
}
