# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <inttypes.h>
# include "defs.h"
# include "sys.h"
# include "global.h"
/****************************************************************************/
static double xt[NBA],yt[NBA],zt[NBA];
static double bx[NBA+1],by[NBA+1],bz[NBA+1];
static double ux[NBA],uy[NBA],uz[NBA];
static double th[3];
static double aO,gO,dO;
static double thH,thO;
static double aH,gH,dH;
static double phCb,thCb;           
static double phHa,thHa,phHa2;
static double aHa,gHa,dHa;
static double aCb,gCb,dCb;
/****************************************************************************/
void agd(double th1, double th2,double b,double *a,double *g,double *d);
void up(int ia,int ic);
void dw(int ia,int ic);
/****************************************************************************/
/***** CHAIN CONSTRUCTION ***************************************************/
/****************************************************************************/
void int2cart_all(int dir) 
{
  int i;

  if (dir<0) 
    for (i=0;i<NCH;i++) int2cart(Ne[i],-1);
  else 
    for (i=0;i<NCH;i++) int2cart(Nb[i],1);
}
/****************************************************************************/
void agd(double th1, double th2,double b,double *a,double *g,double *d) 
{
  (*a)=-b*sin(th2)/sin(th1);
  (*g)=-b*cos(th2);
  (*d)=-b*sin(th2)/tan(th1);
}
/****************************************************************************/
void up(int ia,int ic) 
/* construct backbone N -> C */
{
  int i,k;
  double sph,cph,renorm;
  static double ab[3],gb[3],db[3];
  
  if (ia < 0) {
    for (i=0;i<3;i++) agd(th[(i+1)%3],th[(i+2)%3],1,ab+i,gb+i,db+i);
    return ;
  }
  
  k=3*ia;
  xt[k]=x[iN[ia]]; 
  yt[k]=y[iN[ia]]; 
  zt[k]=z[iN[ia]];
  xt[k+1]=x[iCa[ia]]; 
  yt[k+1]=y[iCa[ia]]; 
  zt[k+1]=z[iCa[ia]];    
  bx[k+1]=xt[k+1]-xt[k];
  by[k+1]=yt[k+1]-yt[k];
  bz[k+1]=zt[k+1]-zt[k];
  renorm=sqrt(bx[k+1]*bx[k+1]+by[k+1]*by[k+1]+bz[k+1]*bz[k+1]);
  bx[k+1]/=renorm; by[k+1]/=renorm; bz[k+1]/=renorm;
  if (ia!=Nb[ic]) {
    bx[k]=x[iN[ia]]-x[iC[ia-1]];
    by[k]=y[iN[ia]]-y[iC[ia-1]];
    bz[k]=z[iN[ia]]-z[iC[ia-1]];
  } else {
    bx[k]=(x[iN[ia]]-x[iN[ia]+1]+(gH-dH)*bx[k+1])/aH;
    by[k]=(y[iN[ia]]-y[iN[ia]+1]+(gH-dH)*by[k+1])/aH;
    bz[k]=(z[iN[ia]]-z[iN[ia]+1]+(gH-dH)*bz[k+1])/aH;
  }
  renorm=sqrt(bx[k]*bx[k]+by[k]*by[k]+bz[k]*bz[k]);
  bx[k]/=renorm; by[k]/=renorm; bz[k]/=renorm;
  
  for (i=3*ia+2;i<=3*Ne[ic]+3;i++) { /* backbone */
    k=i%3;
    ux[i-2]=by[i-2]*bz[i-1]-by[i-1]*bz[i-2]; 
    uy[i-2]=bz[i-2]*bx[i-1]-bz[i-1]*bx[i-2]; 
    uz[i-2]=bx[i-2]*by[i-1]-bx[i-1]*by[i-2];
    cph=cos(ph[i-2]); sph=sin(ph[i-2]);
    bx[i]=ab[k]*cph*bx[i-2]+(gb[k]+db[k]*cph)*bx[i-1]-ab[k]*sph*ux[i-2];
    by[i]=ab[k]*cph*by[i-2]+(gb[k]+db[k]*cph)*by[i-1]-ab[k]*sph*uy[i-2];
    bz[i]=ab[k]*cph*bz[i-2]+(gb[k]+db[k]*cph)*bz[i-1]-ab[k]*sph*uz[i-2];
    renorm=sqrt(bx[i]*bx[i]+by[i]*by[i]+bz[i]*bz[i]);
    bx[i]/=renorm; by[i]/=renorm; bz[i]/=renorm;
    if (i<=3*Ne[ic]+2) {  
      xt[i]=xt[i-1]+b[k]*bx[i]; 
      yt[i]=yt[i-1]+b[k]*by[i]; 
      zt[i]=zt[i-1]+b[k]*bz[i];
    }
  }
  
  return ;
}
/****************************************************************************/
void dw(int ia,int ic) 
/* construct backbone C -> N */
{
  int i,k;
  double sph,cph,renorm; 
  static double ab[3],gb[3],db[3];
  
  if (ia < 0) {
    for (i=0;i<3;i++) agd(th[(i+1)%3],th[i%3],1,ab+i,gb+i,db+i);
    return ;
  }
  
  k=3*ia; 
  xt[k+2]=x[iC[ia]]; 
  yt[k+2]=y[iC[ia]]; 
  zt[k+2]=z[iC[ia]];
  xt[k+1]=x[iCa[ia]]; 
  yt[k+1]=y[iCa[ia]]; 
  zt[k+1]=z[iCa[ia]];    
  bx[k+2]=xt[k+2]-xt[k+1];
  by[k+2]=yt[k+2]-yt[k+1];
  bz[k+2]=zt[k+2]-zt[k+1];
  renorm=sqrt(bx[k+2]*bx[k+2]+by[k+2]*by[k+2]+bz[k+2]*bz[k+2]);
  bx[k+2]/=renorm; by[k+2]/=renorm; bz[k+2]/=renorm;
  if (ia!=Ne[ic]) {
    bx[k+3]=x[iN[ia+1]]-x[iC[ia]];
    by[k+3]=y[iN[ia+1]]-y[iC[ia]];
    bz[k+3]=z[iN[ia+1]]-z[iC[ia]];
  } else {
    bx[k+3]=-(x[iC[ia]]-x[iC[ia]+1]-aO*bx[k+2])/(gO-dO);
    by[k+3]=-(y[iC[ia]]-y[iC[ia]+1]-aO*by[k+2])/(gO-dO);
    bz[k+3]=-(z[iC[ia]]-z[iC[ia]+1]-aO*bz[k+2])/(gO-dO);
  }
  renorm=sqrt(bx[k+3]*bx[k+3]+by[k+3]*by[k+3]+bz[k+3]*bz[k+3]);
  bx[k+3]/=renorm; by[k+3]/=renorm; bz[k+3]/=renorm;
  
  for (i=3*ia+1;i>=3*Nb[ic]-1;i--) { /* backbone */
    ux[i+1]=by[i+1]*bz[i+2]-by[i+2]*bz[i+1]; 
    uy[i+1]=bz[i+1]*bx[i+2]-bz[i+2]*bx[i+1]; 
    uz[i+1]=bx[i+1]*by[i+2]-bx[i+2]*by[i+1];
    k=i%3; /* N 1 ; C 0; Ca 2 */
    if (i>=3*Nb[ic]) {
      cph=cos(ph[i]); sph=sin(ph[i]);
      bx[i]=ab[k]*cph*bx[i+2]+(gb[k]+db[k]*cph)*bx[i+1]-ab[k]*sph*ux[i+1];
      by[i]=ab[k]*cph*by[i+2]+(gb[k]+db[k]*cph)*by[i+1]-ab[k]*sph*uy[i+1];
      bz[i]=ab[k]*cph*bz[i+2]+(gb[k]+db[k]*cph)*bz[i+1]-ab[k]*sph*uz[i+1];
      renorm=sqrt(bx[i]*bx[i]+by[i]*by[i]+bz[i]*bz[i]);
      bx[i]/=renorm; by[i]/=renorm; bz[i]/=renorm;
    }
    if (i>=3*Nb[ic]+1) {
      xt[i-1]=xt[i]-b[k]*bx[i]; 
      yt[i-1]=yt[i]-b[k]*by[i]; 
      zt[i-1]=zt[i]-b[k]*bz[i];
    }
  }
  return ;
}
/****************************************************************************/
void int2cart(int ia,int dir)
{
  /* dir = 0, init */
  /* dir > 0, reconstruct chain from ia to C-term */
  /* dir < 0, reconstruct chain from ia to N-term */
  int i,j,k,ic,ibeg,iend;
  double sph,cph;

  if (NCH == 0) return ;
  
  if (dir==0) {                           
    th[0]=121.7*pi/180; //N
    th[1]=111.0*pi/180; //Ca
    th[2]=116.6*pi/180; //C
    thCb=110.0*pi/180;
    phCb=-120.9*pi/180;
    thHa=109.0*pi/180;
    phHa=118.7*pi/180;
    phHa2=-118.7*pi/180;
    thH=th[0]/2;
    thO=th[2]/2;
    agd(th[0],thHa,bCaHa,&aHa,&gHa,&dHa);
    agd(th[0],thCb,bCaCb,&aCb,&gCb,&dCb);
    agd(th[0],thH,bNH,&aH,&gH,&dH);
    agd(th[2],thO,bCO,&aO,&gO,&dO);
    for (j=0;j<NCH;j++) {
      i=Nb[j];
      x[iN[i]]=y[iN[i]]=z[iN[i]]=0;
      x[iCa[i]]=b[1]; y[iCa[i]]=z[iCa[i]]=0;       
      x[iN[i]+1]=bNH*cos(pi-th[0]/2); 
      y[iN[i]+1]=bNH*sin(pi-th[0]/2); 
      z[iN[i]+1]=0;
    }
    up(-1,0);
    dw(-1,0);
    return;
  }

  for (i=0;i<NBA;i++) xt[i]=yt[i]=zt[i]=0;
  for (i=0;i<NBA+1;i++) bx[i]=by[i]=bz[i]=0;
  for (i=0;i<NBA;i++) ux[i]=uy[i]=uz[i]=0;
  
  ic=a2c[ia];  /* chain number */

  if (dir>0) { 
    up(ia,ic); 
    ibeg=ia; iend=Ne[ic];
  } else { 
    dw(ia,ic); 
    ibeg=Nb[ic]; iend=ia;
  }

  for (i=ibeg;i<=iend;i++) { /* backbone */
    k=3*i; 
    x[iN[i]]=xt[k]; x[iCa[i]]=xt[k+1]; x[iC[i]]=xt[k+2];    
    y[iN[i]]=yt[k]; y[iCa[i]]=yt[k+1]; y[iC[i]]=yt[k+2];     
    z[iN[i]]=zt[k]; z[iCa[i]]=zt[k+1]; z[iC[i]]=zt[k+2];
  }
  
  for (i=ibeg;i<=iend;i++) { 
    j=iCa[i]; k=3*i; /* Ha */
    cph=cos(ph[k]+phHa); sph=sin(ph[k]+phHa);
    x[j+1]=x[j]+aHa*cph*bx[k]+(gHa+dHa*cph)*bx[k+1]-aHa*sph*ux[k];
    y[j+1]=y[j]+aHa*cph*by[k]+(gHa+dHa*cph)*by[k+1]-aHa*sph*uy[k];
    z[j+1]=z[j]+aHa*cph*bz[k]+(gHa+dHa*cph)*bz[k+1]-aHa*sph*uz[k];
    if (seq[i]!='G') { /* Cb */
      cph=cos(ph[k]+phCb); sph=sin(ph[k]+phCb);
      x[j+2]=x[j]+aCb*cph*bx[k]+(gCb+dCb*cph)*bx[k+1]-aCb*sph*ux[k];
      y[j+2]=y[j]+aCb*cph*by[k]+(gCb+dCb*cph)*by[k+1]-aCb*sph*uy[k];
      z[j+2]=z[j]+aCb*cph*bz[k]+(gCb+dCb*cph)*bz[k+1]-aCb*sph*uz[k];
      }
    else { /* Ha2 */
      cph=cos(ph[k]+phHa2); sph=sin(ph[k]+phHa2);
      x[j+2]=x[j]+aHa*cph*bx[k]+(gHa+dHa*cph)*bx[k+1]-aHa*sph*ux[k];
      y[j+2]=y[j]+aHa*cph*by[k]+(gHa+dHa*cph)*by[k+1]-aHa*sph*uy[k];
      z[j+2]=z[j]+aHa*cph*bz[k]+(gHa+dHa*cph)*bz[k+1]-aHa*sph*uz[k];
    }
    j=iN[i]; /* H */
    if (seq[i]=='P' || (i==ibeg && dir>0)) {}
    else {
      x[j+1]=x[j]-aH*bx[k]+(gH-dH)*bx[k+1];
      y[j+1]=y[j]-aH*by[k]+(gH-dH)*by[k+1];
      z[j+1]=z[j]-aH*bz[k]+(gH-dH)*bz[k+1];
    }
    j=iC[i]; /* O */
    if (i==iend && dir<0) {}
    else {
      x[j+1]=x[j]-aO*bx[k+2]+(gO-dO)*bx[k+3];
      y[j+1]=y[j]-aO*by[k+2]+(gO-dO)*by[k+3];
      z[j+1]=z[j]-aO*bz[k+2]+(gO-dO)*bz[k+3];
    }
  }
}
/****************************************************************************/
/*********************Functions to radius of gyration*********************************/
double rgyr(int ic) {
  double xcm,ycm,zcm;
  return cm(&xcm,&ycm,&zcm,ic,0);
}
/****************************************************************************/
double cm(double *xcm,double *ycm,double *zcm,int ic,int offset) 
{
  int i,j,n=0;
  double rg2=0;

  if (ic > NCH-1 || NCH == 0) return 0;

  (*xcm)=(*ycm)=(*zcm)=0;

  for (i=Nb[ic]+offset;i<=Ne[ic]-offset;i++) {
    j=iCa[i];
    (*xcm)+=x[j];
    (*ycm)+=y[j];
    (*zcm)+=z[j];
    n++;
  }

  (*xcm)/=n; (*ycm)/=n; (*zcm)/=n;

  for (i=Nb[ic]+offset;i<=Ne[ic]-offset;i++) {
    j=iCa[i];
    rg2+=(((*xcm)-x[j])*((*xcm)-x[j])+
	  ((*ycm)-y[j])*((*ycm)-y[j])+
	  ((*zcm)-z[j])*((*zcm)-z[j]));
  }

  return sqrt(rg2/n);
}
/*********************************************************************/
/******** Translate coordinates into original periodic box ***********/
/*********************************************************************/
void in2box(double *xb,double *yb,double *zb,double *x,double *y,double *z,
	    int n1,int n2) {
  int i;

  for (i = n1; i <= n2; i++) {
      xb[i]=x[i]; yb[i]=y[i]; zb[i]=z[i];
      while (xb[i]>=BOX) xb[i]-=BOX; while (xb[i]<0) xb[i]+=BOX;
      while (yb[i]>=BOX) yb[i]-=BOX; while (yb[i]<0) yb[i]+=BOX;
      while (zb[i]>=BOX) zb[i]-=BOX; while (zb[i]<0) zb[i]+=BOX;
  }
}
/****************************************************************************/
/*** Apply periodic boundary condition for atoms and crowder agents *********/
/****************************************************************************/
void ch2box()
{
  int i,j,k,nx,ny,nz;

  for (j=0;j<NCH;j++) {
    nx=ny=nz=0;
    k=iBeg[j]+(iEnd[j]-iBeg[j])/2;
    while (x[k]+nx*BOX>BOX) nx--; while (x[k]+nx*BOX<0) nx++;
    while (y[k]+ny*BOX>BOX) ny--; while (y[k]+ny*BOX<0) ny++;
    while (z[k]+nz*BOX>BOX) nz--; while (z[k]+nz*BOX<0) nz++;
    for (i=iBeg[j];i<=iEnd[j];i++) {
      x[i]+=BOX*nx;
      y[i]+=BOX*ny;
      z[i]+=BOX*nz;
    }
  }

  for (j=0;j<NCR;j++) {
    nx=ny=nz=0;
    while (xcr[j] + nx * BOX > BOX) nx--; while (xcr[j] + nx*BOX < 0) nx++;
    while (ycr[j] + ny * BOX > BOX) ny--; while (ycr[j] + ny*BOX < 0) ny++;
    while (zcr[j] + nz * BOX > BOX) nz--; while (zcr[j] + nz*BOX < 0) nz++;
    xcr[j] += BOX*nx;
    ycr[j] += BOX*ny;
    zcr[j] += BOX*nz;
  }
  
  for (i=0;i<NTO;i++) {
    xo[i]=x[i]; 
    yo[i]=y[i]; 
    zo[i]=z[i];
  }
  
  for (i=0;i<NCR;i++) {
    xcro[i]=xcr[i]; 
    ycro[i]=ycr[i]; 
    zcro[i]=zcr[i];
  }
  
}
/****************************************************************************/
void cr2box()
{
  int i,j;

  for (j=0;j<NCR;j++) {
    while (xcr[j] >= BOX) xcr[j] -= BOX;  while (xcr[j] < 0) xcr[j] += BOX;
    while (ycr[j] >= BOX) ycr[j] -= BOX;  while (ycr[j] < 0) ycr[j] += BOX;
    while (zcr[j] >= BOX) zcr[j] -= BOX;  while (zcr[j] < 0) zcr[j] += BOX;
  }

  for (i=0;i<NCR;i++) {
    xcro[i]=xcr[i]; 
    ycro[i]=ycr[i]; 
    zcro[i]=zcr[i];
  }
  
}
/***********************************************************************/
/****Function to calculate the distance between crowders*******/
/***********************************************************************/
double dist(int i1,int i2) {
  double dx,dy,dz;
  
  dx = x[i1] - x[i2]; bc(&dx);
  dy = y[i1] - y[i2]; bc(&dy);
  dz = z[i1] - z[i2]; bc(&dz);

  return sqrt(dx*dx + dy*dy + dz*dz);
}
/****************************************************************************/
double dist_atom_crowd(int i1,int i2) {
  double dx,dy,dz;
  
  dx = x[i1] - xcr[i2]; bc(&dx);
  dy = y[i1] - ycr[i2]; bc(&dy);
  dz = z[i1] - zcr[i2]; bc(&dz);

  return sqrt(dx*dx + dy*dy + dz*dz);
}
/****************************************************************************/
double dist_cr(int ic1,int ic2) {
  double dx,dy,dz;
  
  dx = xcr[ic1] - xcr[ic2]; bc(&dx);
  dy = ycr[ic1] - ycr[ic2]; bc(&dy);
  dz = zcr[ic1] - zcr[ic2]; bc(&dz);

  return sqrt(dx*dx + dy*dy + dz*dz);
}
/****************************************************************************/
/****************Function for boundary condition**************/
void bc(double *x) 
{
  while ((*x)>boxhf) (*x)-=BOX;
  while ((*x)<-boxhf) (*x)+=BOX;
}
/****************************************************************************/
void backup_coord(){
  int i;
    
  for (i=0;i<NTO;i++) {
    xo[i] = x[i]; 
    yo[i] = y[i]; 
    zo[i] = z[i];
  }
  
  for (i=0;i<NCR;i++) {
    xcro[i] = xcr[i];
    ycro[i] = ycr[i];
    zcro[i] = zcr[i];
  }

  return ;
}
/****************************************************************************/
