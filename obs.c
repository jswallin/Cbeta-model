# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include "defs.h"
# include "sys.h"
# include "global.h"
/****************************************************************************/
/***** MEASUREMENTS *********************************************************/
/****************************************************************************/
double mindist_cr_at(){
  int i,j;
  double dx,dy,dz;
  double r2,r2min=1e4;

  for (i=0;i<NCR;i++) {
    for (j=0;j<NTO;j++) {
      dx = xcr[i] - x[j]; bc(&dx);
      dy = ycr[i] - y[j]; bc(&dy);
      dz = zcr[i] - z[j]; bc(&dz);
      r2 = dx *dx + dy *dy + dz *dz ;
      if (r2 < r2min) r2min = r2;
    }
  }

  return sqrt(r2min);
}
/****************************************************************************/
double mindist_cr_cr(){
  int i,j;
  double dx,dy,dz;
  double r2,r2min=1e4;

  for (i=0;i<NCR;i++) {
    for (j=i+1;j<NCR;j++) {
      dx = xcr[i] - xcr[j]; bc(&dx);
      dy = ycr[i] - ycr[j]; bc(&dy);
      dz = zcr[i] - zcr[j]; bc(&dz);
      r2 = dx *dx + dy *dy + dz *dz ;
      if (r2 < r2min) r2min = r2;
    }
  }

  return sqrt(r2min);
}
/****************************************************************************/
void crowd_stat(int iflag) {
  int i;
  static double xav[NCR],yav[NCR],zav[NCR];
  static double x2av[NCR],y2av[NCR],z2av[NCR];
  static double n;
  FILE *fp;
  
  if (iflag < 0) {
    n = 0;
    for (i = 0; i < NCR; i++) {
      xav[i] = yav[i] = zav[i] = 0;
      x2av[i] = y2av[i] = z2av[i] = 0;
    }
    return;
  }


  if (iflag == 0) {
    in2box(xcrb,ycrb,zcrb,xcr,ycr,zcr,0,NCR-1);
    for (i = 0; i < NCR; i++) {
      xav[i] += xcrb[i];
      yav[i] += ycrb[i];
      zav[i] += zcrb[i];
      x2av[i] += xcrb[i] * xcrb[i]; 
      y2av[i] += ycrb[i] * ycrb[i]; 
      z2av[i] += zcrb[i] * zcrb[i]; 
    }
    n++;
    return ;
  }


  if (iflag > 0) {
    fp = fopen("crowdstat","w");

    for (i = 0; i < NCR ; i++) {
      fprintf(fp,"%2i %lf %lf %lf %lf %lf %lf \n",i,
	      xav[i]/n,sqrt(x2av[i]/n - xav[i]*xav[i]/n/n),
	      yav[i]/n,sqrt(y2av[i]/n - yav[i]*yav[i]/n/n),
	      zav[i]/n,sqrt(z2av[i]/n - zav[i]*zav[i]/n/n) );
    }

    fclose(fp);
    return;
  }

  return;
}
/****************************************************************************/
double ree(int ich) {
  /* returns distance between first and last Ca atoms in chain ich */

  if (NCH == 0) return 0;

  int i1 = iCa[Nb[ich]];
  int i2 = iCa[Ne[ich]];
  
  return sqrt( (x[i2] - x[i1]) * (x[i2] - x[i1]) +
	       (y[i2] - y[i1]) * (y[i2] - y[i1]) +
	       (z[i2] - z[i1]) * (z[i2] - z[i1]) );
}
/****************************************************************************/
/**** Functions to calculate maximum size of clustering *********************/
/****************************************************************************/
int depth(int clr_no,int icur,int *alist) 
{
  int i,size=0;  

  alist[icur]=clr_no;
  for (i=0;i<NCR;i++) {
    if (alist[i]<0 && (dist_cr(i,icur) < 2*rcrowd)) size+=depth(clr_no,i,alist);
  }
  return size+1;
}
/****************************************************************************/
int clustering() {
  int i,alist[NCR],size[NCR],n,nmax,smax;
  
  for (i=0;i<NCR;i++) alist[i]=-1;
  
  n=nmax=smax=0;
  for (i=0;i<NCR;i++) {
    if (alist[i]<0) {
      size[n] = depth(n,i,alist);
      if (size[n]>smax) {smax=size[n]; nmax=n;}
      //   printf("cluster %i   # memb %i\n",n+1,size[n]);  
      n++;
    }
  }

  return smax;
}
/****************************************************************************/
/**************** Function to calculate pair correlation function ***********/
/****************************************************************************/
void pair_correlation_crcr(int iflag,int ind){
  double dxcr, dycr, dzcr, Distcr,r,rho;           
  static double his_pair[NTMP][NBIN],n_pair[NTMP];
  static double binsize, corr, MAXDIST,rho_bulk;     /* corr is a coefficient at final calculation */
  int i, j, k, n;
  FILE *fp;

  if (NCR == 0) return ;

  /***************** Initialization ****************/

  if (iflag<0)
    {
      MAXDIST = sqrt(3.)*boxhf;
      binsize = MAXDIST / NBIN;      /******** Size of bins ********/   
      rho_bulk = (double) NCR/(BOX*BOX*BOX);
      printf(">pair_correlation_crcr: MAXDIST %lf binsize %lf rho_bulk %lf\n",MAXDIST,binsize,rho_bulk);
      for (i=0;i<NTMP;i++)          /******* Loop over all temp. index *********/ 
	{
	  for (j=0;j<NBIN;j++)       /******* Loop over all BINS ***********/
	    {
	      his_pair[i][j]=0;   /* Reset file */
	      n_pair[i] = 0;      /* # visits to ind i */
	    }
	}
      fp=fopen("_his_pair_crcr","w");
      fclose(fp);
      return;
    }
  /******************* Update statistics ******************/
  if (iflag == 0) {
    for (i = 1; i < NCR; i++)         /******* Loop over all crwoders ********/
      {
	dxcr = xcr[0] - xcr[i] ; bc(&dxcr);
	dycr = ycr[0] - ycr[i] ; bc(&dycr);
	dzcr = zcr[0] - zcr[i] ; bc(&dzcr);
	Distcr = dxcr * dxcr + dycr * dycr + dzcr * dzcr;
	n =  sqrt(Distcr) / binsize;                      /****** The bin that distance D falls into *********/
	if (n < 0 || n > NBIN-1) printf("shit %d\n",n);
	his_pair[ind][n]++;                               /****** Increment bin n for any temp.index ********/
      }
    n_pair[ind]++;
    return;
  }
  /******************** Print to file *******************/
  if (iflag>0) {
    fp=fopen("_his_pair_crcr","w");
    for (i = 0; i<NTMP; i++)
      {
	for (k = 0; k < NBIN; k++)
	  {
	    r = (k + 0.5) * binsize;
	    corr = 4 * pi * r * r * binsize;
	    rho = his_pair[i][k] / n_pair[i] / corr ;
	    fprintf(fp,"%i %i %lf %lf %lf %lf\n",i,k,r,his_pair[i][k], rho, rho / rho_bulk);
	  }
      }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void pair_correlation_crat(int iflag,int iatom,int ind){
  double dxcr, dycr, dzcr, Distcr,r,rho;           
  static double his_pair[NTMP][NBIN],n_pair[NTMP];
  static double binsize, corr, MAXDIST,rho_bulk;   
  int i, j, k, n;
  FILE *fp;

  if (NCH == 0) return ;
  
  /***************** Initialization ****************/

  if (iflag<0)
    {
      MAXDIST = sqrt(3.) * boxhf;
      binsize = MAXDIST / NBIN;      /******** Size of bins ********/   
      rho_bulk = (double) NCR/(BOX*BOX*BOX);
      printf(">pair_correlation_crat: MAXDIST %lf binsize %lf rho_bulk %lf\n",MAXDIST,binsize,rho_bulk);
      for (i=0;i<NTMP;i++) {
	for (j=0;j<NBIN;j++) his_pair[i][j] = n_pair[i] = 0;  
      }
      fp=fopen("_his_pair_crat","w");
      fclose(fp);
      return;
    }
  
  /******************* Update statistics ******************/

  if (iflag == 0) {
    for (i = 0; i < NCR; i++)         /******* Loop over all crowders ********/
      {
	dxcr = xcr[i] - x[iatom] ; bc(&dxcr);
	dycr = ycr[i] - y[iatom] ; bc(&dycr);
	dzcr = zcr[i] - z[iatom] ; bc(&dzcr);
	Distcr = dxcr * dxcr + dycr * dycr + dzcr * dzcr;
	n =  sqrt(Distcr) / binsize;                      /****** The bin that distance D falls into *********/
	if (n < 0 || n > NBIN-1) printf("shit %d\n",n);
	his_pair[ind][n]++;                                /****** Increment bin n for any temp.index ********/
      }
    n_pair[ind]++;
    return;
  }
  /******************** Print to file *******************/
  if (iflag>0) {
    fp=fopen("_his_pair_crat","w");
    for (i = 0; i<NTMP; i++)
      {
	for(k = 0; k < NBIN; k++)
	  {
	    r = (k+0.5)*binsize;
	    corr = 4 * pi * r * r * binsize;
	    rho = his_pair[i][k] / n_pair[i] / corr ;
	    fprintf(fp,"%i %i %lf %lf %lf %lf\n",i,k,r,his_pair[i][k], rho, rho / rho_bulk);
	  }
      }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void pair_correlation_cm(int iflag,int ich,int ind){
  double dxcr, dycr, dzcr, Distcr,r,rho;           
  static double his_pair[NTMP][NBIN],n_pair[NTMP];
  static double binsize, corr, MAXDIST,rho_bulk;   
  int i, j, k, n;
  FILE *fp;

  if (NCH == 0) return ;
  
  /***************** Initialization ****************/

  if (iflag<0)
    {
      MAXDIST = sqrt(3.) * boxhf;
      binsize = MAXDIST / NBIN;      /******** Size of bins ********/   
      rho_bulk = (double) NCR/(BOX*BOX*BOX);
      printf(">pair_correlation_crat: MAXDIST %lf binsize %lf rho_bulk %lf\n",MAXDIST,binsize,rho_bulk);
      for (i=0;i<NTMP;i++) {
	for (j=0;j<NBIN;j++) his_pair[i][j] = n_pair[i] = 0;  
      }
      fp=fopen("_his_pair_cm","w");
      fclose(fp);
      return;
    }
  
  /******************* Update statistics ******************/

  double xcm,ycm,zcm;
  
  if (iflag == 0) {
    xcm = ycm = zcm = 0;
    for (i = Nb[ich]; i <= Ne[ich]; i++)
      {
	xcm += x[iCa[i]];
	ycm += y[iCa[i]];
	zcm += z[iCa[i]];
      }
    n = Ne[ich] - Nb[ich] + 1;
    xcm /= n; ycm /= n; zcm /= n;
    
    for (i = 0; i < NCR; i++)         /******* Loop over all crowders ********/
      {
	dxcr = xcr[i] - xcm ; bc(&dxcr);
	dycr = ycr[i] - ycm ; bc(&dycr);
	dzcr = zcr[i] - zcm ; bc(&dzcr);
	Distcr = dxcr * dxcr + dycr * dycr + dzcr * dzcr;
	n =  sqrt(Distcr) / binsize;                      /****** The bin that distance D falls into *********/
	if (n < 0 || n > NBIN-1) printf("shit %d\n",n);
	his_pair[ind][n]++;                                /****** Increment bin n for any temp.index ********/
      }
    n_pair[ind]++;
    return;
  }
  /******************** Print to file *******************/
  if (iflag>0) {
    fp=fopen("_his_pair_cm","w");
    for (i = 0; i<NTMP; i++)
      {
	for(k = 0; k < NBIN; k++)
	  {
	    r = (k+0.5)*binsize;
	    corr = 4 * pi * r * r * binsize;
	    rho = his_pair[i][k] / n_pair[i] / corr ;
	    fprintf(fp,"%i %i %lf %lf %lf %lf\n",i,k,r,his_pair[i][k], rho, rho / rho_bulk);
	  }
      }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
#ifdef ___CALC_STRUCT_FACT
/****************************************************************************/
/************* Calculate Structure Factor for crowders ***********/
/****************************************************************************/
# define NSQ 5000
void hist_sq(int iflag,int ind,int icr) {
  static int qx[NSQ],qy[NSQ],qz[NSQ],nq,q2,qmax = 15,qmax2;
  static double his_sqRe[NSQ],minq,nsq;
  static double hisq[NBIN+1],nhis[NBIN+1],epsq;
  int i,j,k,l,m,n,b;
  double dxcr,dycr,dzcr,q,qdr;
  FILE *fp;

  if (NCR == 0 || ind != NTMP-1) return ;
  
  /***************** Initialization ****************/

  if (iflag < 0) {
    /* included q-vectors: minq * (k,l,m) with magnitude q <= qmax */

    minq = 2 * pi / BOX;   /* Minimum value of q, beacuse of periodic boundary condition */
    qmax2 = qmax * qmax;
    
    for (j = 0; j < NSQ; j++) his_sqRe[j] = 0;
    for (j = 0; j < NBIN+1; j++) hisq[j] = nhis[j] = 0;
    epsq = (double) qmax / NBIN;
    nsq = nq = 0;
    
    /* Generate q vectors */
    
    fp = fopen("sq_vectors","w");
    for (k = 0; k <= qmax; k++)  {        // x-component of q-vector
      for (l = 0; l <= qmax; l++)  {      // y-component of q-vector
	for (m = 0; m <= qmax; m++)  {    // z-component of q-vector
	  if (k == 0 && l == 0 && m == 0) continue;
	  if ((q2 = k * k + l * l + m * m) > qmax2) continue;
	  q = sqrt(q2);
	  qx[nq] = k; qy[nq] = l; qz[nq] = m; nq++;
	  b = q / (double) qmax * NBIN;
	  nhis[b]++;
	  fprintf(fp,"%i %i %i %lf %i\n",k,l,m,minq*q,b);
	}
      }
    }
    fclose(fp); 

    printf("<hist_sq> Number of q vectors nq %i\n",nq);
    printf("<hist_sq> qmax %i i\n",qmax);
    
    if (nq >= NSQ) {
      printf("Too many q vectors\n");
      exit(-1);
    }
    
    /* Output file */
    
    fp = fopen("_his_sq","w");
    fclose(fp);

    return;
  }

  /******************* Update statistics ******************/

  if (iflag == 0) {

    for (i = 0; i < NCR; i++) {
      for (j = 0; j <= i; j++) {

	/*Compute the distance between crowders*/
	dxcr = xcr[i] - xcr[j] ; bc(&dxcr);
	dycr = ycr[i] - ycr[j] ; bc(&dycr);
	dzcr = zcr[i] - zcr[j] ; bc(&dzcr);

	for (n = 0; n < nq; n++) {
	  k = qx[n]; l = qy[n]; m = qz[n]; 
	  qdr = minq * (k * dxcr + l * dycr + m * dzcr);
	  if (j == i)
	    his_sqRe[n] += 1;
	  else
	    his_sqRe[n] += 2 * cos(qdr);
	}
      }
    }
    nsq++;
    
    return;
  }

  /******************** Print to file *******************/

  if (iflag > 0) {
    for (b = 0; b < NBIN+1; b++) hisq[b] = 0;
    
    for (n = 0; n < nq; n++) {
      k = qx[n]; l = qy[n]; m = qz[n];
      q = sqrt(k * k + l * l + m * m);
      b = q / qmax * NBIN;
      hisq[b] += his_sqRe[n] / nsq / NCR;
    }

    fp=fopen("_his_sq","w");
    for (n = 0; n < nq; n++) {
      k = qx[n]; l = qy[n]; m = qz[n];
      q = sqrt(k * k + l * l + m * m);
      fprintf(fp,"%lf %lf %2i %2i %2i %lf %lf\n",
	      Phic, q * minq, k,l,m,his_sqRe[n] / nsq / NCR,
	      pow(1-Phic,4.0)/pow(1+2*Phic,2.0));
    }
    fclose(fp);

    fp=fopen("_his_sq_bin","w");
    for (b = 0; b < NBIN+1; b++) {
      q = (.5 + b) * epsq;
      if (nhis[b] > 0)
	fprintf(fp,"%lf %lf %lf %i %lf\n",
		Phic, q * minq, hisq[b] / nhis[b],
		(int) nhis[b],pow(1-Phic,4.0)/pow(1+2*Phic,2.0));
    }
    fclose(fp);

    return;
  }
}
#undef NSQ
#endif
/****************************************************************************/
/*Function to calculate the radius of gyration from the coorddinate( in the case of native state) :*/
/******************************************************************************/
double rgyr_from_coord(double *xt,double *yt,double *zt,int a1,int a2)
{
  int i,n;
  double xcm,ycm,zcm,rg2;

  rg2 = n = xcm = ycm = zcm = 0;
  //Calculate the center of mass
  for (i=a1;i<=a2;i++){
    xcm += xt[i];
    ycm += yt[i];
    zcm += zt[i];
    n++;
  }
  xcm /= n; ycm /= n; zcm /= n;
  
  //Calculate the radius of gyration
  
  for (i=a1;i<=a2;i++){
    rg2 += ((xt[i]-xcm)*(xt[i]-xcm) +
	    (yt[i]-ycm)*(yt[i]-ycm) +
	    (zt[i]-zcm)*(zt[i]-zcm));
  }

  return sqrt(rg2/n);
}
/****************************************************************************/
double calcrmsd_ca(int iflag,double *xt,double *yt,double *zt,int a1,int a2)
{
  int i=0,j;
  double x1[N],y1[N],z1[N];
  double x2[N],y2[N],z2[N];

  if (RMSD==0 || NCH == 0) return 0;

  for (j=a1;j<=a2;j++) {
    x1[i]=x[iCa[j]];
    y1[i]=y[iCa[j]];
    z1[i]=z[iCa[j]];

    x2[i]=xt[i];
    y2[i]=yt[i];
    z2[i]=zt[i];

    i++;
  }


  if (iflag==0) 
    return sqrt(rmsd2(x1,y1,z1,x2,y2,z2,i));
  else 
    return sqrt(rmsd2_nopt(x1,y1,z1,x2,y2,z2,i));
}
/****************************************************************************/
/*double interhp() 
{
  int i,j;
  double e=0;

  if (NCH==1) return 0;

  for (i=0;i<NCH;i++) for (j=i+1;j<NCH;j++) e+=Mhp[i][j];
  return epshp*e;
  }*/
/****************************************************************************/
/*double interhb() 
{
  int i,j;
  double e=0;

  if (NCH==1) return 0;

  for (i=0;i<NCH;i++) for (j=i+1;j<NCH;j++) e+=Mhb[i][j];
  return epshb*e;
  }*/
/****************************************************************************/
int dist2(int x1,int y1,int z1,int x2,int y2,int z2) {
  return ( (x1-x2)*(x1-x2) + (y1-y2)*(y1-y2) + (z1-z2)*(z1-z2) );
}
/****************************************************************************/
int cont_made(int iflag,int i, int j)
/* returns 1 if residues i,j are in contact, and 0 otherwise */
{
  int ica,icb,jca,jcb;
  static int rcut2;
  
  if (iflag<0) {
    printf("cont_made init: ");
    rcut2 = 7.0 * 7.0;
    printf(" rcut = %lf\n",sqrt(rcut2));
    return 0;
  }

  ica = iCa[i];
  icb = iCa[i] + 2;
  jca = iCa[j];
  jcb = iCa[j] + 2;

  if ( dist2(x[ica],y[ica],z[ica],x[jca],y[jca],z[jca]) < rcut2 ) return 1;
  if ( dist2(x[ica],y[ica],z[ica],x[jcb],y[jcb],z[jcb]) < rcut2 ) return 1; 
  if ( dist2(x[icb],y[icb],z[icb],x[jca],y[jca],z[jca]) < rcut2 ) return 1;
  if ( dist2(x[icb],y[icb],z[icb],x[jcb],y[jcb],z[jcb]) < rcut2 ) return 1; 
  
  return 0;
}
/****************************************************************************/
int contacts()
/* returns # native contacts formed */
{
  int i, j, k, n_counted;
	
  for (k=n_counted=0; k<ncontacts; k++) {
     i = contact_pair_i[k];
     j = contact_pair_j[k];
     if ( cont_made(0,i,j) ) n_counted++;
  }

  return n_counted;
}
/****************************************************************************/
int contacts_in_box(int b) 
/* returns # contacts formed in box b */
{
  int i, j, k, n_counted;
	
  for (k=n_counted=0; k<ncontacts; k++) {
    if (boxnum[k] != b) continue;

    i = contact_pair_i[k];
    j = contact_pair_j[k];
    if ( cont_made(0,i,j) ) n_counted++;
  }
  
  return n_counted;
}
/****************************************************************************/
int contacts_all()
/* number of contacts between any two amino acids ij */
{
  int i, j, n_counted = 0;

  for (i=0; i<N; i++) {
    if ( seq[i]=='G' ) continue;
    for (j=0; j<i-2; j++) {
      if ( seq[j]=='G' ) continue;
      if ( cont_made(0,i,j) ) n_counted++;
    }
  }

  return n_counted;
}
/****************************************************************************/
int contacts_res(int iaa)
/* number of contacts between iaa and any other amino acids */
{
  int j, n_counted = 0;

  for (j=0; j<N; j++) {
      if ( seq[iaa] == 'G' || seq[j]== 'G' ) continue;
      if ( j >= iaa-2 && j <= iaa+2 ) continue;
      if ( cont_made(0,iaa,j) ) n_counted++;
    }

  return n_counted;
}
/****************************************************************************/
int contact_profile(int ind,double cont_vec[][N],double count[]) {
  int i,j,k;
  
  for (k=0; k<ncontacts; k++) {
    i = contact_pair_i[k];
    j = contact_pair_j[k];
    if ( cont_made(0,i,j) ) {
      cont_vec[ind][i]++;
      cont_vec[ind][j]++;
    }
  }
  
  count[ind]++;
  
  return 0;
  
} 
/****************************************************************************/
void write_phival(double vec1[][N],double count1[],double vec2[][N],double count2[]) {
  int i,j;
  double Ncont_ave_tse,Ncont_ave_nat,Ncont_ave_nat_Tlow;
  FILE *fp;

  fp = fopen("_phival","w");

  for (i=0;i<NTMP;i++) {
    for (j=0;j<N;j++) {
      if (seq[j] == 'G') continue;

      Ncont_ave_tse = vec1[i][j]/count1[i];
      Ncont_ave_nat = vec2[i][j]/count2[i];
      Ncont_ave_nat_Tlow = vec2[NTMP-1][j]/count2[NTMP-1];

      fprintf(fp,"%i %i %lf %lf %lf %lf\n",i,j,Ncont_ave_tse,Ncont_ave_nat, \
	      Ncont_ave_tse/Ncont_ave_nat,Ncont_ave_tse/Ncont_ave_nat_Tlow);
    }
  }

  fclose(fp);

  return ;
} 
/****************************************************************************/
double relative_contact_order() {
  int i,j;
  double rco=0,n=0;

  //  for (k=n_counted=0; k<ncontacts; k++) {
  //    i = contact_pair_i[k];
  //    j = contact_pair_j[k];

  for (i=0;i<N;i++) {
    for (j=0;j<i-3;j++) {

      if ( cont_made(0,i,j) ) {
	rco += (i-j);
	n++;
      }

    }
  }
  
  return rco/(n*N);
}
/****************************************************************************/
void rco_vs_ncont(int iflag,double rco,int ncont) {
  static double his[MAXCONT+1],n[MAXCONT+1];
  int i;
  FILE *fp;

  if (iflag < 0) {
    for (i=0;i<=MAXCONT;i++) his[i] = n[i] = 0;
    fp=fopen("_rco_vs_ncont","w"); /* reset file */
    fclose(fp);
    return;
  }

  if (iflag == 0) {
    his[ncont] += rco;
    n[ncont]++;
    return;
  }
  
  if (iflag > 0) {
    fp = fopen("_rco_vs_ncont","w");
    for (i=0;i<=MAXCONT;i++) {
      if (n[i] > 0) fprintf(fp,"%i %lf\n",i,his[i]/n[i]);
    }
    fclose(fp);
  }
  
  return;
}
/****************************************************************************/
void his_box_ncont(int iflag,int ncont) {
  static double his[11][MAXCONT+1],n[MAXCONT+1];
  int i,j;
  FILE *fp;

  if (iflag < 0) {
    for (j=0; j<11; j++) for (i=0;i<=MAXCONT;i++) his[j][i] = n[i] = 0;
    fp=fopen("_his_box_cont","w"); /* reset file */
    fclose(fp);
    return;
  }

  if (iflag == 0) {
    for (j=0; j<11; j++) {
      his[j][ncont] += contacts_in_box(j); 
    }
    n[ncont]++;
    return;
  }
  
  if (iflag > 0) {
    fp = fopen("_his_box_cont","w");

    for (j=0; j<11; j++) {
      for (i=0; i<=MAXCONT; i++) {
	if (n[i] > 0 && his[j][i] > 0) {
	  fprintf(fp,"%i %i %lf %lf %lf\n",j,i,\
		  his[j][i]/n[i],his[j][i],n[i]);
	}
      }
    }

    fclose(fp);
  }

  if (iflag > 1) {
    int nat1,nat2;
    int tse1,tse2;
    double nc,sum;
    
    fp = fopen("tse_nat.def","r");
    fscanf(fp,"%d %d %d %d",&tse1,&tse2,&nat1,&nat2);
    fclose(fp);

    fp = fopen("_his_box_cont_nat2","w");

    for (j=0; j<11; j++) {
      nc = sum = 0;
      for (i=nat2; i<=MAXCONT; i++) {
	nc += his[j][i];
	sum += n[i];
      }
      if (nc > 0 && sum > 0 ) {
	fprintf(fp,"%i %lf %lf %lf\n",j,nc/sum,nc,sum);
      }
    }

    fclose(fp);
  }

  return;
}
/****************************************************************************/
void histo_cont(int iflag,int ind,int n) { 
  static double his[NTMP][MAXCONT+1];
  static double ene[NTMP][MAXCONT+1];
  static double sum;
  int i,j;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) 
      for (j=0;j<=MAXCONT;j++) his[i][j]=ene[i][j]=0;
    fp=fopen("_his_cont","w"); /* reset file */
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (n<0 || n>MAXCONT) printf("shit %d\n",n);
    his[ind][n]++;
    ene[ind][n]+=E;
    return;
  }
  
  if (iflag>0) {
    fp=fopen("_his_cont","w");
    for (i=0;i<NTMP;i++) {
      for (j=sum=0;j<=MAXCONT;j++) sum+=his[i][j];
      for (j=0;j<=MAXCONT;j++) {
	if (his[i][j]>0) fprintf(fp,"%i %i %f %f\n",i,j,his[i][j]/sum,ene[i][j]/his[i][j]);
      }
      fprintf(fp,"\n");
    }
    fclose(fp);
  }
  
  return;
}
/****************************************************************************/
void histo_contmap(int iflag, int ind)
{ 
  static double his[NTMP][N][N],n[NTMP];
  int k, i, j;
  FILE *fp;

  if (iflag<0) {
    for (k=0;k<NTMP;k++) 
      for (i=0;i<N;i++) for (j=0;j<N;j++) his[k][i][j]=n[k]=0;
    fp=fopen("_his_contmap","w");
    fclose(fp);
    return;
  }

  if (iflag==0) { 
    for (i=0;i<N;i++) {
      for (j=0;j<i;j++) {
	if (cont_made(0,i,j) > 0) {
	  his[ind][i][j]++;	
	  his[ind][j][i]++;
	}
      }
    }
    n[ind]++;	
    return;		
  }	

   if (iflag>0) {
    fp=fopen("_his_contmap","w");
    for (k=0;k<NTMP;k++) {
      if (n[k]>0) {
	for (i=0;i<N;i++) {
	  for (j=0;j<N;j++) {
            fprintf(fp,"%i %i %i %f\n", k, i, j, his[k][i][j]/n[k]);
          }
	}	  
	fprintf(fp,"\n");
      }
    }
    fclose(fp);
    return;
  }	
}
/****************************************************************************/
/*void histo_cont_rg(int iflag,int itmp,int cont,double x)
{ 
  static double low2=0.0,high2=35.0;
  static double his[NTMP][MAXCONT+1][NBIN],out[NTMP];
  static double eps2,sum;
  int i,j,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) 
      for (k=0;k<MAXCONT+1;k++) for (j=0;j<NBIN;j++) his[i][k][j]=out[i]=0;
    eps2=(high2-low2)/NBIN;
    fp=fopen("_his_cont_rg","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (cont>=0 && cont<=MAXCONT && x>low2 && x<high2) { 
      j=(x-low2)/eps2; 
      his[itmp][cont][j]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_his_cont_rg","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<MAXCONT+1;k++) for (j=0;j<NBIN;j++) sum+=his[i][k][j];
      if (sum>0) {
	for (k=0;k<MAXCONT+1;k++) {
	  for (j=0;j<NBIN;j++) {
	    fprintf(fp,"%i %i %f %f\n",i,k,
		    low2+eps2*(j+.5),his[i][k][j]/sum/eps2);
	  }
	}	  
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i cont_rg out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
  } */
/****************************************************************************/
void helixprof(int iflag,int itmp)
{ 
  static double his[NTMP][N],n[NTMP];
  int i,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) for (k=0;k<N;k++) his[i][k]=n[i]=0;
    fp=fopen("_helixprof","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    for (i=0;i<N;i++) his[itmp][i]+=helix_aa(i);
    n[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_helixprof","w");
    for (i=0;i<NTMP;i++) {
      for (k=0;k<N && n[i]>0;k++) 
	fprintf(fp,"%i %i %e\n",i,k,his[i][k]/n[i]);
      fprintf(fp,"\n");
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histoe(int iflag,int itmp,double x)
{ 
  static double low=-80.0;      
  static double high=120.0;
  static double his[NTMP][NBIN],out[NTMP];
  static double eps,sum;
  int i,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) for (k=0;k<NBIN;k++) his[i][k]=out[i]=0;
    eps=(high-low)/NBIN;
    fp=fopen("_hise","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x>low && x<high) { 
      k=(x-low)/eps; his[itmp][k]++;
    }
    else out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hise","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBIN;k++) sum+=his[i][k];
      if (sum>0) {
	for (k=0;k<NBIN;k++) 
	  fprintf(fp,"%i %f %f\n",i,low+eps*(k+.5),his[i][k]/sum/eps);
      }
      fprintf(fp,"\n");
      if (iflag>1 && out[i]>0) printf("temp %i energy out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histod1(int iflag,int itmp,double x)
{ 
  static double low=0.0;      
  static double high=30.0;
  static double his[NTMP][NBINL],out[NTMP];
  static double eps,sum;
  int i,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) for (k=0;k<NBINL;k++) his[i][k]=out[i]=0;
    eps=(high-low)/NBINL;
    fp=fopen("_hisd1","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x>low && x<high) { 
      k=(x-low)/eps; his[itmp][k]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hisd1","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBINL;k++) sum+=his[i][k];
      if (sum>0) {
	for (k=0;k<NBINL;k++) 
	  fprintf(fp,"%i %lf %lf\n",i,low+eps*(k+.5),his[i][k]/sum/eps);
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i rmsd out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histod2(int iflag,int itmp,double x)
{ 
  static double low=0.0;      
  static double high=30.0;
  static double his[NTMP][NBINL],out[NTMP];
  static double eps,sum;
  int i,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) for (k=0;k<NBINL;k++) his[i][k]=out[i]=0;
    eps=(high-low)/NBINL;
    fp=fopen("_hisd2","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x>low && x<high) { 
      k=(x-low)/eps; his[itmp][k]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hisd2","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBINL;k++) sum+=his[i][k];
      if (sum>0) {
	for (k=0;k<NBINL;k++) 
	  fprintf(fp,"%i %lf %lf\n",i,low+eps*(k+.5),his[i][k]/sum/eps);
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i rmsd out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histod3(int iflag,int itmp,double x)
{ 
  static double low=0.0;      
  static double high=30.0;
  static double his[NTMP][NBINL],out[NTMP];
  static double eps,sum;
  int i,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) for (k=0;k<NBINL;k++) his[i][k]=out[i]=0;
    eps=(high-low)/NBINL;
    fp=fopen("_hisd3","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x>low && x<high) { 
      k=(x-low)/eps; his[itmp][k]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hisd3","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBINL;k++) sum+=his[i][k];
      if (sum>0) {
	for (k=0;k<NBINL;k++) 
	  fprintf(fp,"%i %lf %lf\n",i,low+eps*(k+.5),his[i][k]/sum/eps);
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i rmsd out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histod4(int iflag,int itmp,double x)
{ 
  static double low=0.0;      
  static double high=30.0;
  static double his[NTMP][NBINL],out[NTMP];
  static double eps,sum;
  int i,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) for (k=0;k<NBINL;k++) his[i][k]=out[i]=0;
    eps=(high-low)/NBINL;
    fp=fopen("_hisd4","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x>low && x<high) { 
      k=(x-low)/eps; his[itmp][k]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hisd4","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBINL;k++) sum+=his[i][k];
      if (sum>0) {
	for (k=0;k<NBINL;k++) 
	  fprintf(fp,"%i %lf %lf\n",i,low+eps*(k+.5),his[i][k]/sum/eps);
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i rmsd out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histoed1(int iflag,int itmp,double x1,double x2)
{ 
  static double low1=-100.0,low2=0.0;      
  static double high1=100.0,high2=30.0;
  static double his[NTMP][NBIN][NBIN],out[NTMP];
  static double eps1,eps2,sum;
  int i,j,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) 
      for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) his[i][k][j]=out[i]=0;
    eps1=(high1-low1)/NBIN;
    eps2=(high2-low2)/NBIN;
    fp=fopen("_hised1","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x1>low1 && x1<high1 && x2>low2 && x2<high2) { 
      k=(x1-low1)/eps1; j=(x2-low2)/eps2; 
      his[itmp][k][j]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hised1","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) sum+=his[i][j][k];
      if (sum>0) {
	for (k=0;k<NBIN;k++) {
	  for (j=0;j<NBIN;j++) {
	    fprintf(fp,"%i %f %f %f\n",i,low1+eps1*(k+.5),
		    low2+eps2*(j+.5),his[i][k][j]/sum/eps1/eps2);
	  }
	}	  
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i ermsd out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histoed2(int iflag,int itmp,double x1,double x2)
{ 
  static double low1=-100.0,low2=0.0;      
  static double high1=100.0,high2=30.0;
  static double his[NTMP][NBIN][NBIN],out[NTMP];
  static double eps1,eps2,sum;
  int i,j,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) 
      for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) his[i][k][j]=out[i]=0;
    eps1=(high1-low1)/NBIN;
    eps2=(high2-low2)/NBIN;
    fp=fopen("_hised2","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x1>low1 && x1<high1 && x2>low2 && x2<high2) { 
      k=(x1-low1)/eps1; j=(x2-low2)/eps2; 
      his[itmp][k][j]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hised2","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) sum+=his[i][j][k];
      if (sum>0) {
	for (k=0;k<NBIN;k++) {
	  for (j=0;j<NBIN;j++) {
	    fprintf(fp,"%i %f %f %f\n",i,low1+eps1*(k+.5),
		    low2+eps2*(j+.5),his[i][k][j]/sum/eps1/eps2);
	  }
	}	  
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i ermsd out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histoed3(int iflag,int itmp,double x1,double x2)
{ 
  static double low1=0.0,low2=0.0;      
  static double high1=30.0,high2=30.0;
  static double his[NTMP][NBIN][NBIN],out[NTMP];
  static double eps1,eps2,sum;
  int i,j,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) 
      for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) his[i][k][j]=out[i]=0;
    eps1=(high1-low1)/NBIN;
    eps2=(high2-low2)/NBIN;
    fp=fopen("_hised3","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x1>low1 && x1<high1 && x2>low2 && x2<high2) { 
      k=(x1-low1)/eps1; j=(x2-low2)/eps2; 
      his[itmp][k][j]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hised3","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) sum+=his[i][j][k];
      if (sum>0) {
	for (k=0;k<NBIN;k++) {
	  for (j=0;j<NBIN;j++) {
	    fprintf(fp,"%i %f %f %f\n",i,low1+eps1*(k+.5),
		    low2+eps2*(j+.5),his[i][k][j]/sum/eps1/eps2);
	  }
	}	  
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i ermsd out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histodd1(int iflag,int itmp,double x1,double x2)
{ 
  static double low1=0.0,low2=0.0;      
  static double high1=80.0,high2=30.0;
  static double his[NTMP][NBIN][NBIN],out[NTMP];
  static double eps1,eps2,sum;
  int i,j,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) 
      for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) his[i][k][j]=out[i]=0;
    eps1=(high1-low1)/NBIN;
    eps2=(high2-low2)/NBIN;
    fp=fopen("_hisdd1","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x1>low1 && x1<high1 && x2>low2 && x2<high2) { 
      k=(x1-low1)/eps1; j=(x2-low2)/eps2; 
      his[itmp][k][j]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hisdd1","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) sum+=his[i][j][k];
      if (sum>0) {
	for (k=0;k<NBIN;k++) {
	  for (j=0;j<NBIN;j++) {
	    fprintf(fp,"%i %f %f %f\n",i,low1+eps1*(k+.5),
		    low2+eps2*(j+.5),his[i][k][j]/sum/eps1/eps2);
	  }
	}	  
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i dd1 out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
} 
/****************************************************************************/
void histodd2(int iflag,int itmp,double x1,double x2)
{ 
  static double low1=0.0,low2=0.0;      
  static double high1=80.0,high2=30.0;
  static double his[NTMP][NBIN][NBIN],out[NTMP];
  static double eps1,eps2,sum;
  int i,j,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) 
      for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) his[i][k][j]=out[i]=0;
    eps1=(high1-low1)/NBIN;
    eps2=(high2-low2)/NBIN;
    fp=fopen("_hisdd2","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (x1>low1 && x1<high1 && x2>low2 && x2<high2) { 
      k=(x1-low1)/eps1; j=(x2-low2)/eps2; 
      his[itmp][k][j]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hisdd2","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<NBIN;k++) for (j=0;j<NBIN;j++) sum+=his[i][j][k];
      if (sum>0) {
	for (k=0;k<NBIN;k++) {
	  for (j=0;j<NBIN;j++) {
	    fprintf(fp,"%i %f %f %f\n",i,low1+eps1*(k+.5),
		    low2+eps2*(j+.5),his[i][k][j]/sum/eps1/eps2);
	  }
	}	  
	fprintf(fp,"\n");
      }
      if (iflag>1 && out[i]>0) printf("temp %i dd2 out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
void histonn1(int iflag,int itmp,int m,int n)
{ 
  static int low1=0,low2=0;      
  static int high1=100,high2=100;
  static double his[NTMP][101][101],out[NTMP];
  static double sum;
  int i,j,k;
  FILE *fp;

  if (iflag<0) {
    for (i=0;i<NTMP;i++) 
      for (k=0;k<101;k++) for (j=0;j<101;j++) his[i][k][j]=out[i]=0;
    fp=fopen("_hisnn1","w");
    fclose(fp);
    return;
  }

  if (iflag==0) {
    if (m>=low1 && m<=high1 && n>=low2 && n<=high2) { 
      his[itmp][m][n]++;
    } else 
      out[itmp]++;
    return;
  }

  if (iflag>0) {
    fp=fopen("_hisnn1","w");
    for (i=0;i<NTMP;i++) {
      sum=0; for (k=0;k<101;k++) for (j=0;j<101;j++) sum+=his[i][j][k];
      if (sum>0) {
	for (k=0;k<101;k++) {
	  for (j=0;j<101;j++) {
	    fprintf(fp,"%i %i %i %f\n",i,k,j,his[i][k][j]/sum);
	  }
	  fprintf(fp,"\n");
	}	  
      }
      if (iflag>1 && out[i]>0) printf("temp %i nn1 out %f\n",i,out[i]);
    }
    fclose(fp);
    return;
  }
}
/****************************************************************************/
int sheet(int ic)
{
  int i,j,ns=0;
  double phi,psi;
  
  for (j=0;j<NCH;j++) {
    if (j==ic || ic<0) {
      for (i=Nb[j];i<=Ne[j];i++) {
	phi=ph[3*i]; psi=ph[3*i+1];
	while (phi>pi) phi-=pi2; while (phi<-pi) phi+=pi2;
	while (psi>pi) psi-=pi2; while (psi<-pi) psi+=pi2;
	phi*=180/pi; psi*=180/pi;
	if (phi>-160 && phi<-50 && psi>100 && psi<160) ns++;
      }
    }
  }

  return ns;
}
/****************************************************************************/
int helix(int ic)
{
  int i,j,nh=0;
  
  for (j=0;j<NCH;j++) {
    if (j==ic || ic<0) {      
      for (i=Nb[j];i<=Ne[j];i++) {
	nh+=helix_aa(i);
      }
    }
  }

  return nh;
}
/****************************************************************************/
int helix_aa(int i)
{
  double phi,psi;
  
  phi=ph[3*i]; 
  psi=ph[3*i+1];

  while (phi>pi) phi-=pi2; while (phi<-pi) phi+=pi2;
  while (psi>pi) psi-=pi2; while (psi<-pi) psi+=pi2;

  phi*=180/pi; 
  psi*=180/pi;

  if (phi>-90 && phi<-30 && psi>-77 && psi<-17) return 1;

  return 0;
}
/****************************************************************************/
double cmdist(int iflag,int ic1,int ic2) {
  double xcm1,ycm1,zcm1;
  double xcm2,ycm2,zcm2;
  double xt,yt,zt;
  double d2,tmp;
  int i,j,k;
  int i0,j0,k0;

  static int noff1,noff2;

  if (ic1>NCH-1 || ic2>NCH-1) return 0;

  if (iflag<0) {
    noff1=0;
    noff2=0;
    //    noff2=(Ne[1]-Nb[1]+1)/2-8;
    printf("CM calculation: noff1 %i noff2 %i\n",noff1,noff2); 
    return 0;
  }

  cm(&xcm1,&ycm1,&zcm1,ic1,noff1);
  cm(&xcm2,&ycm2,&zcm2,ic2,noff2);

  d2=((xcm1-xcm2)*(xcm1-xcm2)+
      (ycm1-ycm2)*(ycm1-ycm2)+
      (zcm1-zcm2)*(zcm1-zcm2)); 

  i0=j0=k0=0;

  for (i=-1;i<=1;i++) {
    for (j=-1;j<=1;j++) {
      for (k=-1;k<=1;k++) {
	if (i==0 && j==0 && k==0) continue;
	xt=xcm2+i*BOX;
	yt=ycm2+j*BOX;
	zt=zcm2+k*BOX;
	if ((tmp=((xcm1-xt)*(xcm1-xt)+
		  (ycm1-yt)*(ycm1-yt)+
		  (zcm1-zt)*(zcm1-zt)))<d2) 
	  {
	    d2=tmp;
	    i0=i;
	    j0=j;
	    k0=k;
	  }
      }
    }
  }
  
  if (iflag>0) {

    for (i=iBeg[ic2];i<=iEnd[ic2];i++) {
      x[i]=x[i]+i0*BOX;
      y[i]=y[i]+j0*BOX;
      z[i]=z[i]+k0*BOX;
      xo[i]=x[i];
      yo[i]=y[i];
      zo[i]=z[i];
    }

  }

  return sqrt(d2);
}
/****************************************************************************/
void rama() {
  double phi,psi;
  int i,j;
  FILE *fp;
  
  fp = fopen("_rama","a");
  
  for (i=0;i<NCH;i++) {
    for (j=Nb[i]+1;j<Ne[i];j++) {
	phi=ph[3*j]; psi=ph[3*j+1];
	while (phi>pi) phi-=pi2; while (phi<-pi) phi+=pi2;
	while (psi>pi) psi-=pi2; while (psi<-pi) psi+=pi2;
	phi*=180/pi; psi*=180/pi;
	
    }
  }
  
  fclose(fp);
}
/****************************************************************************/
