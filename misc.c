# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include <sys/types.h>
# include <sys/stat.h>
# include "mt64.h"
# include "defs.h"
# include "sys.h"
# include "global.h"
/****************************************************************************/
void write_xyz(double *x,double *y,double *z,FILE *fp);
void write_ph(double *ph1,double *ph2,FILE *fp);
void read_xyz(double *x,double *y,double *z,FILE *fp);
void read_ph(double *ph1,double *ph2,FILE *fp);
/****************************************************************************/
/***** INPUT/OUTPUT *********************************************************/
/****************************************************************************/
//double settings_read(char str_target[]) {
  //FILE *fp;
  //char str[100];
  //double val;

  //fp = fopen(SETTINGS,"r");

  //if (fp == NULL) return 0;
  
  //while (feof(fp) == 0) {
    //strcpy(str,"");
    //fscanf(fp,"%s",str);
    //if (strcmp(str,str_target) == 0) {
      //fscanf(fp,"%lf",&val);
      //printf("settings file> assigning %s=%lf\n",str_target,val);
      //}
    //while (fgetc(fp) != '\n' && feof(fp) == 0);
    //}
  
  //fclose(fp);

  //return val;
  //}
/****************************************************************************/
void printheading(double frcrd[],double frpiv[],double frbgs[],double frtrl[])
{
  int i;
  printf("\nProgram chp.c (050414)\n");
  printf("\n");
  printf("ind, T, g, fraction crd (of total), fraction piv, fraction bgs, fraction trl:\n");
  for (i=0;i<NTMP;i++) {
    printf("%i %10.6f %10.6f %10.6f %10.6f %10.6f %10.6f\n",i,1./beta[i],g[i],
	   frcrd[i],frpiv[i],frbgs[i],frtrl[i]);
  }
  printf("\n");
  printf("Bond lengths: NCa %f  CaC %f  CN %f\n",b[1],b[2],b[0]);
  printf("              NH %f  CaHa %f  CaCb %f  CO %f\n",bNH,bCaHa,bCaCb,bCO);
  printf("H bonds: powa %f  powb %f  epshb %f\n",
	 POWA,POWB,epshb);
  printf("         sighb %f  cuthb %f\n",
          sighb,cuthb);
  printf("Hydrophobicity: epshp %f sighp %f\n",epshp,sighp);
  printf("Bias: kbias %f\n",kbias);
  printf("Self-avoidance: ksa %f  cut %f\n",ksa,cut);
  printf("Radii: C %f  N %f  O %f  H %f\n",
	  sigsa[0],sigsa[1],sigsa[2],sigsa[3]);
  printf("       ALPHA1 %f ALPHA2 %f\n",ALPHA1,ALPHA2);
  printf("Native structure from files: \n");  
  printf("  %s\n",NATIVE1);  
  printf("  %s\n",NATIVE2);  
  printf("  RMSD %i NR1 %i NR2 %i\n",RMSD,NR1,NR2);  
  printf("MC parameters:\n");
  printf("  N %i NCH %i NCR %i BOX %f\n",N,NCH,NCR,(double)BOX);
  printf("  MCCYC %i  NTHERM %i  MCCON %i\n",
	 MCCYC,NTHERM,MCCON);
  printf("  PIVSTP %f  TRLSTP %f  RGRSTP %f\n",PIVSTP,TRLSTP,RGRSTP);
  printf("  abgs %f  bbgs %f\n",abgs,bbgs);
  printf("  INORM %i  IRT %i  IRAMA %i  ICONF %i\n",
          INORM,IRT,IRAMA,ICONF);
  printf("\n-----------------------------\n\n");
  fflush(0);
} 
/****************************************************************************/
void write_checkpnt(long *imc,double so[][NOBS],double mino[]) 
{ 
  int i,j;
  double sw[NTMP*NOBS];
  FILE *fp;

  int2cart_all(1);
  for (i=0;i<NTO;i++) {xo[i]=x[i]; yo[i]=y[i]; zo[i]=z[i];}
  ecalc();

  for (i=0;i<NTMP;i++) for (j=0;j<NOBS;j++) sw[NOBS*i+j]=so[i][j];

  fp=fopen(CHECKDATA,"w");
  fwrite(imc,sizeof(long),1,fp);
  fwrite(&E,sizeof(double),1,fp);
  fwrite(&ind,sizeof(int),1,fp);
  fwrite(sw,sizeof(double),NTMP*NOBS,fp);
  fwrite(mino,sizeof(double),NOBS,fp);
  if (RAN_RAN3) {
    fwrite(&seed,sizeof(long),1,fp);
    fwrite(&orig_seed,sizeof(long),1,fp);
  }
  fclose(fp); 

  if (RAN_MARSENNE_TWISTER)
    save_twister("CHECKRAND");

  write_conf_ver2(0,CHECKCONF,"w");
  
  /*  historot(3,0);
      histoe(3,0,0);
      histohb(3,0,0);
      histoq(3,0,0);
      histormsd(3,0,0);
      histormsd(3,0,0);
      helix_prof(3,0); */

  return;
}
/****************************************************************************/
void read_checkpnt(long *imc,double so[][NOBS],double mino[]) 
{ 
  int i,j;
  double sr[NTMP*NOBS],e;
  long orig;
  FILE *fp;
  
  fp=fopen(CHECKDATA,"r");

  fread(imc,sizeof(long),1,fp);
  fread(&e,sizeof(double),1,fp);
  fread(&ind,sizeof(int),1,fp);
  fread(sr,sizeof(double),NTMP*NOBS,fp);
  fread(mino,sizeof(double),NOBS,fp);
  if (RAN_RAN3) {
    fread(&seed,sizeof(long),1,fp);
    fread(&orig,sizeof(long),1,fp);
  }
  fclose(fp); 

  for (i=0;i<NTMP;i++) for (j=0;j<NOBS;j++) so[i][j]=sr[NOBS*i+j];

  (*imc)++;
  if (RAN_RAN3) {
    orig_seed = orig;
    while (orig<seed) ran3n(&orig);
  } else if (RAN_MARSENNE_TWISTER) {
    restore_twister("CHECKRAND");
  }

  read_conf_ver2(-1,CHECKCONF,NCH);
  if (0 == read_conf_ver2(0,CHECKCONF,NCH)) {printf("Exiting...\n"); exit(-1);}
  read_conf_ver2(1,CHECKCONF,NCH);

  int2cart_all(1);
  for (i=0;i<NTO;i++) {xo[i] = x[i]; yo[i] = y[i]; zo[i] = z[i];}
  for (i=0;i<NCR;i++) {xcro[i] = xcr[i]; ycro[i] = ycr[i]; zcro[i] = zcr[i];} 
  ecalc();
  
  printf("\n************************************");
  printf("\n**** Restarting from checkpoint ****");
  printf("\n************************************");
  printf("\nReading checkpoint file %s",CHECKDATA);
  printf("\nReading checkpoint file %s",CHECKCONF);
  printf("\nimc %li ind %i seed %li orig_seed %li",*imc,ind,seed,orig_seed);
  printf("\nEnergy: calc %e read %e diff %e",E,e,E-e);
  printf("\nRestarting from MC cycle imc %li",*imc);
  printf("\n************************************\n\n");
  
  /*  historot(-2,0);
      histoe(-2,0,0);
      histohb(-2,0,0);
      histoq(-2,0,0);
      histormsd(-2,0,0);
      histormsd(-2,0,0);
      helix_prof(-2,0); */

  fflush(stdout);

  return;
}
/****************************************************************************/
void write_xyz(double *x,double *y,double *z,FILE *fp) {
  fwrite(x,sizeof(double),1,fp);
  fwrite(y,sizeof(double),1,fp);
  fwrite(z,sizeof(double),1,fp);
}
/****************************************************************************/
void write_ph(double *ph1,double *ph2,FILE *fp) {
  fwrite(ph1,sizeof(double),1,fp);
  fwrite(ph2,sizeof(double),1,fp);
}
/****************************************************************************/
void read_xyz(double *x,double *y,double *z,FILE *fp) {
  fread(x,sizeof(double),1,fp);
  fread(y,sizeof(double),1,fp);
  fread(z,sizeof(double),1,fp);
}
/****************************************************************************/
void read_ph(double *ph1,double *ph2,FILE *fp) {
  fread(ph1,sizeof(double),1,fp);
  fread(ph2,sizeof(double),1,fp);
}
/****************************************************************************/
void write_conf_ver2(int iflag, char *fn,char *fmode) 
{
  int i,j,nch=NCH,ncr=NCR,naa;
  static FILE *fps = NULL;
  FILE *fpw = NULL;
  char fnw[100];
  
  if (iflag < 0) {
    sprintf(fnw,"%s_ver2",fn);
    fps = fopen(fnw,fmode);
    return ;
  }
  
  if (iflag > 0) {
    if (fps != NULL)
      fclose(fps);
    fps = NULL;
    return ;
  }

  if (iflag == 0) {
    fpw = fps;

    if (fps == NULL) {
      sprintf(fnw,"%s_ver2",fn);
      fpw = fopen(fn,fmode);
    }
    
    fwrite(&ind,sizeof(int),1,fpw);
    fwrite(&nch,sizeof(int),1,fpw);
    fwrite(&ncr,sizeof(int),1,fpw);
    fwrite(&E,sizeof(double),1,fpw);
    
    for (i=0; i<NCH; i++) {
      naa = Ne[i] - Nb[i] + 1; 
      fwrite(&naa,sizeof(int),1,fpw);
      j = iN[Nb[i]];     write_xyz(x + j,y + j,z + j,fpw);
      j = iN[Nb[i]] + 1; write_xyz(x + j,y + j,z + j,fpw);
      j = iCa[Nb[i]];    write_xyz(x + j,y + j,z + j,fpw);
      for (j=Nb[i]; j<=Ne[i]; j++) write_ph(ph + 3*j,ph + 3*j + 1,fpw);
    }
    
    for (j=0; j<NCR; j++) write_xyz(xcr + j,ycr + j,zcr + j,fpw);
    
    if (fps == NULL)
      fclose(fpw);
  }
}
/****************************************************************************/
void write_conf(char *fn,char *fmode) 
{
  int i,j,nch=NCH,ncr=NCR,naa;
  FILE *fp;

  fp=fopen(fn,fmode);

  fwrite(&ind,sizeof(int),1,fp);
  fwrite(&nch,sizeof(int),1,fp);
  fwrite(&ncr,sizeof(int),1,fp);

  for (i=0; i<NCH; i++) {
    naa = Ne[i] - Nb[i] + 1; 
    fwrite(&naa,sizeof(int),1,fp);
    j = iN[Nb[i]];     write_xyz(x + j,y + j,z + j,fp);
    j = iN[Nb[i]] + 1; write_xyz(x + j,y + j,z + j,fp);
    j = iCa[Nb[i]];    write_xyz(x + j,y + j,z + j,fp);
    for (j=Nb[i]; j<=Ne[i]; j++) write_ph(ph + 3*j,ph + 3*j + 1,fp);
  }

  for (j=0; j<NCR; j++) write_xyz(xcr + j,ycr + j,zcr + j,fp);

  fclose(fp);
}
/****************************************************************************/
int read_conf_ver2(int iflag,char *fn,int m)
/* iflag < 0   open file  */
/* iflag > 0   close file */
/* iflag = 0   read next conformation (first m chains) */
{
  int i,j,nch,ncr,naa;
  static FILE *fp;

  if (iflag < 0) {
    fp = fopen(fn,"r");
    return 0;
  }

  if (iflag > 0) {
    if (fp != NULL)
      fclose(fp);
    return 0;
  }

  if (iflag == 0) {
    if (1 != fread(&ind,sizeof(int),1,fp) ) return 0;
    fread(&nch,sizeof(int),1,fp);
    fread(&ncr,sizeof(int),1,fp);
    fread(&E,sizeof(double),1,fp);
    
    //  printf("nch %i ncr %i %lf\n",nch,ncr,E);
    
    for (i=0; i<m; i++) {
      if (1 != fread(&naa,sizeof(int),1,fp)) return 0;
      j = iN[Nb[i]];     read_xyz(x + j,y + j,z + j,fp);
      j = iN[Nb[i]] + 1; read_xyz(x + j,y + j,z + j,fp);
      j = iCa[Nb[i]];    read_xyz(x + j,y + j,z + j,fp);
      for (j=Nb[i]; j<=Ne[i]; j++) read_ph(ph + 3*j,ph + 3*j + 1,fp);
    }
    
    for (j=0; j<NCR; j++) read_xyz(xcr + j,ycr + j,zcr + j,fp);
  }
  
  return 1;
}
/****************************************************************************/
int read_conf(int iflag,char *fn,int m)
/* iflag < 0   open file  */
/* iflag > 0   close file */
/* iflag = 0   read next conformation (first m chains) */
{
  int i,j,nch,ncr,naa;
  static FILE *fp;

  if (iflag<0) {
    fp=fopen(fn,"r");
    return 0;
  }

  if (iflag>0) {
    fclose(fp);
    return 0;
  }

  if (1 != fread(&ind,sizeof(int),1,fp) ) return 0;
  fread(&nch,sizeof(int),1,fp);
  fread(&ncr,sizeof(int),1,fp);

  for (i=0; i<m; i++) {
    if (1 != fread(&naa,sizeof(int),1,fp)) return 0;
    j = iN[Nb[i]];     read_xyz(x + j,y + j,z + j,fp);
    j = iN[Nb[i]] + 1; read_xyz(x + j,y + j,z + j,fp);
    j = iCa[Nb[i]];    read_xyz(x + j,y + j,z + j,fp);
    for (j=Nb[i]; j<=Ne[i]; j++) read_ph(ph + 3*j,ph + 3*j + 1,fp);
  }

  for (j=0; j<NCR; j++) read_xyz(xcr + j,ycr + j,zcr + j,fp);

  return 1;
}
/****************************************************************************/
void plot_rama(int ia,int ch)
{
  FILE *fp;
  while (ph[3*ia]>pi) ph[3*ia]-=pi2;
  while (ph[3*ia]<-pi) ph[3*ia]+=pi2;
  while (ph[3*ia+1]>pi) ph[3*ia+1]-=pi2;
  while (ph[3*ia+1]<-pi) ph[3*ia+1]+=pi2;
  fp=fopen(RAMA,"a"); 
  fprintf(fp,"%f %f\n",ph[3*ia]*180/pi,ph[3*ia+1]*180/pi);
  fclose(fp); 
}
/****************************************************************************/
void runtime(long it,double o[],int n)
{
  int i;
  FILE *fp;
  fp=fopen(RT,"a");
  fprintf(fp,"%li %i ",it,ind);
  for (i=1;i<n;i++) fprintf(fp,"%f ",o[i]);
  fprintf(fp,"\n");
  fclose(fp);
}
/****************************************************************************/
void runtime_allcols(long it,double o[],int n)
{
  int i;
  FILE *fp;
  fp=fopen(RT,"a");
  fprintf(fp,"%li %i ",it,ind);
  for (i=1;i<n;i++) {if (i!=7) fprintf(fp,"%f ",o[i]);}
  fprintf(fp,"\n");
  fclose(fp);
}
/****************************************************************************/
void relax_conf() { 
  int i,nstep = 0;
  long int tmpl[NTMP];
  double tmpf[NTMP];
  
  printf("Relaxing...\n");

  cellcalc("ecalc","",0,0,0,&Eev,&Ehb,&Ecc,&Ecp,0);
  printf("   start: Ecc + Ecp = %lf \n", Ecc + Ecp );
  printf("          Ehb       = %lf \n", Ehb );
  printf("          Eev       = %lf \n", Eev );
  fflush(stdout);
  
  relax = 1;
  
  while ( Ecc + Ecp > eclash || Eev > 10 * N || Ehb > 10 * N) {

    for (i = 0; i < MCCON; i++) {
      if ((NCR > 0 && genrand() < 0.5) || NCH == 0) { 
	crowd_translate();
     } else {
	pivot(0,tmpl+ind,tmpf+ind,tmpl+ind,tmpf+ind);
      }
    }

    /* normalize */

    int2cart_all(1);
    backup_coord();
    ecalc();
    
    nstep++;

    if (nstep % 10000 == 0) 
      printf("nstep %i: %lf %lf %lf\n",nstep, Ecc + Ecp, Ehb, Eev);

  }

  printf("     end: Ecc + Ecp = %lf \n", Ecc + Ecp );
  printf("          Ehb       = %lf \n", Ehb );
  printf("          Eev       = %lf \n", Eev );

  printf("...done after %i steps\n",nstep);

  relax = 0;

  return;
}
/****************************************************************************/
/***** INITIALIZATION *******************************************************/
/****************************************************************************/
void init(void)
{
  extern void init_aa(void);
  int i,j,k,c,b;
  int checkpnt_start=0;
  long *ldum=0;
  double *fdum=0;
  FILE *fp;

  /* various constants */   

  pi=acos(-1.);
  pi2=2*pi;
  pid2=pi/2;

  /* read amino acid sequences */

  if (NCH > 0) { 
    fp=fopen(INPUT,"r");
    for (i=j=0;i<NCH;i++) {
      Nb[i]=j;
      while ((c=getc(fp))!='\n') seq[j++]=c; 
      Ne[i]=j-1;
    }
    if (Ne[NCH-1]!=N-1) {printf("N wrong %i! %i\n",Ne[NCH-1],N); exit(0);}
    fclose(fp);
  }

  printf("Chains:\n");
  if (NCH > 0) {
    printf("Chain#   First      Last     Sequence \n");
    for (i=0;i<NCH;i++) {
      printf("%3i %9i %9i     ",i,Nb[i],Ne[i]); 
      for (j=0;j<=Ne[i]-Nb[i];j++) {
	printf("%c",seq[j+Nb[i]]); 
	if ((j+1)%5==0) putchar(' ');
      }
      printf("\n");
    }
  }
  printf("\n");

  printf("Crowders:\n");
  if (NCR > 0) {
    Phic = (double) NCR * 4./3. * pi * rcrowd * rcrowd * rcrowd / vbox;
    printf("  Number of crowders NCR %d\n",NCR);
    printf("  Radius %lf\n",rcrowd);
    printf("  Softness sigcr %lf\n",sigcr);
    printf("  Packing fraction Phi %lf\n",Phic);
    if (ISOFT) printf("  Attraction strength coefficient %f\n",epsilonatt);
    printf("  Repulsion strength coefficient %f\n",epsilonrep);
    printf("  Soft attraction %d\n",ISOFT);
  }
  printf("\n");

  /* temperatures */
  
  for (i=0;i<NTMP;i++)    
    beta[i]=1./TMAX*pow(TMAX/TMIN,(double)i/max(NTMP-1,1));

  /* g parameters */

  if (IREADG == 1 && NTMP > 1) {
    fp = fopen(INPUTG,"r");
    for (i=0; i<NTMP; i++) fscanf(fp,"%i %lf",&i,&g[i]);
    fclose(fp);
  }

  /* read native contact list */

  if (NCH > 0) {
    fp=fopen(CONTACTS, "r");
    ncontacts = k = 0;
    if (fp != NULL) {
      while (fscanf(fp, "%d %d %d",&i,&j,&b) == 3) { 
	contact_pair_i[k] = i;	
	contact_pair_j[k] = j;
	boxnum[k] = b;
	k++; 
      }
      ncontacts = k;
      if (ncontacts > MAXCONT) {
	printf("ncontacts too large, %i\n\n\n",ncontacts);
	exit(-1);
      }
      printf("Read file %s, ncontacts %i\n", CONTACTS, ncontacts);
      fclose(fp); 
    }
  }
  
  /* calculate atom indices */

  for (i=0;i<N;i++) {  
    iN[i]=(i==0) ? 0:iC[i-1]+2;
    iCa[i]=iN[i]+2;
    iC[i]=iCa[i]+3;    
    atyp[iN[i]]=1;    
    atyp[iN[i]+1]=4;
    atyp[iCa[i]+1]=3;
    atyp[iCa[i]]=atyp[iC[i]]=0;
    atyp[iC[i]+1]=2;
    atyp[iCa[i]+2]=5;
    if (seq[i]=='G') atyp[iCa[i]+2]=3;
    //    if (seq[i]=='P') atyp[iN[i]+1]=0;
  }
  if (NCH > 0) {
    if (iC[N-1]+2!=NTO) {printf("in init: %i -- NTO\n",iC[N-1]+2);exit(-1);}
  }
  
  for (j=0;j<NCH;j++) {
    iBeg[j]=iN[Nb[j]];
    iEnd[j]=iC[Ne[j]]+1;
  }

  /* i2a[],a2c[],idonacc[] */

  for (i = 0; i < NTO; i++)
    idonacc[i] = 0;

  for (i = 0; i < NCH; i++) { 
    for (j = Nb[i]; j <= Ne[i]; j++) {
      for (k = iN[j]; k <= iC[j] + 1; k++) i2a[k] = j;
      a2c[j] = i;
      idonacc[iN[j] + 1] = +1;
      idonacc[iC[j] + 1] = -1;
    }
  }

  /* gam[][] */

  for (i=0;i<N;i++) {
    for (j=0;j<N;j++) {
      if (seq[i]=='G' || seq[j]=='G')
	gam[i][j]=0.75;
      else
	gam[i][j]=1.00;
    }
  }

  /* set seed */ 

  if (ISEED==1){
    FILE *devrandom = fopen("/dev/urandom", "r");
    fread(&seed, sizeof(seed), 1, devrandom);
    seed=-labs(seed%100000000);
    fclose(devrandom);

    devrandom = fopen("/dev/urandom", "r");
    fread(&seed_mt, sizeof(seed_mt), 1, devrandom);
    fclose(fp);
  } 

  if (RAN_RAN3) {
    orig_seed=seed;
    printf("random number generator: ran3() -- See Numerical Recipes\n");
    printf("orig_seed = %li\n",seed);
  } else if (RAN_MARSENNE_TWISTER) {
    printf("random number generator: genran64_real3() -- Marsenne Twister\n");
    printf("orig_seed_mt = %lu\n",seed_mt);
  }

  /* initialize various subroutines */
  
  init_aa();
  int2cart(0,0);

  localsa("init","",0,0,&Eloc);
  bias("init","",0,0,&Ebias);
  hp("init","",0,0,&Ehp);
  cellcalc("init","",0,0,0,&Eev,&Ehb,&Ecc,&Ecp,0);
  
  sel_aa("init");
  sel_aa2(-1,0,0);
  pivot(-1,ldum,fdum,ldum,fdum);
  bgs(-1,ldum,fdum,ldum,fdum);
  trans(-1);
  
  cont_made(-1,0,0);
  cmdist(-1,0,0);
  //  hist_sq(-1,0,0);
  //  histo_cont(-1,0,0);
  //  histo_contmap(-1,0);
  histoe(-1,0,0);
  //  helixprof(-1,0);
  //  histod1(-1,0,0);
  //  histod2(-1,0,0);
  //  histod3(-1,0,0);
  //  histod4(-1,0,0);
  //  histoed1(-1,0,0,0);
  //  histoed2(-1,0,0,0);
  //  histoed3(-1,0,0,0);
  
  //  crowd_stat(-1);
  
  for (i=0;i<N;i++) ph[3*i+2]=pi; 

  if (RMSD > 0) {
    printf("init rmsd...\n");
    readpdb(1,NATIVE1,xnat,ynat,znat,NR1,NR2);
    if (RMSD > 1) readpdb(1,NATIVE2,xnat2,ynat2,znat2,NR1,NR2);
    printf("    rmsd DONE.\n");
  }

  for (j=0;j<NOBS;j++) minobs[j]=1e4;

  /* initialize dynamical variables */

  if ( (fp = fopen(CHECKCONF,"r")) != NULL ) {
    printf("Reading checkpoint\n");
    fclose(fp);
    read_checkpnt(&imc0,so,minobs);
    checkpnt_start=1;
  } else {

    if (ISTART == 0) { 
      printf("*** Disordered start ***\n\n");
      
      for (i=0;i<N;i++) { 
	ph[3*i]=pi2*genrand(); 
	ph[3*i+1]=pi2*genrand(); 
      }
      
      int2cart_all(1);
      for (j=0;j<NCH;j++) {
	movech(BOX*genrand(),BOX*genrand(),BOX*genrand(),j);
      }
      
      for (i=0;i<NCR;i++) {
	xcr[i] = BOX * genrand();
	ycr[i] = BOX * genrand();
	zcr[i] = BOX * genrand();
      }
      
      ind = NTMP*genrand(); 
      
    } else if (ISTART == 1) {
      printf("*** Helix start ***\n\n");
      
      for (i=0;i<N;i++) {
	ph[3*i]=-58.0*pi/180.; 
	ph[3*i+1]=-47.0*pi/180.; 
      }
      int2cart_all(1);
      for (j=0;j<NCH;j++) {
	movech(BOX*genrand(),BOX*genrand(),BOX*genrand(),j);
      }
      
      for (i=0;i<NCR;i++) { /* crowders */
	xcr[i] = BOX * genrand();
	ycr[i] = BOX * genrand();
	zcr[i] = BOX * genrand();
      }
      
      ind = NTMP*genrand(); 
    } else if (ISTART == 2) {
      printf("*** Start from conf file: %s ***\n\n",STARTCONF);
      
      read_conf_ver2(-1,STARTCONF,NCH);
      if (0 == read_conf_ver2(0,STARTCONF,NCH)) {printf("Exiting...\n"); exit(-1);}
      read_conf_ver2(1,STARTCONF,NCH);
      int2cart_all(1);
    }  else {
      printf("ISTART wrong %i\n",ISTART);
      exit(-1);
    }    
  }

  /* backup coordinates */

  backup_coord();
  ecalc();

  //   exit(-1);
  
  if (!checkpnt_start) {
    relax_conf();
  }
  
  int2cart_all(1);
  backup_coord();
  ecalc();

  printf("Initial temperature index: ind = %i\n\n",ind);
  printf("Initial conformation:\n");
  printf("  E   %8.15lf\n",E);
  if (NCH > 0) {
    printf("  Eev %8.15lf Ebias %8.15lf Ehb  %8.15lf\n",Eev,Ebias,Ehb);
    printf("  Ehp %8.15lf Eloc  %8.15lf \n",Ehp,Eloc);
  }
  if (NCR > 0) {
    printf("  Ecc %8.15lf Ecp %8.15lf\n",Ecc,Ecp);
  }
  //  if (NCH > 0) printf(" Rg in the native state: %8.15lf\n ",RGYR);

  dumppdb("start.pdb",0,0,fdum,"ALL");
}
/****************************************************************************/
void write_averages(char *fn,int ntmp,int nobs,double so[][NOBS]) {
  int i,j;
  FILE *fp;

  fp = fopen(fn,"w");
  for (i = 0; i < ntmp; i++) {
    fprintf(fp,"%f %i ",1./beta[i],i);
    for (j = 1; j < nobs && so[i][0] > 0; j++) 
      fprintf(fp,"%lf ",so[i][j] / so[i][0]);
    fprintf(fp,"\n");
  }
  fclose(fp);
  
}
/****************************************************************************/
int acceptance(long *nflip,double *accflip,
	       long *ncrd,double *acccrd,
	       long *npivup,double *accpivup,
	       long *npivdw,double *accpivdw,
	       long *nbgsup,double *accbgsup,
	       long *nbgsdw,double *accbgsdw,
	       long *ntrl,double *acctrl) {
  int i;
  FILE *fp;

  fp = fopen(ACCSTAT,"w");
  
  /* write acceptance probabilities */

  fprintf(fp,"Acceptance\n"); 
  fprintf(fp,"==========\n\n"); 
  for (i = 0; i < NTMP; i++) {
    fprintf(fp,"temp %i  ",i);
    fprintf(fp,"nflip %10ld  acc flip ",nflip[i]);
    if (nflip[i] > 0) fprintf(fp,"%5.3f       ", accflip[i] / nflip[i]);
    else fprintf(fp,"            ");
    fprintf(fp,"\n");
  }
  for (i = 0; i < NTMP; i++) {
    fprintf(fp,"temp %i  ",i);
    fprintf(fp,"ncrd %10ld  acc crd ",ncrd[i]);
    if (ncrd[i] > 0) fprintf(fp,"%5.3f       ", acccrd[i] / ncrd[i]);
    else fprintf(fp,"            ");
    fprintf(fp,"\n");
  }
  fprintf(fp,"\n");
  for (i = 0; i < NTMP; i++) {
    fprintf(fp,"temp %i  ",i);
    fprintf(fp,"npivup %10ld  acc pivup ",npivup[i]);
    if (npivup[i] > 0) fprintf(fp,"%5.3f       ", accpivup[i] / npivup[i]);
    else fprintf(fp,"            ");
    fprintf(fp,"npivdw %10ld  acc pivdw ",npivdw[i]);
    if (npivdw[i] > 0) fprintf(fp,"%5.3f       ", accpivdw[i] / npivdw[i]);
    else fprintf(fp,"            ");
    fprintf(fp,"\n");
  }
  for (i = 0; i < NTMP; i++) {
    fprintf(fp,"temp %i  ",i);
    fprintf(fp,"nbgsup %10ld  acc bgsup ",nbgsup[i]);
    if (nbgsup[i]>0) fprintf(fp,"%5.3f       ",accbgsup[i]/nbgsup[i]);
    else fprintf(fp,"            ");
    fprintf(fp,"nbgsdw %10ld  acc bgsdw ",nbgsdw[i]);
    if (nbgsdw[i]>0) fprintf(fp,"%5.3f       ",accbgsdw[i]/nbgsdw[i]);
    else fprintf(fp,"            ");
    fprintf(fp,"\n");
  }
  for (i = 0;i < NTMP; i++) {
    fprintf(fp,"temp %i  ",i);
    fprintf(fp,"ntrl   %10ld  acc trl   ",ntrl[i]);
    if (ntrl[i]>0) fprintf(fp,"%5.3f       ",acctrl[i]/ntrl[i]);
    else fprintf(fp,"            ");
    fprintf(fp,"\n");
  }


  fclose(fp);
    
  return 0;
}
/****************************************************************************/
void init_aa(void)
{
  int i,j,n;

  for (i=0;i<NTO;i++) {for (j=i+1;j<NTO;j++) aa[i][j]=1;} 
  for (i=0;i<NTO;i++) aa[i][i]=-1;

  for (j=0;j<NCH;j++) {
    for (i=Nb[j];i<=Ne[j];i++) {
      if (seq[i]!='P') aa[iN[i]][iN[i]+1]=-1;     /* NH */

      aa[iN[i]][iCa[i]]=-1;                       /* NCa */
      aa[iN[i]][iCa[i]+1]=-1;                     /* NHa */
      aa[iN[i]][iCa[i]+2]=-1;                     /* NCb/Ha2 */
      aa[iN[i]][iC[i]]=-1;                        /* NC */
      if (seq[i]!='P') aa[iN[i]+1][iCa[i]]=-1;    /* HCa */
      aa[iCa[i]][iCa[i]+1]=-1;                    /* CaHa */
      aa[iCa[i]][iCa[i]+2]=-1;                    /* CaCb/Ha2 */
      aa[iCa[i]][iC[i]]=-1;                       /* CaC */
      aa[iCa[i]][iC[i]+1]=-1;                     /* CaO */
      if (i<Ne[j]) aa[iCa[i]][iN[i+1]]=-1;        /* CaN */
      if (i>Nb[j] && seq[i]!='P') 
	aa[iCa[i-1]][iN[i]+1]=-1;                 /* CaH */
      if (i<Ne[j]) aa[iCa[i]][iCa[i+1]]=-1;       /* CaCa */
      aa[iCa[i]+1][iCa[i]+2]=-1;                  /* HaCb/Ha2 */
      aa[iCa[i]+1][iC[i]]=-1;                     /* HaC */
      aa[iCa[i]+2][iC[i]]=-1;                     /* Cb/Ha2C */
      aa[iC[i]][iC[i]+1]=-1;                      /* CO */
      if (i<Ne[j]) aa[iC[i]][iN[i+1]]=-1;         /* CN */
      if (i>Nb[j] && seq[i]!='P') 
	aa[iC[i-1]][iN[i]+1]=-1;                  /* CH */
    
      if (i<Ne[j]) aa[iC[i]][iCa[i+1]]=-1;        /* CCa*/
      if (i<Ne[j]) aa[iC[i]+1][iN[i+1]]=-1;       /* ON */
      if (i>Nb[j] && seq[i]!='P')
	aa[iC[i-1]+1][iN[i]+1]=-1;                /* OH */
      if (i<Ne[j]) aa[iC[i]+1][iCa[i+1]]=-1;      /* OCa */
    }
  }


  /* distances fixed because of fixed Pro phi angle */

  for (j=0;j<NCH;j++) {  
    for (i=Nb[j]+1;i<=Ne[j];i++) {
      // check
      if (seq[i]=='P') {
	aa[iCa[i-1]][iCa[i]+1]=-1;   /* CaHa */
	aa[iCa[i-1]][iCa[i]+2]=-1;   /* CaCb/Ha2 */
	aa[iCa[i-1]][iC[i]]=-1;      /* CaC */
	aa[iC[i-1]][iCa[i]+1]=-1;    /* CHa */
	aa[iC[i-1]][iCa[i]+2]=-1;    /* CCb/Ha2 */
	aa[iC[i-1]][iC[i]]=-1;       /* CC */
	aa[iC[i-1]+1][iCa[i]+1]=-1;  /* OHa */
	aa[iC[i-1]+1][iCa[i]+2]=-1;  /* OCb/Ha2 */
	aa[iC[i-1]+1][iC[i]]=-1;     /* OO */
      }
    }
  }

  /* phi pairs separated by 3 bonds */

  for (j=0;j<NCH;j++) {  
    for (i=Nb[j];i<Ne[j];i++) {
      if (seq[i+1]!='P') {
	aa[iC[i]][iCa[i+1]+2]=2;          /* CpCb/Ha2 */
	aa[iC[i]][iC[i+1]]=2;             /* CpCp */
	aa[iC[i]][iCa[i+1]+1]=2;          /* CpHa */
      }
    }
  }

  for (j=0;j<NCH;j++) {  
    for (i=Nb[j];i<=Ne[j];i++) {
      if (seq[i]!='P') {
	aa[iN[i]+1][iCa[i]+2]=2;          /* HCb/Ha2 */
	aa[iN[i]+1][iC[i]]=2;             /* HCp */
	aa[iN[i]+1][iCa[i]+1]=2;          /* HHa */
      }
    }
  }

  /* psi pairs separated by 3 bonds */

  for (j=0;j<NCH;j++) {  
    for (i=Nb[j];i<=Ne[j];i++) {
      aa[iCa[i]+2][iC[i]+1]=2;            /* Cb/Ha2O */
      aa[iN[i]][iC[i]+1]=2;               /* NO */
      aa[iCa[i]+1][iC[i]+1]=2;            /* HaO */
    }
  }

  for (j=0;j<NCH;j++) {  
    for (i=Nb[j];i<Ne[j];i++) {
      aa[iCa[i]+2][iN[i+1]]=2;            /* Cb/Ha2N */
      aa[iN[i]][iN[i+1]]=2;               /* NN */
      aa[iCa[i]+1][iN[i+1]]=2;            /* HaN */
    }
  }

  /* HH and OO repulsions */ 

  //  for (i=0;i<N-1;i++) {
  //    if (seq[i]!='P' && seq[i+1]!='P') {
  //    aa[iN[i]+1][iN[i+1]+1]=2;         /* HH */
  //  }
  //  aa[iC[i]+1][iC[i+1]+1]=2;           /* OO */
  // }

  /* extra helix breaking term for proline */ 

  /*  for (i=1;i<N;i++) {
      if (seq[i]=='P') {
      aa[iN[i-1]][iCa[i]+4]=3;                // NCd 
      }
      }*/
  
  n=0;
  for (i=0;i<NTO;i++) {
    for (j=i+1;j<NTO;j++) {
      if (aa[i][j]==1) {
        aa[i][j]=atyp[i]+atyp[j]*NA;
        aa[j][i]=atyp[j]+atyp[i]*NA;
        n++;
      }
      else if (aa[i][j]==2) {
        aa[i][j]=atyp[i]+atyp[j]*NA+NA*NA;
        aa[j][i]=atyp[j]+atyp[i]*NA+NA*NA;
        n++;
      }
      else if (aa[i][j]==-1) {
        aa[j][i]=-1;
      }
      else {
        printf("aa wrong\n");exit(-1);
      }
    }
  }

  printf("no of atom pairs %i   no of distances %i\n",NTO*(NTO-1)/2,n);
}
/****************************************************************************/
