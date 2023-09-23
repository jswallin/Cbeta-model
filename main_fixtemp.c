# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include "defs.h"
# include "sys.h"
# include "global.h"
# include "main.h"
/****************************************************************************/
int main (int argc,char *argv[])
{
  int i,j;
  long nflip[NTMP],npivup[NTMP],npivdw[NTMP],nbgsup[NTMP],nbgsdw[NTMP],ncrd[NTMP];
  long ntrl[NTMP],imc,nnorm=0;
  double acctrl[NTMP],accflip[NTMP];
  double accpivup[NTMP],accpivdw[NTMP],accbgsup[NTMP],accbgsdw[NTMP],acccrd[NTMP];
  double o[NOBS],so[NTMP][NOBS];
  double tmp,eold,roff=0,tmpx=0,tmpy=0,tmpz=0;
  double rg1,rmsd1=0,ncont=0;
  double minobs[NOBS];
  //  int maxc;
  /* Set fractions of MC move types: */
  /* frcrd = fraction of *all moves* that are crowder moves (remaining are protein moves). */
  /* frbgs = fraction of *protein moves* that are BGS moves. */
  /* frtrl = fraction of *protein moves* that are translation moves. */
  /* frpiv = fraction of *protein moves* that are pivot moves */
  double frcrd[NTMP]={0.5};
  double frbgs[NTMP]={0.0};
  double frtrl[NTMP]={0.0};
  double frpiv[NTMP]={1.0};

  printf("\nCommand line: "); 
  for (i=0;i<argc;i++) printf("%s ",argv[i]);
  printf("\n\n");
  fflush(stdout);

  if (NTMP > 1) {printf("main_fixtemp: NTMP>1\n"); exit(-1);}

  for (i=0;i<NTMP;i++) {
    nflip[i]=npivup[i]=npivdw[i]=nbgsup[i]=nbgsdw[i]=ntrl[i]=ncrd[i]=0;
    accflip[i]=accpivup[i]=accpivdw[i]=accbgsup[i]=accbgsdw[i]=acccrd[i]=0; 
    acctrl[i]=0;
    if (NCH == 0) frcrd[i] = 1; /* only crowders */
    if (NCR == 0) frcrd[i] = 0; /* no crowders */
    frpiv[i] = 1 - frtrl[i] - frbgs[i]; /* normalize */
  }
  
  init();
  printheading(frcrd,frpiv,frbgs,frtrl);

  for (imc = imc0; imc<MCCYC; imc++) {
    
    for (i=0; i<MCCON;i++) {
      if ( genrand() < frcrd[ind] ) { 
	ncrd[ind]++;
	acccrd[ind] += crowd_translate(0);
      } else { 
	if ((tmp = genrand())<frtrl[ind]) 
	  {ntrl[ind]++; acctrl[ind]+=trans(0);} 
	else if (tmp<frtrl[ind]+frbgs[ind]) 
	  {bgs(0,nbgsup+ind,accbgsup+ind,nbgsdw+ind,accbgsdw+ind);} 
	else if (tmp<frtrl[ind]+frbgs[ind]+frpiv[ind]) 
	  {pivot(0,npivup+ind,accpivup+ind,npivdw+ind,accpivdw+ind);}
      }

      if ( (i + 1) % INORM == 0 && NCH > 0 ) {
	nnorm++;
	int2cart_all(1);
	tmpx=tmpy=tmpz=0;
	for (j=0;j<NCH;j++) {
	  tmpx+=fabs(xo[iEnd[j]]-x[iEnd[j]]);
	  tmpy+=fabs(yo[iEnd[j]]-y[iEnd[j]]);
	  tmpz+=fabs(zo[iEnd[j]]-z[iEnd[j]]);
	}
	roff+=(tmp=(tmpx+tmpy+tmpz)/NCH);
	backup_coord();
	eold=E;
	Eevo=Eev; Eloco=Eloc; Ebiaso=Ebias; Ehbo=Ehb; Ehpo=Ehp; Ecco=Ecc; Ecpo=Ecp;
	ecalc();
	if (tmp>1e-6 || fabs(eold-E)>1e-6) {
          printf("normalization %li\n",imc); 
          printf("tmp: %e %e %e\n",tmpx,tmpy,tmpz);
          printf("e:   %lf   %lf \n",eold,E);
          printf("eev:   %lf   %lf \n",Eevo,Eev);
          printf("eloc:   %lf   %lf \n",Eloco,Eloc);
          printf("ebias:   %lf   %lf \n",Ebiaso,Ebias);
          printf("ehb:   %lf   %lf \n",Ehbo,Ehb);
          printf("ehp:   %lf   %lf \n",Ehpo,Ehp);
	  //          printf("ecc:   %lf   %lf \n",Ecco,Ecc);
	  //          printf("ecp:   %lf   %lf \n",Ecpo,Ecp);
	  fflush(stdout);
	}
      }
    }

    ch2box();
    cr2box();
    
    rmsd1 = calcrmsd_ca(0,xnat,ynat,znat,Nb[0],Ne[0]); 
    rg1 = rgyr(0);    
    ncont = contacts();

    o[1] = E; o[2] = Eloc; o[3] = Ebias; o[4] = Eev; o[5] = Ehb; o[6] = Ehp;
    o[7] = Ecc; o[8] = Ecp;
    o[9] = helix(0); o[10] = sheet(0);
    o[10] = rmsd1;
    o[11] = ncont;
    o[12] = rg1;

    if (imc+1>NTHERM) {
      so[ind][0]++; for (i=1;i<NOBS;i++) so[ind][i]+=o[i];

      /* histograms */

      histoe(0,ind,E);
    }

    if ((imc+1)%ICHECK==0) {    
      ch2box();
      cr2box();
      
      /* write averages */

      write_averages(AVERAGES,NTMP,NOBS,so);

      /* write histograms */
      
      histoe(1,0,0);

      /* MC statistics */

      pivot(1,npivup,accpivup,npivdw,accpivdw);
      bgs(1,nbgsup,accbgsup,nbgsdw,accbgsdw);
      trans(1);
      
      acceptance(nflip,accflip,ncrd,acccrd,npivup,accpivup,npivdw,accpivdw,
		 nbgsup,accbgsup,nbgsdw,accbgsdw,ntrl,acctrl);

    }
    
    if ((imc+1)%IRT==0) {
      dumppdb(OUTPDB,imc,NOBS,o,"ALL"); 
      write_conf_ver2(0,OUTCONF,"w"); 
      runtime(imc+1,o,NOBS);
    } 

    if ((imc+1)%ICONF==0) write_conf_ver2(0,CONF,"a"); 

    /* minimum conformations */

    if (E<minobs[0]) {
      minobs[0]=E; 
      dumppdb("mine.pdb",imc,NOBS,o,"ALL");
      write_conf_ver2(0,"mine.conf","w");
    } 

    /* write checkpoint data */

    if ((imc+1)%ICHECK==0) write_checkpnt(&imc,so,minobs);
  }

  /* end run */

  printf("run over\n");
  printf("average error in last atom %e\n\n",roff/nnorm);

  /* write averages */

  write_averages(AVERAGES,NTMP,NOBS,so);
  
  /* write histograms */

  printf("Writing histograms \n");
  
  histoe(2,0,0);

  /* MC statistics */
  
  printf("Writing acceptance rates \n");

  trans(1);
  pivot(1,npivup,accpivup,npivdw,accpivdw);
  bgs(1,nbgsup,accbgsup,nbgsdw,accbgsdw);
  
  acceptance(nflip,accflip,ncrd,acccrd,npivup,accpivup,npivdw,accpivdw,
	     nbgsup,accbgsup,nbgsdw,accbgsdw,ntrl,acctrl);
  
  return 0;
}
/****************************************************************************/
