# include <time.h>
# include <stdio.h>
# include <ctype.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include "mt64.h"
# include "defs.h"
# include "sys.h"
# include "global.h"
# include "main.h"
/****************************************************************************/
int main (int argc,char *argv[])
{
  int i,j;
  double acctrl[NTMP],accflip[NTMP];
  double o[NOBS];
  double tmp,eold,roff=0,tmpx,tmpy,tmpz;
  double rg1,rmsd1=0,rmsd2=0,ncont=0;
  double g2[NTMP];
  long nflip[NTMP],npivup[NTMP],npivdw[NTMP],nbgsup[NTMP],nbgsdw[NTMP],ncrd[NTMP];
  long ntrl[NTMP],imc,nnorm=0;
  double accpivup[NTMP],accpivdw[NTMP],accbgsup[NTMP],accbgsdw[NTMP],acccrd[NTMP];
  //int maxc;

  /* Set fractions of MC move types: */
  /* frcrd = fraction of *all moves* that are crowder moves (remaining are protein moves). */
  /* frbgs = fraction of *protein moves* that are BGS moves. */
  /* frtrl = fraction of *protein moves* that are translation moves. */
  /* frpiv = fraction of *protein moves* that are pivot moves */
  double frcrd[NTMP]={0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50,0.50};
  double frbgs[NTMP]={0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00};
  double frpiv[NTMP]={1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00,1.00};
  double frtrl[NTMP]={0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00,0.00};  

  int unf1,unf2,nat1,nat2,tse1,tse2;
  FILE *fp;
  
  printf("\nCommand line: "); 
  for (i=0;i<argc;i++) printf("%s ",argv[i]);
  printf("\n\n");
  fflush(stdout);

  fp = fopen("unf_tse_nat.def","r");
  if (fp != NULL) {
    fscanf(fp,"%d %d %d %d %d %d",&unf1,&unf2,&tse1,&tse2,&nat1,&nat2);
    printf("unf: %d-%d; tse: %d-%d; nat: %d-%d\n\n",unf1,unf2,tse1,tse2,nat1,nat2);
    fclose(fp);
  }
  
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

  for (imc = imc0; imc < MCCYC; imc++) { 
    nflip[ind]++; accflip[ind] += flip(E);
    
    for (i = 0; i < MCCON; i++) {

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
    
    if (imc + 1 > NTHERM) {
      
      /* observables */      

      so[ind][0]++; for (i=1;i<NOBS;i++) so[ind][i] += o[i];
      if (ncont >= nat1) {
	son[ind][0]++; for (i=1;i<NOBS;i++) son[ind][i] += o[i];
      }
      if (ncont < nat1) {
	sou[ind][0]++; for (i=1;i<NOBS;i++) sou[ind][i] += o[i];
      }      
      
      /* histograms */
      
      histo_contmap(0,ind);
      histo_cont(0,ind,ncont);
      histoe(0,ind,E);
      histod1(0,ind,rg1);
      if (ncont >= nat1)  histod2(0,ind,rg1);
      if (ncont < nat1)  histod3(0,ind,rg1); 
    }

    if ((imc+1)%ICHECK==0) {    
      ch2box();
      cr2box();
      
      /* new g parameters */
      
      fp=fopen(OUTPUTG,"w");
      for (i=0;i<NTMP;i++) 
	g2[i] = g[i] + log(max(so[i][0]/(imc+1-NTHERM),0.01/NTMP));
      for (i=0;i<NTMP;i++) fprintf(fp,"%i %f\n",i,g2[i]-g2[NTMP-1]);
      fclose(fp);
	
      /* write averages */

      write_averages(AVERAGES,NTMP,NOBS,so);
      write_averages("averages_nat",NTMP,NOBS,son);
      write_averages("averages_unf",NTMP,NOBS,sou);
      write_averages("averages_trans",NTMP,NOBS,sot);

      /* write histograms */
      
      histo_contmap(1,0);
      histo_cont(1,0,0);
      histoe(1,0,0);
      histod1(1,0,0);
      histod2(1,0,0);
      histod3(1,0,0);
      histod4(1,0,0);

      crowd_stat(1);

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

    if ((imc+1)%IRAMA==0 && ind==NTMP-1) plot_rama(5,0);  
    if ((imc+1)%ICONF==0) write_conf_ver2(0,CONF,"a"); 

    /* minimum conformations */

    if (E<minobs[0]) {
      minobs[0]=E; 
      dumppdb("mine.pdb",imc,NOBS,o,"ALL");
      write_conf_ver2(0,"mine.conf","w");
    } 
    if (rmsd1<minobs[1]) {minobs[1]=rmsd1; dumppdb("rmin1.pdb",imc,NOBS,o,"ALL");}
    if (rmsd2<minobs[2]) {minobs[2]=rmsd2; dumppdb("rmin2.pdb",imc,NOBS,o,"ALL");}

    /* write checkpoint data */

    if ( (imc + 1) % ICHECK == 0) write_checkpnt(&imc,so,minobs);

  }

  /* end run */

  printf("run over\n");
  printf("average error in last atom %e\n\n",roff/nnorm);
  printf("\n\n========\n\n");

  printf("Writing averages\n");

  /* write averages */

  write_averages(AVERAGES,NTMP,NOBS,so);
  write_averages("averages_nat",NTMP,NOBS,son);
  write_averages("averages_unf",NTMP,NOBS,sou);
  write_averages("averages_trans",NTMP,NOBS,sot);
  
  /* write histograms */

  printf("Writing histograms \n");
  
  histo_contmap(2,0);
  histo_cont(2,0,0);
  histoe(2,0,0);
  histod1(2,0,0);
  histod2(2,0,0);
  histod3(2,0,0);
  histod4(2,0,0);

  /* new g parameters */

  printf("Writing new g parameters \n");
  
  fp=fopen(OUTPUTG,"w");
  for (i=0;i<NTMP;i++) g[i] += log(max(so[i][0]/(MCCYC-NTHERM),0.01/NTMP));
  for (i=0;i<NTMP;i++) fprintf(fp,"%i %f\n",i,g[i]-g[NTMP-1]);
  fclose(fp);

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
