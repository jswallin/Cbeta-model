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
/****************************************************************************/
/***** ENERGIES *************************************************************/
/****************************************************************************/
void ecalc() {
  cellcalc("ecalc","",0,0,0,&Eev,&Ehb,&Ecc,&Ecp,0);
  localsa("ecalc","",0,0,&Eloc);
  bias("ecalc","",0,0,&Ebias);
  hp("ecalc","",0,0,&Ehp);
  E = Eev + Ehb + Ecc + Ecp + Eloc + Ebias + Ehp;
}
/****************************************************************************/
/* Crowders */
/****************************************************************************/
double crowd_crowd_ecalc(int iflag,double r2) {
  double rcalc,rcalc3,rcalc6,e;
  static double Dcrowd,sref,inf2;
  
  if (NCR == 0) return 0;
  
  if (iflag < 0) {
    Dcrowd = 2 * rcrowd; 
    sref = 2 * sigcr;
    inf2 = (Dcrowd - sref) * (Dcrowd - sref);
    printf("<crowd_crowd> rcrowd %lf hard core %lf cutoff %lf\n",
	   rcrowd,sqrt(inf2),cut_crowd_crowd);
    crowd_crowd_ecalc(1,0);
    return 0;
  }

  if (iflag == 0) {

    if (r2 < inf2)
      return eclash;
    
    rcalc = sref / (sqrt(r2) - Dcrowd + sref);
    rcalc3 = rcalc * rcalc * rcalc ;
    rcalc6 = rcalc3 * rcalc3;
    e = rcalc6 * rcalc6;
    
    return min(e, eclash);
  }


  if (iflag > 0) {
    double r, del = cut_crowd_crowd / 200;
    FILE *fp;

    fp = fopen("_crowd_crowd_ecalc.pot","w");
    for (r = 0; r <= cut_crowd_crowd; r += del) 
      fprintf(fp,"%lf %lf\n",r,crowd_crowd_ecalc(0,r*r));
    fclose(fp);

    return 0;
  }

  return 0;
}
/****************************************************************************/
double crowd_atom_ecalc(int iflag,int atmtyp,double r2) {
  int i;
  double rcalc,rcalc3,rcalc6,e = 0;
  static double Dac[NA],srefat[NA],inf2; 

  if (NCR == 0) return 0;
  
  if (iflag < 0) {
    for (i = 0; i < NA; i++) {
      Dac[i] = rcrowd + sigsa[i];    /* range of interaction */
      srefat[i] = sigcr + sigsa[i];  /* softness of crowder + atom */
    }
    inf2 = (rcrowd - sigcr) * (rcrowd - sigcr);
    printf("<crowd_atom> epsilonrep %lf hard core %lf cutoff %lf  \n",epsilonrep,sqrt(inf2),cut_crowd_atom);
    for (i=0; i<NA; i++) printf("<crowd_atom>  Dac[%2d]    = %lf \n",i,Dac[i]);
    for (i=0; i<NA; i++) printf("<crowd_atom>  srefat[%2d] = %lf \n",i,srefat[i]);
    crowd_atom_ecalc(1,0,0);
    return 0;
  }

  if (iflag == 0) {

    if (r2 < inf2)
      return eclash;
    
    rcalc = srefat[atmtyp] / (sqrt(r2) - Dac[atmtyp] + srefat[atmtyp]);
    rcalc3 = rcalc * rcalc * rcalc ;
    rcalc6 = rcalc3 * rcalc3;
    
    e = rcalc6 * rcalc6;  
    
    return min(e, eclash);
  }

  if (iflag > 0) {
    double r, del = cut_crowd_atom / 200;
    FILE *fp;
    
    fp = fopen("_crowd_atom_ecalc.pot","w");
    for (i = 0; i < NA; i++) {
      for (r = 0; r <= cut_crowd_atom; r += del) 
	fprintf(fp,"%i %lf %lf\n",i,r,crowd_atom_ecalc(0,i,r*r));
    }
    fclose(fp);
    return 0;
  }

  return 0;
}
/****************************************************************************/
double crowd_atom_attr_ecalc(int iflag,int atmtyp,char aatyp,double r2) {
  int i;
  double e = 0,sigma;
  static double Dac[NA];
  
  if (NCR == 0 || ISOFT == 0) return 0;
  
  if (iflag < 0) {
    for (i = 0; i < NA; i++) 
      Dac[i] = rcrowd + sigsa[i];    /* range of interaction */
    printf("<crowd_atom_attr> epsilonatt %lf  \n",epsilonatt);
    crowd_atom_attr_ecalc(1,0,0,0);
    return 0;
  }

  if (iflag == 0) {

    //    if (atmtyp == 5 && aatyp == 'L') {  /* hydrophobic crowder */
    if (atmtyp == 5) {         /* nonspecific crowder */
      sigma = sqrt(r2) - Dac[atmtyp];
      e -= exp(- sigma * sigma / 2);
    }
    
    return epsilonatt * e;
  }

  if (iflag > 0) {
    double r, del = cut_crowd_atom / 200;
    FILE *fp;
    
    fp = fopen("_crowd_atom_attr_ecalc.pot","w");
    for (r = 0; r <= cut_crowd_atom; r += del) 
      fprintf(fp,"%lf %lf %lf %lf\n",r,
	      crowd_atom_attr_ecalc(0,5,'L',r*r),
	      crowd_atom_ecalc(0,5,r*r),
	      crowd_atom_ecalc(0,5,r*r) + crowd_atom_attr_ecalc(0,5,'L',r*r));
    fclose(fp);
    
    return 0;
  }

  return 0;
}
/****************************************************************************/
/* Chains */
/****************************************************************************/
void localsa(char *move,char *mc,int a1,int a2,double *esa)
{
  int a,i,j,k,kbeg,kend;
  double r2,r6,e = 0;
  static int i1[MAXLOC],i2[MAXLOC],id[MAXLOC],ip1[N],ip2[N],npair;
  static double cut2;
  static double sig2[2*NA*NA],asa[2*NA*NA],bsa[2*NA*NA];

  if (NCH == 0) return ;
    
  if (!strcmp(move,"init")) {
    npair = 0;
    for (k = 0; k < N; k++) {
      ip1[k] = npair;
      for (i = iN[k]; i <= iC[k] + 1; i++) {
	for (j = i + 1; j < NTO; j++) {
	  if ( (a = aa[i][j]) >= NA *NA ) {
	    i1[npair] = i;
	    i2[npair] = j;
	    id[npair++] = a;
	  }
	}
      }
      ip2[k] = npair - 1;
    }
    
    //    npair=0;
    //    for (i=0;i<NTO;i++) {
    //      for (j=i+1;j<NTO;j++) {
    //	if ((a=aa[i][j])>=NA*NA) {
    //          i1[npair]=i;
    //          i2[npair]=j;
    //          id[npair++]=a;
    //	}
    //      }
    //    } 
    
    cut2 = cut * cut;
    for (i = 0; i < NA; i++) {
      for (j = 0; j < NA; j++) {
	sig2[i + j * NA + NA * NA] = pow(sigsa[i] + sigsa[j], 2.0);
      }
    }
    for (i = NA * NA; i < 2 * NA * NA; i++) {
      asa[i] = -7 * pow(sig2[i] / cut2,6.0);
      bsa[i] = 6 * pow(sig2[i] / cut2,6.0) / cut2;
    } 
    printf("# local pairs %i cut %f\n",npair,sqrt(cut2));
    if (npair > MAXLOC) {printf("%i -- too many local pairs\n",npair);exit(-1);}

    return ;
  }

  if (!strcmp(move,"ecalc")) {

    for (k = 0; k < npair; k++) {
      i = i1[k]; j = i2[k];
      r2 = ( (x[i] - x[j]) * (x[i]-x[j]) +
	     (y[i] - y[j]) * (y[i]-y[j]) +
	     (z[i] - z[j]) * (z[i]-z[j]) );
      if (r2 > cut2) continue;
      r6 = sig2[(a = id[k])] / r2; r6 = r6 * r6 * r6;
      e += r6 * r6 + asa[a] + bsa[a]*r2;
    }

    (*esa) = ksa * e;

    return ;
  }

  if (!strcmp(move,"piv") || !strcmp(move,"bgs")) {

    kbeg = (a1 > Nb[a2c[a1]] ? ip1[a1 - 1] : ip1[a1]);
    kend = ip2[a2];
    
    for (k = kbeg; k <= kend; k++) {
      i = i1[k]; j = i2[k];
      r2 = ( (x[i] - x[j]) * (x[i]-x[j]) +
	     (y[i] - y[j]) * (y[i]-y[j]) +
	     (z[i] - z[j]) * (z[i]-z[j]) );
      if (r2 > cut2) continue;
      r6 = sig2[(a = id[k])] / r2; r6 = r6 * r6 * r6;
      e += r6 * r6 + asa[a] + bsa[a]*r2;
    }

    if (!strcmp(mc,"before")) 
      (*esa) -= ksa * e;

    if (!strcmp(mc,"after")) 
      (*esa) += ksa * e;

    return ;
  }

}
/****************************************************************************/
void bias(char *move,char *mc,int a1,int a2,double *ebi) {
  int i,j,k,l,m;
  double r2,e = 0;
  static double q[4];

  if (NCH == 0) return ;

  if (!strcmp(move,"init")) {
    q[0] = +0.42; // C
    q[1] = -0.42; // O  
    q[2] = -0.20; // N  
    q[3] = +0.20; // H
    return ;
  }

  if (!strcmp(move,"ecalc")) {
    for (i = 0; i < N; i++) {
      for (j = 0; j < 2; j++) {
	l = iN[i] + j;
	for (k = 0; k < 2; k++) {
	  m = iC[i] + k;
	  r2 = ((x[l] - x[m]) * (x[l] - x[m]) +
		(y[l] - y[m]) * (y[l] - y[m]) +
		(z[l] - z[m]) * (z[l] - z[m]));
	  e += q[j + 2] * q[k] / sqrt(r2);
	}
      }
    }

    (*ebi) = kbias * e;

    return ;
  }


  if (!strcmp(move,"piv") || !strcmp(move,"bgs")) {

    for (i = a1; i <= a2; i++) {
      for (j = 0; j < 2; j++) {
	l = iN[i] + j;
	for (k = 0; k < 2; k++) {
	  m = iC[i] + k;
	  r2 = ((x[l] - x[m]) * (x[l] - x[m]) +
		(y[l] - y[m]) * (y[l] - y[m]) +
		(z[l] - z[m]) * (z[l] - z[m]));
	  e += q[j + 2] * q[k] / sqrt(r2);
	}
      }
    }
    
    if (!strcmp(mc,"before")) 
      (*ebi) -= kbias * e;

    if (!strcmp(mc,"after")) 
      (*ebi) += kbias * e;

    return ;
  }

}
/****************************************************************************/
double hp_ecalc(int iflag,int l,int m) {
  int li,mi;
  double dx,dy,dz,r2,dr,e = 0;
  static double aco,bco,cuthp2;             

  if (iflag < 0) {
    cuthp2 = (sighp+dchp)*(sighp+dchp);
    aco = -(1+dchp*dchp)*exp(-dchp*dchp/2);
    bco = dchp*exp(-dchp*dchp/2);

    printf("hp: sighp %f dchp %f\n",sighp,dchp);
    printf("hp: aco %e bco %e\n",aco,bco);

    return 0;
  }

  if (a2c[l] == a2c[m] && abs(m - l) < 3)
    return 0;
  
  li = iCa[l] + 2; 
  mi = iCa[m] + 2;
  
  dx = x[li] - x[mi]; bc(&dx); 
  dy = y[li] - y[mi]; bc(&dy);
  dz = z[li] - z[mi]; bc(&dz);
  
  if ((r2 = dx * dx + dy * dy + dz * dz) > cuthp2) return 0;
  
  dr = sqrt(r2) - sighp;
  e -= exp(- dr * dr / 2) + aco + bco * dr;
  
  return e;
}
/****************************************************************************/
void hp(char *move,char *mc,int a1,int a2,double *ehp)  {
  int i,j,ip,hp0,hp1;
  double e = 0;
  static int lhp[N],nhp;

  if (NCH == 0) return ;

  if (!strcmp(move,"init")) {
    nhp = 0;
    for (i = 0; i < N; i++) 
      if (seq[i] == 'L') lhp[nhp++] = i;
    printf("<hp> nhp %i\n",nhp);
    hp_ecalc(-1,0,0); 
    return ;
  }

  if (!strcmp(move,"ecalc")) {

    for (i = 0; i < nhp; i++) {
      ip = lhp[i];
      for (j = i + 1; j < nhp; j++) 
	e += hp_ecalc(0,ip,lhp[j]);
    } 

    (*ehp) = epshp * e;

    return ;
  }

  if (!strcmp(move,"piv")) {

    if (a1 > a2) return ;

    hp0 = 0; while (hp0 <= nhp - 1 && lhp[hp0] < a1) ++hp0; // first hp index with lhp >= a1
    hp1 = nhp - 1; while (hp1 >= 0 && lhp[hp1] > a2) --hp1; // last hp index with lhp <= a2

    for (i = hp0; i <= hp1; i++) {
      ip = lhp[i];
      
      for (j = 0; j < hp0; j++) 
	e += hp_ecalc(0,ip,lhp[j]);
      
      for (j = hp1 + 1; j < nhp; j++) 
	e += hp_ecalc(0,ip,lhp[j]);
    } 

    if (!strcmp(mc,"before"))
      (*ehp) -= epshp * e;
    
    if (!strcmp(mc,"after"))
      (*ehp) += epshp * e;
    
    return ;
  }

  if (!strcmp(move,"bgs")) {
    
    hp0 = 0; while (hp0 <= nhp - 1 && lhp[hp0] < a1) ++hp0;
    hp1 = nhp - 1; while (hp1 >= 0 && lhp[hp1] > a2) --hp1;

    for (i = hp0; i <= hp1; i++) {
      ip = lhp[i];
      
      for (j = 0; j < i; j++) 
	e += hp_ecalc(0,ip,lhp[j]);
      
      for (j = hp1 + 1; j < nhp; j++) 
	e += hp_ecalc(0,ip,lhp[j]);
    }
    
    if (!strcmp(mc,"before"))
      (*ehp) -= epshp * e;
    
    if (!strcmp(mc,"after"))
      (*ehp) += epshp * e;
    
    return ;
  }

  if (!strcmp(move,"trans") || !strcmp(move,"rotate") ) {

    if (!strcmp(mc,"before"))
      (*ehp) -= Ehp;

    if (!strcmp(mc,"after")) {
      hp("ecalc","",0,0,&e);
      (*ehp) += e;
    }

    return ;
  }
}
/****************************************************************************/
double hbonds_ecalc(int iflag,int i,int j,double r2)
{
  int d1,d2,a1,a2,ia,ja,c1,c2;
  double r4,r6;
  double dx,dy,dz,dxd,dyd,dzd,dxa,dya,dza;
  double ca,cb;
  static double ahb,bhb,sighb2,cut2;  
  static double halfpowa,halfpowb;
  static double cst;

  if (iflag < 0) {
    cut2=cuthb*cuthb;
    bhb=-30*(pow(sighb/cuthb,10.0)-pow(sighb/cuthb,12.0))/cut2;
    ahb=-(5*pow(sighb/cuthb,12.0)-6*pow(sighb/cuthb,10.0))-bhb*cut2;
    sighb2=sighb*sighb;
    halfpowa=0.5*POWA;
    halfpowb=0.5*POWB;
    cst=pow(1.0/bNH,POWA)*pow(1.0/bCO,POWB);
    return 0.0;
  }

  d2 = i; d1 = i - 1;
  a1 = j; a2 = j - 1;
  ia = i2a[i]; c1 = a2c[ia];
  ja = i2a[j]; c2 = a2c[ja];
  //  d1 = iN[i]; d2 = d1 + 1; 
  //  a2 = iC[j]; a1 = a2 + 1;

  if (c1 == c2 && ja >= ia - 2 && ja <= ia + 1)
    return 0; 

  dx = x[a1] - x[d2]; bc(&dx); 
  dy = y[a1] - y[d2]; bc(&dy);
  dz = z[a1] - z[d2]; bc(&dz);
  if ( (r2 = dx*dx+dy*dy+dz*dz) > cut2) return 0;

  dxd = x[d2]-x[d1];
  dyd = y[d2]-y[d1];
  dzd = z[d2]-z[d1];
  if (POWA>0 && (ca=dxd*dx+dyd*dy+dzd*dz)<0) return 0;

  dxa = x[a2]-x[a1];
  dya = y[a2]-y[a1];
  dza = z[a2]-z[a1];
  if (POWB>0 && (cb=dxa*dx+dya*dy+dza*dz)<0) return 0; 

  r6 = sighb2 / r2; r4 = r6 * r6; r6 = r6 * r4;   

  return cst * gam[ia][ja]*pow(ca*ca/r2,halfpowa)*pow(cb*cb/r2,halfpowb)*
    (r6*(5*r6-6*r4)+ahb+bhb*r2);
}
/****************************************************************************/
double exvol_ecalc(int iflag,int a,double r2) {
  int i,j;
  double r6;
  static double sig2[NA*NA],asa[NA*NA],bsa[NA*NA];
  static double cutg2;
  
  if (iflag < 0) {
    cutg2 = ALPHA1 * ALPHA1 * cut * cut;
    for (i = 0; i < NA; i++) {
      for (j = 0; j < NA; j++) {
	if (i == 4 && j == 4)
	  sig2[i+j*NA] = pow(ALPHA2*(sigsa[i]+sigsa[j]),2.0);
	else if (i == 2 && j == 2)
	  sig2[i+j*NA]=pow(sigsa[i]+sigsa[j],2.0);
	else
	  sig2[i+j*NA]=pow(ALPHA1*(sigsa[i]+sigsa[j]),2.0);
      }
    }
    sig2[5+5*NA] = pow(ALPHA1*sighp,2.0);

    for (i = 0; i < NA*NA; i++) {
      asa[i] = -7*pow(sig2[i]/cutg2,6.0);
      bsa[i] = 6*pow(sig2[i]/cutg2,6.0)/cutg2;
    }

    exvol_ecalc(1,0,0);
    
    return 0;
  }

  if (iflag == 0) {    
    r6 = sig2[a] / r2; r6 = r6 * r6 * r6;
    return r6 * r6 + asa[a] + bsa[a] * r2;
  }

  if (iflag > 0) {
    double r, del = sqrt(cutg2) / 200;
    FILE *fp;

    fp = fopen("_exvol_ecalc.pot","w");
    for (r = 0; r <= sqrt(cutg2); r += del) 
      fprintf(fp,"%lf %lf\n",r,exvol_ecalc(0,0,r*r));
    fclose(fp);

    return 0;
  }

  return 0;
}
/****************************************************************************/
/* Cell division */
/****************************************************************************/
void cellcalc(char *move,char *mc,int icr,int a1,int a2,double *eev,double *ehb,
	      double *ecc,double *ecp,double emax) {
  int i,j,ic;
  int in_cell,nl,nlc,a,an;
  short ix,iy,iz,list[NTO],listc[NTO];
  double et1 = 0,et2 = 0,et3 = 0, et4 = 0;
  double dx,dy,dz,r2;
  static short cell[MAXCELL];    /* division into cells */
  static short cell2[MAXCELL];   /* division into cells */
  static short lx[NTO],ly[NTO],lz[NTO];
  static int lc[NTO],nec;
  static int lc2[NTO],nec2;
  static int pnt[NTO];           /* pointers to atoms */
  static int pnt2[NTO];          /* pointers to atoms */
  static int pnto[NTO];        
  static int ns;
  static double cutg;

  /* crowders */
  
  static double mcp[NCR],mcp_new[NCR];
  
  if (NCH == 0) return ;
  
  if (!strcmp(move,"init")) { /* init */
    nec = 0;    
    for (i = 0; i < MAXCELL; i++) cell[i] = cell2[i] = -1;
    for (i = 0; i < NTO; i++) pnt[i] = pnt2[i] = pnto[i] = -1;

    ns = BOX / (ALPHA1*cut);
    cutg = BOX / (double) ns;

    listcalc(-1,0,0,list,&et1,&et2);
    listcalc_crowd(-1,0,list,0,list,mcp,&et1,0);

    if (ns*ns*ns > MAXCELL) {
      printf("# cells > MAXCELL\n");
      exit(-1);
    } else
      printf("# cells %i ns %i cutg %f\n",ns*ns*ns,ns,cutg);

    return ;
  }

  if (!strcmp(move,"ecalc")) { /* full energy calc */
    for (i = 0; i < NCR; i++) mcp[i] = 0;
    for (i = 0; i < NTO; i++) if ((a = pnt[i]) < -1) cell[-a-2] = -1;
    nec = 0;

    et4 = 0;
    for (i = 0; i < NCR; i++) {
      for (j = 0; j < i; j++) {
	dx = xcr[i] - xcr[j]; bc(&dx); 
	dy = ycr[i] - ycr[j]; bc(&dy); 
	dz = zcr[i] - zcr[j]; bc(&dz); 
	r2 = dx * dx + dy * dy + dz * dz;

	if (r2 > cut_crowd_crowd2) continue;

	et4 += crowd_crowd_ecalc(0,r2);
      }
    }
    (*ecc) = epsilonrep * et4;
    
    in2box(xb,yb,zb,x,y,z,0,NTO-1);
    for (i = NTO-1; i >= 0; i--) {
      ix = xb[i] / cutg; iy = yb[i] / cutg; iz = zb[i] / cutg;
      ic = ix + ns * (iy + ns * iz);
      pnt[i] = cell[ic];
      if (cell[ic] < 0) {
	lx[nec] = ix; ly[nec] = iy; lz[nec] = iz; lc[nec++] = ic;
	pnt[i] = -ic-2;
      }
      cell[ic] = i;
    }

    i = 0;
    et1 = et2 = et3 = et4 = 0;
    while (i < nec) {
      nl = 0;
      build_list_inc(cell[(ic = lc[i])],pnt,&nl,list,0,0,0);
      in_cell = nl;

      build_list_fwd(ic,pnt,cell,lx[i],ly[i],lz[i],&nl,list,ns,0,0,0);
      listcalc(0,in_cell,nl,list,&et1,&et2);

      if (NCR > 0) {
	nlc = 0;
	build_list_crowd(lx[i],ly[i],lz[i],&nlc,listc,cutg);
	listcalc_crowd(0,in_cell,list,nlc,listc,mcp,&et3,0);
      }
      
      i++;
    }

    (*eev) = ksa * et1;
    (*ehb) = epshb * et2;
    (*ecp) = epsilonrep * et3;
    
    return ;
  }

  if ( !strcmp(move,"crowd") ) {
    
    if ( !strcmp(mc,"before") ) {

      for (i = 0; i < NCR; i++) {
	if (icr == i) continue;
	dx = xcr[icr] - xcr[i]; bc(&dx); 
	dy = ycr[icr] - ycr[i]; bc(&dy); 
	dz = zcr[icr] - zcr[i]; bc(&dz); 
	r2 = dx * dx + dy * dy + dz * dz;

	if (r2 > cut_crowd_crowd2) continue;

	et1 += crowd_crowd_ecalc(0,r2);
      }
      (*ecc) -= epsilonrep * et1;

      (*ecp) -= epsilonrep * mcp[icr];

      return ;
    }

    if ( !strcmp(mc,"after") ) {

      for (i = 0; i < NCR; i++) {
	if (icr == i) continue;
	dx = xcr[icr] - xcr[i]; bc(&dx); 
	dy = ycr[icr] - ycr[i]; bc(&dy); 
	dz = zcr[icr] - zcr[i]; bc(&dz); 
	r2 = dx * dx + dy * dy + dz * dz;

	if (r2 > cut_crowd_crowd2) continue;

	et1 += crowd_crowd_ecalc(0,r2);
      }
      (*ecc) += epsilonrep * et1;
      
      if (relax == 0 && (*ecc) + (*ecp) > emax) return;

      mcp_new[icr] = 0;
      for (j = 0; j < NTO; j++) {
	dx = xcr[icr] - x[j]; bc(&dx); 
	dy = ycr[icr] - y[j]; bc(&dy); 
	dz = zcr[icr] - z[j]; bc(&dz); 
	r2 = dx * dx + dy * dy + dz * dz;

	if (r2 > cut_crowd_atom2) continue;
	
	mcp_new[icr] += crowd_atom_ecalc(0,atyp[j],r2);
	
	if (ISOFT > 0)
	  mcp_new[icr] += crowd_atom_attr_ecalc(0,atyp[j],seq[i2a[j]],r2);
      }

      (*ecp) += epsilonrep * mcp_new[icr]; 
    }
    
    if ( !strcmp(mc,"acc") ) {
      mcp[icr] = mcp_new[icr];
      return ;
    }

    if ( !strcmp(mc,"rej") ) {
      return ;
    }
    
    return ;
  }

  if ( !strcmp(move,"trans") || !strcmp(move,"rotate") ) {
    
    if ( !strcmp(mc,"before") ) {
      (*eev) -= Eev;
      (*ehb) -= Ehb;
      (*ecc) -= Ecc;
      (*ecp) -= Ecp;
      return ;
    }

    if ( !strcmp(mc,"after") ) {
      cellcalc("ecalc","",0,0,0,&et1,&et2,&et3,&et4,0);
      (*eev) += et1;
      (*ehb) += et2;
      (*ecc) += et3;
      (*ecp) += et4;
      return ;
    }

    return ;
  }

  if ( !strcmp(move,"bgs") ) {
    
    if ( !strcmp(mc,"before") ) {
      for (i = a1; i <= a2; i++) pnto[i] = -1;
      nec2 = 0;

      for (i = 0; i < NCR; i++) mcp_new[i] = mcp[i];
      
      i = a2;
      while (i >= a1) {
	a = cell[(ic = get_cell(i,pnt))]; 
	lc2[nec2++] = ic;
	while ((a = pnto[a] = pnt[a]) >= 0); 	
	while (i >= a1 && pnto[i] != -1) i--;
      }

      for (i = 0; i < nec2; i++) {
	ic = lc2[i];
	get_cell_xyz(&ix,&iy,&iz,ic,ns);
	
	nl = 0;
	build_list_inc(cell[ic],pnt,&nl,list,1,a1,a2);
	in_cell = nl;
	
	build_list_inc(cell[ic],pnt,&nl,list,-1,a1,a2); 
	build_list_fwd(ic,pnt,cell,ix,iy,iz,&nl,list,ns,0,0,0);
	build_list_bwd(ic,pnt,cell,ix,iy,iz,&nl,list,ns,-1,a1,a2);
	listcalc(0,in_cell,nl,list,&et1,&et2);

	if (NCR > 0) {
	  nlc = 0;
	  build_list_crowd(ix,iy,iz,&nlc,listc,cutg);
	  listcalc_crowd(1,in_cell,list,nlc,listc,mcp_new,&et3,0);
	}
      }

      (*eev) -= ksa * et1;
      (*ehb) -= epshb * et2;
      (*ecp) -= epsilonrep * et3;
      
      return ;
    }

    if ( !strcmp(mc,"after") ) {
      nec = 0;

      in2box(xb,yb,zb,x,y,z,a1,a2);

      for (i = a2; i >= a1; i--) {
	ix = xb[i] / cutg; iy = yb[i] / cutg; iz = zb[i] / cutg;
	ic = ix + ns * (iy + ns * iz);
	pnt2[i] = cell2[ic];
	if (cell2[ic] < 0) {
	  lx[nec] = ix; ly[nec] = iy; lz[nec] = iz; lc[nec++] = ic;
	  pnt2[i] = -ic-2;
	}
	cell2[ic] = i;
      }

      for (i = 0; i < nec; i++) {
	et1 = et2 = et3 = 0;
	ic = lc[i]; ix = lx[i]; iy = ly[i]; iz = lz[i]; 

	nl = 0;
	build_list_inc(cell2[ic],pnt2,&nl,list,0,0,0);
	in_cell = nl;	
	build_list_fwd(ic,pnt2,cell2,ix,iy,iz,&nl,list,ns,0,0,0);

	build_list_inc(cell[ic],pnt,&nl,list,-1,a1,a2);
	build_list_fwd(ic,pnt,cell,ix,iy,iz,&nl,list,ns,-1,a1,a2);
	build_list_bwd(ic,pnt,cell,ix,iy,iz,&nl,list,ns,-1,a1,a2);	
	listcalc(0,in_cell,nl,list,&et1,&et2);
	
	if (NCR > 0) {
	  nlc = 0;
	  build_list_crowd(ix,iy,iz,&nlc,listc,cutg);
	  listcalc_crowd(0,in_cell,list,nlc,listc,mcp_new,&et3,0);
	}
	
	(*eev) += ksa * et1;
	(*ehb) += epshb * et2;
	(*ecp) += epsilonrep * et3;
	
	if (relax == 0 && (*eev) + (*ehb) + (*ecp) > emax) return;
      }

      return ;
    }

    if ( !strcmp(mc,"acc") ) {

      for (i = 0; i < NCR; i++)
	mcp[i] = mcp_new[i];
      
      for (i = 0; i < nec2; i++) {
        a = cell[(ic = lc2[i])];
	while ((an = a) >= 0) {
	  a = pnt[a];
	  if (an >= a1 && an <= a2)
	    remove_from_list(ic,pnt,cell,an);
	}
      }

      for (i = 0; i < nec; i++) {
	a = cell2[(ic = lc[i])];
	while ((an = a) >= 0) {
	  a = pnt2[a];
	  insert_into_list(ic,pnt,cell,an);
	} 
      }
      
      for (i = 0; i < nec; i++)
	cell2[lc[i]] = -1;

      return ;
    }

    if ( !strcmp(mc,"rej") ) {
      for (i = 0; i < nec; i++)
	cell2[lc[i]] = -1;

      return ;
    }

    return ;
  }

  if ( !strcmp(move,"piv") ) {

    if ( !strcmp(mc,"before") ) {
      for (i = a1; i <= a2; i++) pnto[i] = -1;
      nec2 = 0;
      
      for (i = 0; i < NCR; i++) mcp_new[i] = mcp[i];

      i = a2;
      while (i >= a1) {
	a = cell[(ic = get_cell(i,pnt))]; 
	lc2[nec2++] = ic;
	while ((a = pnto[a] = pnt[a]) >= 0); 	
	while (i >= a1 && pnto[i] != -1) i--;
      }

      for (i = 0; i < nec2; i++) {
	ic = lc2[i];
	get_cell_xyz(&ix,&iy,&iz,ic,ns);

	nl = 0;
	build_list_inc(cell[ic],pnt,&nl,list,1,a1,a2);
	in_cell = nl;
	
	build_list_inc(cell[ic],pnt,&nl,list,-1,a1,a2); 
	build_list_fwd(ic,pnt,cell,ix,iy,iz,&nl,list,ns,-1,a1,a2);
	build_list_bwd(ic,pnt,cell,ix,iy,iz,&nl,list,ns,-1,a1,a2);
	listcalc(1,in_cell,nl,list,&et1,&et2);

	if (NCR > 0) {
	  nlc = 0;
	  build_list_crowd(ix,iy,iz,&nlc,listc,cutg);
	  listcalc_crowd(1,in_cell,list,nlc,listc,mcp_new,&et3,0);
	}
      }

      (*eev) -= ksa * et1;
      (*ehb) -= epshb * et2;
      (*ecp) -= epsilonrep * et3;
      
      return ;
    }

    if ( !strcmp(mc,"after") ) {

      in2box(xb,yb,zb,x,y,z,a1,a2);

      nec = 0;
      for (i = a2; i >= a1; i--) {
	ix = xb[i] / cutg; iy = yb[i] / cutg; iz = zb[i] / cutg;
	ic = ix + ns * (iy + ns * iz);
	pnt2[i] = cell2[ic];
	if (cell2[ic] < 0) {
	  lx[nec] = ix; ly[nec] = iy; lz[nec] = iz; lc[nec++] = ic;
	  pnt2[i] = -ic-2;
	}
	cell2[ic] = i;
      }
      
      for (i = 0; i < nec; i++) {
	et1 = et2 = et3 = 0;
	ic = lc[i]; ix = lx[i]; iy = ly[i]; iz = lz[i]; 
	
	nl = 0;
	build_list_inc(cell2[ic],pnt2,&nl,list,1,a1,a2);
	in_cell = nl;
	
	build_list_inc(cell[ic],pnt,&nl,list,-1,a1,a2); 
	build_list_fwd(ic,pnt,cell,ix,iy,iz,&nl,list,ns,-1,a1,a2);
	build_list_bwd(ic,pnt,cell,ix,iy,iz,&nl,list,ns,-1,a1,a2);	
	listcalc(1,in_cell,nl,list,&et1,&et2);

	if (NCR > 0) {
	  nlc = 0;
	  build_list_crowd(ix,iy,iz,&nlc,listc,cutg);
	  listcalc_crowd(0,in_cell,list,nlc,listc,mcp_new,&et3,0);
	}
	
	(*eev) += ksa * et1;
	(*ehb) += epshb * et2;
	(*ecp) += epsilonrep * et3;

	if (relax == 0 && (*eev) + (*ehb) + (*ecp) > emax) return;
      }
      
      return ;
    }

    if ( !strcmp(mc,"acc") ) {
      for (i = 0; i < NCR; i++) mcp[i] = mcp_new[i];
      
      for (i = 0; i < nec2; i++) {
        a = cell[(ic = lc2[i])];
	while ((an = a) >= 0) {
	  a = pnt[a];
	  if (an >= a1 && an <= a2)
	    remove_from_list(ic,pnt,cell,an);
	}
      }

      for (i = 0; i < nec; i++) {
	a = cell2[(ic = lc[i])];
	while ((an = a) >= 0) {
	  a = pnt2[a];
	  insert_into_list(ic,pnt,cell,an);
	} 
      }
      for (i = 0; i < nec; i++) cell2[lc[i]] = -1;
      return ;
    }

    if ( !strcmp(mc,"rej") ) {
      for (i = 0; i < nec; i++) cell2[lc[i]] = -1;
      return ;
    }

    return ;
  }

}
/****************************************************************************/
void listcalc(int iflag,int in_cell,int nl,short list[],double *eev,double *ehb) {
  int i,j,ihb,jhb,n,m,a;
  double r2,dx,dy,dz;
  static double cutg2;

  if (iflag < 0) {
    cutg2 = ALPHA1 * ALPHA1 * cut * cut;
    exvol_ecalc(-1,0,0);
    hbonds_ecalc(-1,0,0,0);
    return ;
  }
  
  if (iflag == 0) {

    for (n = 0; n < in_cell; n++) {
      i = list[n]; ihb = idonacc[i];
      
      for (m = n + 1; m < nl; m++) {
	j = list[m]; jhb = idonacc[j];

	if ((a = aa[i][j]) < 0 || a >= NA*NA) continue;
	
	dx = x[i] - x[j]; bc(&dx); 
	dy = y[i] - y[j]; bc(&dy); 
	dz = z[i] - z[j]; bc(&dz); 
	if ( (r2 = dx * dx + dy * dy + dz * dz) > cutg2) continue;
	
	(*eev) += exvol_ecalc(0,a,r2);
	
	if (ihb > 0 && jhb < 0)
	  (*ehb) += hbonds_ecalc(0,i,j,r2);
	else if (ihb < 0 && jhb > 0)
	  (*ehb) += hbonds_ecalc(0,j,i,r2);
      }
    }

    return ;
  }
    
  if (iflag == 1) {
    
    for (n = 0; n < in_cell; n++) {
      i = list[n]; ihb = idonacc[i];
      
      for (m = in_cell; m < nl; m++) {
	j = list[m]; jhb = idonacc[j];
	
	if ((a = aa[i][j]) < 0 || a >= NA*NA) continue;
	
	dx = x[i] - x[j]; bc(&dx); 
	dy = y[i] - y[j]; bc(&dy); 
	dz = z[i] - z[j]; bc(&dz); 
	if ( (r2 = dx * dx + dy * dy + dz * dz) > cutg2) continue;

	(*eev) += exvol_ecalc(0,a,r2);
	
	if (ihb > 0 && jhb < 0)
	  (*ehb) += hbonds_ecalc(0,i,j,r2);
	else if (ihb < 0 && jhb > 0)
	  (*ehb) += hbonds_ecalc(0,j,i,r2);
      }

    }

    return ;
  }

}
/****************************************************************************/
void listcalc_crowd(int iflag,int nl,short list[],int nlc,short listc[],
		    double mcp[],double *ecp,double emax) {
  int i,j,m,n,at;
  double dx,dy,dz,r2;
  double et1,et2 = 0;
  
  if (iflag < 0) {
    crowd_crowd_ecalc(-1,0);
    crowd_atom_ecalc(-1,0,0);
    if (ISOFT) crowd_atom_attr_ecalc(-1,0,0,0);
    return ;
  }

  if (iflag == 0) {

    for (n = 0; n < nl; n++) {
      i = list[n]; at = atyp[i];

      for (m = 0; m < nlc; m++) {
	j = listc[m];
	
	dx = xcr[j] - x[i]; bc(&dx); 
	dy = ycr[j] - y[i]; bc(&dy); 
	dz = zcr[j] - z[i]; bc(&dz); 
	r2 = dx * dx + dy * dy + dz * dz;

	if (r2 > cut_crowd_atom2) continue;
	
	et2 += (et1 = crowd_atom_ecalc(0,at,r2));
	mcp[j] += et1;

	if (ISOFT > 0) {
	  et2 += (et1 = crowd_atom_attr_ecalc(0,at,seq[i2a[i]],r2));
	  mcp[j] += et1;
	}
	
      }
    }
    (*ecp) += et2;
    
    return ;
  }

  if (iflag == 1) {

    for (n = 0; n < nl; n++) {
      i = list[n]; at = atyp[i];
      
      for (m = 0; m < nlc; m++) {
	j = listc[m];

	dx = xcr[j] - x[i]; bc(&dx); 
	dy = ycr[j] - y[i]; bc(&dy); 
	dz = zcr[j] - z[i]; bc(&dz); 
	r2 = dx * dx + dy * dy + dz * dz;

	if (r2 > cut_crowd_atom2) continue;

	et2 += (et1 = crowd_atom_ecalc(0,at,r2));
	mcp[j] -= et1;

	if (ISOFT > 0) {
	  et2 += (et1 = crowd_atom_attr_ecalc(0,at,seq[i2a[i]],r2));
	  mcp[j] -= et1;
	}

      }
    } 
   (*ecp) += et2;

    return ;
  }
  
}
/****************************************************************************/
/***** UPDATES CROWDERS *****************************************************/
/****************************************************************************/
int crowd_translate() {
  int icr;
  double dx,dy,dz;
  double de = 0,de_cc = 0,de_ca = 0,de_acc;
  double dummy;
  
  if (NCR == 0) return 0;
  
  /* pick a crowder */

  icr = NCR * genrand();

  /* energy calc */

  cellcalc("crowd","before",icr,0,0,&dummy,&dummy,&de_cc,&de_ca,0);

  /* pick a direction */
  
  dx = CRSTEP * 2 * (.5 - genrand()); 
  dy = CRSTEP * 2 * (.5 - genrand()); 
  dz = CRSTEP * 2 * (.5 - genrand()); 

  /* translate */
  
  xcr[icr] += dx;
  ycr[icr] += dy;
  zcr[icr] += dz;

  de_acc = - log(genrand()) / beta[ind];
  
  cellcalc("crowd","after",icr,0,0,&dummy,&dummy,&de_cc,&de_ca,e_max);
  de += de_cc + de_ca;
  
  if (de < de_acc) {
    /* accept */    
    xcro[icr] = xcr[icr];
    ycro[icr] = ycr[icr];
    zcro[icr] = zcr[icr];
    E = ( (Eloc + Ebias + Ehp + Ehb + Eev) +
	  (Ecc   += de_cc) +
	  (Ecp   += de_ca) );
    cellcalc("crowd","acc",icr,0,0,&dummy,&dummy,&de_cc,&de_ca,0);
    return 1;
  } else {
    /* reject */
    xcr[icr] = xcro[icr];
    ycr[icr] = ycro[icr];
    zcr[icr] = zcro[icr];

    cellcalc("crowd","rej",icr,0,0,&dummy,&dummy,&de_cc,&de_ca,0);
    return 0;
  }
}
/****************************************************************************/
/***** UPDATES CHAINS *******************************************************/
/****************************************************************************/
int sel_aa2(int iflag,int a1,int a2) {
 
  if (NCH == 0) return 0;

  if (iflag < 0) {
    printf("init sel_aa...\n");
    printf("    DONE sel_aa.\n");
    return -1;
  }

  return a1 + (a2 - a1 + 1) * genrand();
}
/****************************************************************************/
int sel_aa(char *mov) {
  int i,j,n,len,mid;
  static int ibgsup[N],ibgsdw[N],ipivup[N],ipivdw[N];
  static int nbgsup,nbgsdw,npivup,npivdw;
  FILE *fp;

 
  if (NCH == 0) return 0;

  if (!strcmp(mov,"init")) {
    printf("init sel_aa...\n");
    nbgsup = nbgsdw = npivup = npivdw = 0;
    
    for (i = 0; i < NCH; ++i) {
      len = Ne[i] - Nb[i] + 1;
      mid = Nb[i] + len / 2;
      for (j = Nb[i]; j < mid; ++j)  ipivdw[npivdw++] = j;
      for (j = mid; j <= Ne[i]; ++j) ipivup[npivup++] = j;
      for (j = Nb[i] + 3; j < mid; ++j)  ibgsdw[nbgsdw++] = j;
      for (j = mid; j <= Ne[i] - 3; ++j) ibgsup[nbgsup++] = j;
    } 
    
    fp = fopen("_updates_sel_aa","w");
    fprintf(fp,"nbgsup %i nbgsdw %i \n",nbgsup,nbgsdw);
    fprintf(fp,"npivup %i npivdw %i \n\n",npivup,npivdw);
    for (i = 0; i < nbgsup; ++i) fprintf(fp,"bgsup: ch %i iaa %i \n",a2c[ibgsup[i]],ibgsup[i]);
    for (i = 0; i < nbgsdw; ++i) fprintf(fp,"bgsdw: ch %i iaa %i \n",a2c[ibgsdw[i]],ibgsdw[i]);
    for (i = 0; i < npivup; ++i) fprintf(fp,"pivup: ch %i iaa %i \n",a2c[ipivup[i]],ipivup[i]);
    for (i = 0; i < npivdw; ++i) fprintf(fp,"pivdw: ch %i iaa %i \n",a2c[ipivdw[i]],ipivdw[i]);
    fclose(fp);

    if (nbgsup == 0) printf("*** nbgsup == 0 ***\n");
    if (nbgsdw == 0) printf("*** nbgsdw == 0 ***\n");

    printf("    DONE sel_aa.\n");
    return -1;
  }

  if (!strcmp(mov,"bgsup")) {
    n = nbgsup * genrand();
    //    printf("bgsup %s %2i %2i\n)",mov,n,ibgsup[n]);
    return ibgsup[n];
  }
  
  if (!strcmp(mov,"bgsdw")) {
    n = nbgsdw * genrand();
    //    printf("bgsdw %s %2i %2i\n)",mov,n,ibgsdw[n]);
    return ibgsdw[n];
  }
  
  if (!strcmp(mov,"pivup")) {
    n = npivup * genrand();
    //    printf("pivup %s %2i %2i\n)",mov,n,ipivup[n]);
    return ipivup[n];
  }
  
  if (!strcmp(mov,"pivdw")) {
    n = npivdw * genrand();
    //    printf("pivdw %s %2i %2i\n)",mov,n,ipivdw[n]);
    return ipivdw[n];
  }
  
  return 0;
}
/****************************************************************************/
void qsi(long *nup,double *accup,long *ndw,double *accdw) {
  
  if (genrand() < 0.5) {
    (*nup)++; (*accup) += qsiup(0,sel_aa("bgsup"));
  } else {
    (*ndw)++; (*accdw) += qsidw(0,sel_aa("bgsdw"));
  }
  
}
/****************************************************************************/
int qsiup(int iflag,int ia1)
{
  int i,k,i1,i2,ic,ia2;
  double dph[8];
  double de = 0,de_cc = 0,de_ca = 0;
  double de_exvol = 0,de_hbonds = 0;
  double de_local = 0,de_bias = 0,de_hp = 0;
  
  if (NCH == 0) return 0;
  
  if (iflag < 0) {
    /* init */
    return -1;
  }

  ic = a2c[ia1];
  ia2 = ia1 + 3;

  i1 = iCa[ia1] + 1;
  i2 = iEnd[ic];

  /* energy calc */

  cellcalc("bgs","before",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
  localsa("bgs","before",ia1,ia2,&de_local);
  bias("bgs","before",ia1,ia2,&de_bias);
  hp("bgs","before",ia1,Ne[ic],&de_hp);
  
  /* move */

  for (k = 0, i = ia1; i <= ia2; i++) {
    ph[3 * i]     += ( dph[k++] = (2. / 360.) * pi2 * 2. * (.5 - genrand()) ); 
    ph[3 * i + 1] += ( dph[k++] = (2. / 360.) * pi2 * 2. * (.5 - genrand()) );
  }

  //  printf("qsiup: ia1 %d ia2 %d %lf %lf\n",ia1,ia2,dph[0],dph[1]);
  
  int2cart(ia1,1);

  /* energy calc */

  cellcalc("bgs","after",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,e_max);
  de += de_exvol + de_hbonds + de_cc + de_ca;

  if (de < e_max) {
    localsa("bgs","after",ia1,ia2,&de_local);
    bias("bgs","after",ia1,ia2,&de_bias);
    hp("bgs","after",ia1,Ne[ic],&de_hp);
    de += de_local + de_bias + de_hp;
  }

  if ( de < -log(genrand()) / beta[ind] ) {
    for (i = i1; i <= i2; i++) {
      xo[i] = x[i]; 
      yo[i] = y[i]; 
      zo[i] = z[i];
    }

    E = ( (Ecc   += de_cc) +
	  (Ecp   += de_ca) +
	  (Eev   += de_exvol) + 
	  (Ehb   += de_hbonds) +
	  (Eloc  += de_local) + 
	  (Ebias += de_bias) +
	  (Ehp   += de_hp) );
    
    cellcalc("bgs","acc",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);

    return 1;
  }
  else {
    for (k = 0, i = ia1; i <= ia2; i++) {
      ph[3 * i]     -= dph[k++]; 
      ph[3 * i + 1] -= dph[k++]; 
    }
    for (i = i1; i <= i2; i++) {
      x[i] = xo[i]; 
      y[i] = yo[i]; 
      z[i] = zo[i];
    }
    cellcalc("bgs","rej",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    return 0;
  }
}
/****************************************************************************/
int qsidw(int iflag,int ia2)
{
  int i,k,i1,i2,ic,ia1;
  double dph[8];
  double de = 0,de_cc = 0,de_ca = 0;
  double de_exvol = 0,de_hbonds = 0;
  double de_local = 0,de_bias = 0,de_hp = 0;
  
  if (NCH == 0) return 0;
  
  if (iflag < 0) {
    return -1;
  }

  ic = a2c[ia2];
  ia1 = ia2 - 3;

  i1 = iBeg[ic];
  i2 = iC[ia2] - 1;

  /* energy calc */

  cellcalc("bgs","before",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
  localsa("bgs","before",ia1,ia2,&de_local);
  bias("bgs","before",ia1,ia2,&de_bias);
  hp("bgs","before",Nb[ic],ia2,&de_hp);
  
  /* move */

  for (k = 0, i = ia1; i <= ia2; i++) {
    ph[3 * i]     += ( dph[k++] = (2. / 360.) * pi2 * 2. * (.5 - genrand()) ); 
    ph[3 * i + 1] += ( dph[k++] = (2. / 360.) * pi2 * 2. * (.5 - genrand()) ); 
  }

  int2cart(ia2,-1);

  /* energy calc */

  cellcalc("bgs","after",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,e_max);
  de += de_exvol + de_hbonds + de_cc + de_ca;

  if (de < e_max) {
    localsa("bgs","after",ia1,ia2,&de_local);
    bias("bgs","after",ia1,ia2,&de_bias);
    hp("bgs","after",Nb[ic],ia2,&de_hp);
    de += de_local + de_bias + de_hp;
  }
  
  if ( de < -log(genrand()) / beta[ind] ) {
    for (i = i1; i <= i2; i++) {
      xo[i] = x[i]; 
      yo[i] = y[i]; 
      zo[i] = z[i];
    }

    E = ( (Ecc   += de_cc) +
	  (Ecp   += de_ca) +
	  (Eev   += de_exvol) + 
	  (Ehb   += de_hbonds) +
	  (Eloc  += de_local) + 
	  (Ebias += de_bias) +
	  (Ehp   += de_hp) );
    
    cellcalc("bgs","acc",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);

    return 1;
  }
  else {
    for (k = 0, i = ia1; i <= ia2; i++) {
      ph[3 * i]     -= dph[k++]; 
      ph[3 * i + 1] -= dph[k++]; 
    }
    for (i = i1;i <= i2;i++) {
      x[i] = xo[i]; 
      y[i] = yo[i]; 
      z[i] = zo[i];
    }
    cellcalc("bgs","rej",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    return 0;
  }
}
/****************************************************************************/
void bgs(int iflag,long *nup,double *accup,long *ndw,double *accdw) {
  
  int i,ia,acc,ia0,ic,n;
  static long int natt[N],nacc[N];
  static int nmov,imov[N];
  FILE *fp;

  if (iflag < 0) {
    printf("init bgs\n");
    fp = fopen("_bgs_acc","w");
    fclose(fp);

    nmov = 0;
    for (i = 0; i < N; i++) {
      ic = a2c[i];
      if (i >= Nb[ic] + 3 && i <= Ne[ic] - 3)
	imov[nmov++] = i;
    }
    for (i = 0; i < N; i++) natt[i] = nacc[i] = 0;
    printf("nmov = %i \n",nmov);
    printf("  .. DONE\n");
    return ;
  }

  if (iflag == 0) {
    if (nmov == 0)
      return ;
    
    n = nmov * genrand();
    ia = imov[n];
    ic = a2c[ia]; 

    ia0 = Nb[ic] + (Ne[ic] - Nb[ic] + 1) / 2;

    if (ia > ia0 + 1) {
      (*nup)++; (*accup) += ( acc = bgsup("mc",ia) );
    }

    if (ia < ia0 - 1) {
      (*ndw)++; (*accdw) += ( acc = bgsdw("mc",ia) );
    }

    if (ia >= ia0 - 1 && ia <= ia0 + 1) {    
      if (genrand() < 0.5) {
	(*nup)++; (*accup) += ( acc = bgsup("mc",ia) );
      } else {
	(*ndw)++; (*accdw) += ( acc = bgsdw("mc",ia) );
      }
    }
    
    natt[ia] += 1;
    nacc[ia] += acc;

    return ;
  }

  if (iflag > 0) {
    fp = fopen("_bgs_acc","w");

    for (i = 0; i < N; i++) {
      fprintf(fp,"%i %2li %2li %2.8lf\n",i,natt[i],nacc[i],(double) nacc[i] / (double) natt[i]);
    }
    
    fclose(fp);

    return ;
  }

}
/****************************************************************************/
int bgsup(char *move,int ia1)
{
  int i,j,k,ia2,ic,i1,i2,iph[8],nph = 0;
  double de = 0,de_cc = 0,de_ca = 0,de_acc;
  double de_exvol = 0,de_hbonds = 0;
  double de_local = 0,de_bias = 0,de_hp = 0;
  double ix[8],iy[8],iz[8],bx[8],by[8],bz[8];
  double rx[3],ry[3],rz[3];
  double dx[3][8],dy[3][8],dz[3][8];
  double A[8][8],p[8],psi[8],dph[8],sum,r1,r2;
  double wfw,wbw;
  double renorm;
  
  ic = a2c[ia1];  
  ia2 = ia1 + 3;

  i1 = iCa[ia1] + 1;
  i2 = iEnd[ic];

  if (ia1 < 0 || ia1 > N-1 || ia2 < 0 || ia2 > N-1) {
    printf("bgsup ia1 %i ia2 %i\n",ia1,ia2);
    fflush(stdout);
    exit(-1);
  }
  
  /* energy calc */

  if (!strcmp(move,"mc")) {
    cellcalc("bgs","before",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    localsa("bgs","before",ia1,ia2,&de_local);
    bias("bgs","before",ia1,ia2,&de_bias);
    hp("bgs","before",ia1,Ne[ic],&de_hp);
  }

  /** forward move **/

  for (i = 0; i < 4; i++) {
    ix[nph] = x[iCa[ia1 + i]];
    iy[nph] = y[iCa[ia1 + i]];
    iz[nph] = z[iCa[ia1 + i]];
    bx[nph] = x[iCa[ia1 + i]] - x[iN[ia1 + i]];
    by[nph] = y[iCa[ia1 + i]] - y[iN[ia1 + i]];
    bz[nph] = z[iCa[ia1 + i]] - z[iN[ia1 + i]];
    renorm = sqrt(bx[nph]*bx[nph] + by[nph]*by[nph] + bz[nph]*bz[nph]);
    bx[nph] /= renorm; by[nph] /= renorm; bz[nph] /= renorm;
    iph[nph] = 3 * (ia1 + i);
    nph++;

    ix[nph] = x[iC[ia1 + i]];
    iy[nph] = y[iC[ia1 + i]];
    iz[nph] = z[iC[ia1 + i]];
    bx[nph] = x[iC[ia1 + i]] - x[iCa[ia1 + i]];
    by[nph] = y[iC[ia1 + i]] - y[iCa[ia1 + i]];
    bz[nph] = z[iC[ia1 + i]] - z[iCa[ia1 + i]];
    renorm = sqrt(bx[nph]*bx[nph] + by[nph]*by[nph] + bz[nph]*bz[nph]);
    bx[nph] /= renorm; by[nph] /= renorm; bz[nph] /= renorm;
    iph[nph] = 3 * (ia1 + i) + 1;
    nph++;
  }

  i = iCa[ia2];     rx[0] = x[i]; ry[0] = y[i]; rz[0] = z[i];
  i = iC[ia2];      rx[1] = x[i]; ry[1] = y[i]; rz[1] = z[i];
  i = iC[ia2] + 1;  rx[2] = x[i]; ry[2] = y[i]; rz[2] = z[i];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < nph; j++) {
      dx[i][j] = by[j] * (rz[i] - iz[j]) - bz[j] * (ry[i] - iy[j]) ;
      dy[i][j] = bz[j] * (rx[i] - ix[j]) - bx[j] * (rz[i] - iz[j]) ;
      dz[i][j] = bx[j] * (ry[i] - iy[j]) - by[j] * (rx[i] - ix[j]) ;
    }
  }

  /* Matrix A (symmetric, values in upper triangle incl diagonal) */

  for (i = 0; i < nph; i++) {
    for (j = i; j < nph; j++) {

      A[i][j] = 0;
      for (k = 0; k < 3; k++) {
	A[i][j] += dx[k][i] * dx[k][j] + dy[k][i] * dy[k][j] + dz[k][i] * dz[k][j];
      }
      A[i][j] *= bbgs;
      if (i == j) A[i][j] += 1;
      A[i][j] *= abgs / 2;

    }
  }

  /* Cholesky decomp of A=LL^T (put in lower triangle A and p[]) */

  for (i = 0; i < nph; i++) {
    for (j = i; j < nph; j++) {
      for (sum = A[i][j] , k = i - 1; k >= 0; k--) sum -= A[i][k] * A[j][k];
      if (i == j)
	p[i] = sqrt(sum);
      else
	A[j][i] = sum / p[i];
    }
  }

  /* Gaussian random numbers, psi */

  for (i = 0; i < 8; i += 2) {
    r1 = sqrt( -log(genrand()) );
    r2 = genrand();
    psi[i] = r1 * cos(pi2 * r2);
    psi[i + 1] = r1 * sin(pi2 * r2);
  }

  /* solve L^T * dph = psi */

  for (i = nph - 1; i >= 0; i--) {
    for (sum = psi[i] , k = i + 1; k < nph; k++) if (k < nph) sum -= A[k][i] * dph[k];
    dph[i] = sum / p[i];
  }
  
  //  printf("%lf %lf %lf %lf %lf %lf %lf %lf \n",
  // 	 dph[0],dph[1],dph[2],dph[3],dph[4],dph[5],dph[6],dph[7]);	 

  /* detailed balance */

  sum = 0; for (i = 0; i < nph; i++) sum += psi[i] * psi[i];
  wfw = exp(-sum); 
  for (i = 0; i < nph; i++) wfw *= p[i]; 

  /** backward probability **/

  for (i = 0; i < nph; i++) ph[iph[i]] += dph[i];
  int2cart(ia1,1);

  if (!strcmp(move,"move")) 
    return 1;
  
  nph = 0;
  for (i = 0; i < 4; i++){
    ix[nph] = x[iCa[ia1 + i]];
    iy[nph] = y[iCa[ia1 + i]];
    iz[nph] = z[iCa[ia1 + i]];
    bx[nph] = x[iCa[ia1 + i]] - x[iN[ia1 + i]];
    by[nph] = y[iCa[ia1 + i]] - y[iN[ia1 + i]];
    bz[nph] = z[iCa[ia1 + i]] - z[iN[ia1 + i]];
    renorm = sqrt(bx[nph]*bx[nph] + by[nph]*by[nph] + bz[nph]*bz[nph]);
    bx[nph] /= renorm; by[nph] /= renorm; bz[nph] /= renorm;
    nph++;

    ix[nph] = x[iC[ia1 + i]];
    iy[nph] = y[iC[ia1 + i]];
    iz[nph] = z[iC[ia1 + i]];
    bx[nph] = x[iC[ia1 + i]] - x[iCa[ia1 + i]];
    by[nph] = y[iC[ia1 + i]] - y[iCa[ia1 + i]];
    bz[nph] = z[iC[ia1 + i]] - z[iCa[ia1 + i]];
    renorm = sqrt(bx[nph]*bx[nph] + by[nph]*by[nph] + bz[nph]*bz[nph]);
    bx[nph] /= renorm; by[nph] /= renorm; bz[nph] /= renorm;
    nph++;
  } 

  i = iCa[ia2];      rx[0] = x[i]; ry[0] = y[i]; rz[0] = z[i];
  i = iC[ia2];       rx[1] = x[i]; ry[1] = y[i]; rz[1] = z[i];
  i = iC[ia2] + 1;   rx[2] = x[i]; ry[2] = y[i]; rz[2] = z[i];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < nph; j++) {
      dx[i][j] = by[j] * (rz[i]-iz[j]) - bz[j] * (ry[i]-iy[j]) ;
      dy[i][j] = bz[j] * (rx[i]-ix[j]) - bx[j] * (rz[i]-iz[j]) ;
      dz[i][j] = bx[j] * (ry[i]-iy[j]) - by[j] * (rx[i]-ix[j]) ;
    }
  }

  for (i=0;i<nph;i++) {
    for (j=i;j<nph;j++) {
      A[i][j]=0;
      for (k=0;k<3;k++) {
	A[i][j]+=dx[k][i]*dx[k][j]+dy[k][i]*dy[k][j]+dz[k][i]*dz[k][j];
      }
      A[i][j]*=bbgs;
      if (i==j) A[i][j]+=1;
      A[i][j]*=abgs/2;
    }
  }

  for (i=0;i<nph;i++) {
    for (j=i;j<nph;j++) {
      for (sum=A[i][j],k=i-1;k>=0;k--) sum-=A[i][k]*A[j][k];
      if (i==j)
	p[i]=sqrt(sum);
      else
	A[j][i]=sum/p[i];
    }
  }

  for (i=0;i<nph;i++) {
    psi[i] = p[i] * dph[i]; for (j=i + 1;j < nph; j++) if (j < nph) psi[i] += A[j][i] * dph[j];  
  }

  sum = 0; for (i = 0; i < nph; i++) sum += psi[i] * psi[i];
  wbw = exp(-sum); 
  for (i = 0; i < nph; i++) wbw *= p[i]; 			  

  /** accept/reject **/

  de_acc = - log( (wfw / wbw) * genrand() ) / beta[ind];

  cellcalc("bgs","after",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,e_max);
  de += de_exvol + de_hbonds + de_cc + de_ca;

  if (de < e_max) {
    localsa("bgs","after",ia1,ia2,&de_local);
    bias("bgs","after",ia1,ia2,&de_bias);
    hp("bgs","after",ia1,Ne[ic],&de_hp);
    de += de_local + de_bias + de_hp;
  }
  
  if (de < de_acc) {
    for (i = i1; i <= i2; i++) {
      xo[i] = x[i]; 
      yo[i] = y[i]; 
      zo[i] = z[i];
    }

    E = ( (Ecc   += de_cc) +
	  (Ecp   += de_ca) +
	  (Eev   += de_exvol) + 
	  (Ehb   += de_hbonds) +
	  (Eloc  += de_local) + 
	  (Ebias += de_bias) +
	  (Ehp   += de_hp) );

    cellcalc("bgs","acc",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    return 1;
  }
  else {
    for (i = 0; i < nph; i++) ph[iph[i]] -= dph[i];
    for (i = i1; i <= i2; i++) {
      x[i] = xo[i]; 
      y[i] = yo[i]; 
      z[i] = zo[i];
    }
    cellcalc("bgs","rej",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    return 0;
  }
}
/****************************************************************************/
int bgsdw(char *move,int ia2)
{
  int i,j,k,ia1,ic,i1,i2,iph[8],nph=0;
  double de = 0;
  double de_cc = 0,de_ca = 0,de_acc;
  double de_exvol = 0,de_hbonds = 0;
  double de_local = 0,de_bias = 0,de_hp = 0;
  double ix[8],iy[8],iz[8],bx[8],by[8],bz[8];
  double rx[3],ry[3],rz[3];
  double dx[3][8],dy[3][8],dz[3][8];
  double A[8][8],p[8],psi[8],dph[8],sum,r1,r2;
  double wfw,wbw;
  double renorm;
  
  if (NCH == 0) return 0;

  ic = a2c[ia2];
  ia1 = ia2 - 3;

  i1 = iBeg[ic];
  i2 = iC[ia2] - 1;


  if (ia1 < 0 || ia1 > N-1 || ia2 < 0 || ia2 > N-1) {
    printf("bgsdw ia1 %i ia2 %i\n",ia1,ia2);
    fflush(stdout);
    exit(-1);
  }


  /* energy calc */

  if (!strcmp(move,"mc")) {
    cellcalc("bgs","before",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    localsa("bgs","before",ia1,ia2,&de_local);
    bias("bgs","before",ia1,ia2,&de_bias);
    hp("bgs","before",Nb[ic],ia2,&de_hp);
  }

  /** forward move **/

  for (i = 0; i < 4; i++){
    ix[nph] = x[iCa[ia1 + i]];
    iy[nph] = y[iCa[ia1 + i]];
    iz[nph] = z[iCa[ia1 + i]];
    bx[nph] = x[iCa[ia1 + i]] - x[iN[ia1 + i]];
    by[nph] = y[iCa[ia1 + i]] - y[iN[ia1 + i]];
    bz[nph] = z[iCa[ia1 + i]] - z[iN[ia1 + i]];
    renorm = sqrt(bx[nph]*bx[nph] + by[nph]*by[nph] + bz[nph]*bz[nph]);
    bx[nph] /= renorm; by[nph] /= renorm; bz[nph] /= renorm;
    iph[nph] = 3 * (ia1 + i);
    nph++;

    ix[nph] = x[iC[ia1 + i]];
    iy[nph] = y[iC[ia1 + i]];
    iz[nph] = z[iC[ia1 + i]];
    bx[nph] = x[iC[ia1 + i]] - x[iCa[ia1 + i]];
    by[nph] = y[iC[ia1 + i]] - y[iCa[ia1 + i]];
    bz[nph] = z[iC[ia1 + i]] - z[iCa[ia1 + i]];
    renorm = sqrt(bx[nph]*bx[nph] + by[nph]*by[nph] + bz[nph]*bz[nph]);
    bx[nph] /= renorm; by[nph] /= renorm; bz[nph] /= renorm;
    iph[nph] = 3 * (ia1 + i) + 1;
    nph++;
  }
  
  i = iCa[ia2];     rx[0] = x[i]; ry[0] = y[i]; rz[0] = z[i];
  i = iC[ia2];      rx[1] = x[i]; ry[1] = y[i]; rz[1] = z[i];
  i = iC[ia2] + 1;  rx[2] = x[i]; ry[2] = y[i]; rz[2] = z[i];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < nph; j++) {
      dx[i][j] = by[j] * (rz[i] - iz[j]) - bz[j] * (ry[i] - iy[j]) ;
      dy[i][j] = bz[j] * (rx[i] - ix[j]) - bx[j] * (rz[i] - iz[j]) ;
      dz[i][j] = bx[j] * (ry[i] - iy[j]) - by[j] * (rx[i] - ix[j]) ;
    }
  }

  /* Matrix A (symmetric, values in upper triangle incl diagonal) */

  for (i = 0; i < nph; i++) {
    for (j = i; j < nph; j++) {

      A[i][j] = 0;
      for (k = 0; k < 3; k++) {
	A[i][j] += dx[k][i] * dx[k][j] + dy[k][i] * dy[k][j] + dz[k][i] * dz[k][j];
      }
      A[i][j] *= bbgs;
      if (i == j) A[i][j] += 1;
      A[i][j] *= abgs / 2;

    }
  }

  /* Cholesky decomp of A=LL^T (put in lower triangle A and p[]) */

  for (i = 0; i < nph; i++) {
    for (j = i; j < nph; j++) {
      for (sum = A[i][j] , k = i - 1; k >= 0; k--) sum -= A[i][k] * A[j][k];
      if (i == j)
	p[i] = sqrt(sum);
      else
	A[j][i] = sum / p[i];
    }
  }

  /* Gaussian random numbers, psi */

  for (i = 0; i < 8; i += 2) {
    r1 = sqrt( -log(genrand()) );
    r2 = genrand();
    psi[i] = r1 * cos(pi2 * r2);
    psi[i + 1] = r1 * sin(pi2 * r2);
  }

  /* solve L^T * dph = psi */

  for (i = nph - 1; i >= 0; i--) {
    for (sum = psi[i] , k = i + 1; k < nph; k++) if (k < nph) sum -= A[k][i] * dph[k];
    dph[i] = sum / p[i];
  }

  /* detailed balance */

  sum = 0; for (i = 0; i < nph; i++) sum += psi[i] * psi[i];
  wfw = exp(-sum); 
  for (i = 0; i < nph; i++) wfw *= p[i]; 

  /** backward probability **/

  for (i = 0; i < nph; i++) ph[iph[i]] += dph[i];
  int2cart(ia2,-1);

  if (!strcmp(move,"move")) 
    return 1;
  
  nph = 0;
  for (i = 0; i < 4; i++){
    ix[nph] = x[iCa[ia1 + i]];
    iy[nph] = y[iCa[ia1 + i]];
    iz[nph] = z[iCa[ia1 + i]];
    bx[nph] = x[iCa[ia1 + i]] - x[iN[ia1 + i]];
    by[nph] = y[iCa[ia1 + i]] - y[iN[ia1 + i]];
    bz[nph] = z[iCa[ia1 + i]] - z[iN[ia1 + i]];
    renorm = sqrt(bx[nph]*bx[nph] + by[nph]*by[nph] + bz[nph]*bz[nph]);
    bx[nph] /= renorm; by[nph] /= renorm; bz[nph] /= renorm;
    nph++;

    ix[nph] = x[iC[ia1 + i]];
    iy[nph] = y[iC[ia1 + i]];
    iz[nph] = z[iC[ia1 + i]];
    bx[nph] = x[iC[ia1 + i]] - x[iCa[ia1 + i]];
    by[nph] = y[iC[ia1 + i]] - y[iCa[ia1 + i]];
    bz[nph] = z[iC[ia1 + i]] - z[iCa[ia1 + i]];
    renorm = sqrt(bx[nph]*bx[nph] + by[nph]*by[nph] + bz[nph]*bz[nph]);
    bx[nph] /= renorm; by[nph] /= renorm; bz[nph] /= renorm;
    nph++;
  } 

  i = iCa[ia2];     rx[0] = x[i]; ry[0] = y[i]; rz[0] = z[i];
  i = iC[ia2];      rx[1] = x[i]; ry[1] = y[i]; rz[1] = z[i];
  i = iC[ia2] + 1;  rx[2] = x[i]; ry[2] = y[i]; rz[2] = z[i];
  for (i = 0; i < 3; i++) {
    for (j = 0; j < nph; j++) {
      dx[i][j] = by[j] * (rz[i]-iz[j]) - bz[j] * (ry[i]-iy[j]) ;
      dy[i][j] = bz[j] * (rx[i]-ix[j]) - bx[j] * (rz[i]-iz[j]) ;
      dz[i][j] = bx[j] * (ry[i]-iy[j]) - by[j] * (rx[i]-ix[j]) ;
    }
  }

  for (i=0;i<nph;i++) {
    for (j=i;j<nph;j++) {
      A[i][j]=0;
      for (k=0;k<3;k++) {
	A[i][j]+=dx[k][i]*dx[k][j]+dy[k][i]*dy[k][j]+dz[k][i]*dz[k][j];
      }
      A[i][j]*=bbgs;
      if (i==j) A[i][j]+=1;
      A[i][j]*=abgs/2;
    }
  }

  for (i=0;i<nph;i++) {
    for (j=i;j<nph;j++) {
      for (sum=A[i][j],k=i-1;k>=0;k--) sum-=A[i][k]*A[j][k];
      if (i==j)
	p[i]=sqrt(sum);
      else
	A[j][i]=sum/p[i];
    }
  }

  for (i=0;i<nph;i++) {
    psi[i] = p[i] * dph[i]; for (j = i + 1; j < nph; j++) if (j < nph) psi[i] += A[j][i] * dph[j];  
  }

  sum = 0; for (i = 0; i < nph; i++) sum += psi[i] * psi[i];
  wbw = exp(-sum); 
  for (i = 0; i < nph; i++) wbw *= p[i]; 			  

  /** accept/reject **/

  de_acc = - log( (wfw / wbw) * genrand() ) / beta[ind];
  
  cellcalc("bgs","after",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,e_max);
  de += de_exvol + de_hbonds + de_cc + de_ca;

  if (de < e_max) {
    localsa("bgs","after",ia1,ia2,&de_local);
    bias("bgs","after",ia1,ia2,&de_bias);
    hp("bgs","after",Nb[ic],ia2,&de_hp);
    de += de_local + de_bias + de_hp;
  }

  if (de < de_acc) {
    for (i = i1; i <= i2; i++) {
      xo[i] = x[i]; 
      yo[i] = y[i]; 
      zo[i] = z[i];
    }
    
    E = ( (Ecc   += de_cc) +
	  (Ecp   += de_ca) +
	  (Eev   += de_exvol) + 
	  (Ehb   += de_hbonds) +
	  (Eloc  += de_local) + 
	  (Ebias += de_bias) +
	  (Ehp   += de_hp) );

    cellcalc("bgs","acc",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);

    return 1;
  }
  else {
    for (i = 0; i < nph; i++) ph[iph[i]] -= dph[i];
    for (i = i1; i <= i2; i++) {
      x[i] = xo[i]; 
      y[i] = yo[i]; 
      z[i] = zo[i];
    }
    cellcalc("bgs","rej",0,i1,i2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    return 0;
  }
}
/****************************************************************************/
void pivot(int iflag,long *nup,double *accup,long *ndw,double *accdw) {
  int i,ia,acc,ia0,ic;
  static long int natt[N],nacc[N];
  FILE *fp;
  
  if (iflag < 0) {
    printf("Init pivot...\n");
    pivup(-1,0);
    pivdw(-1,0);
    for (i = 0; i < N; i++) natt[i] = nacc[i] = 0;
    fp = fopen("_piv_acc","w");
    fclose(fp);
    printf("\n ...DONE\n");
    return ;
  }

  if (iflag == 0) {

    ia = N * genrand();
    ic = a2c[ia];
    ia0 = Nb[ic] + (Ne[ic] - Nb[ic] + 1) / 2;

    if (ia > ia0 + 1) {
      (*nup)++; (*accup) += (acc = pivup(0,ia));
    }

    if (ia < ia0 - 1) {
      (*ndw)++; (*accdw) += (acc = pivdw(0,ia));
    }

    if (ia >= ia0 - 1 && ia <= ia0 + 1) {
      if (genrand() < 0.5) {
	(*nup)++; (*accup) += (acc = pivup(0,ia));
      } else {
	(*ndw)++; (*accdw) += (acc = pivdw(0,ia));
      }
    }
    
    
    natt[ia] += 1;
    nacc[ia] += acc;

    return ;
  }

  if (iflag > 0) {
    fp = fopen("_piv_acc","w");

    for (i = 0; i < N; i++) {
      fprintf(fp,"%i %2li %2li %2.8lf\n",i,natt[i],nacc[i],(double) nacc[i] / (double) natt[i]);
    }
    
    fclose(fp);
  }

  return;
}
/****************************************************************************/
int pivup(int iflag,int ia)
{
  int i,j,i1,i2,ia1,ia2,iup,ic;
  double phnew;
  double de = 0,de_cc = 0,de_ca = 0,de_acc;
  double de_exvol = 0,de_hbonds = 0;
  double de_local = 0,de_bias = 0,de_hp = 0;
  double ex,ey,ez;
  double dx,dy,dz; 
  double Axx,Axy,Axz,Ayy,Ayz,Azz;
  double Bxx,Bxy,Bxz,Byy,Byz,Bzz;
  double Cxx,Cxy,Cxz,Cyx,Cyy,Cyz,Czx,Czy,Czz;
  double dph,cdph,sdph,renorm;
  
  if (NCH == 0) return 0;
  
  if (iflag < 0) {
    return -1;
  }

  iup = (genrand() < 0.5 ? 3 * ia : 3 * ia + 1);
  ic = a2c[ia];
  ia2 = Ne[ic];

  //  printf("pivup ia %i\n",ia);
  
  if (iup%3 == 0) {
    i1 = iN[ia]; i2 = iCa[ia]; ia1 = ia;
  }
  else if (iup%3 == 1) {
    i1 = iCa[ia]; i2 = iC[ia]; ia1 = ia + 1;
  }
  else {
    printf("%i -- error in pivot\n",iup); exit(-1);
  }

  /* energy calc */

  cellcalc("piv","before",0,i2 + 1,iEnd[ic],&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
  localsa("piv","before",ia,ia,&de_local);
  bias("piv","before",ia,ia,&de_bias);
  hp("piv","before",ia1,ia2,&de_hp);

  /* move */
  
  dph = PIVSTP * (genrand() - 0.5);
  phnew = ph[iup] + dph;

  cdph = cos(dph); sdph = sin(dph);
  ex = x[i2] - x[i1];
  ey = y[i2] - y[i1];
  ez = z[i2] - z[i1];
  renorm=sqrt(ex*ex+ey*ey+ez*ez);
  ex/=renorm; ey/=renorm; ez/=renorm;

  Axx=ex*ex; Axy=ex*ey; Axz=ex*ez; Ayy=ey*ey; Ayz=ey*ez; Azz=ez*ez;
  Bxx=Byy=Bzz=cdph; Bxy=-ez*sdph; Byz=-ex*sdph; Bxz=ey*sdph; 
  Cxx=Bxx+(1-Bxx)*Axx-Bxy*Axy-Bxz*Axz;
  Cxy=Bxy+(1-Bxx)*Axy-Bxy*Ayy-Bxz*Ayz;
  Cxz=Bxz+(1-Bxx)*Axz-Bxy*Ayz-Bxz*Azz;
  Cyx=-Bxy+Bxy*Axx+(1-Byy)*Axy-Byz*Axz;
  Cyy=Byy+Bxy*Axy+(1-Byy)*Ayy-Byz*Ayz;
  Cyz=Byz+Bxy*Axz+(1-Byy)*Ayz-Byz*Azz;
  Czx=-Bxz+Bxz*Axx+Byz*Axy+(1-Bzz)*Axz;
  Czy=-Byz+Bxz*Axy+Byz*Ayy+(1-Bzz)*Ayz;
  Czz=Bzz+Bxz*Axz+Byz*Ayz+(1-Bzz)*Azz;

  for (j = i2 + 1; j <= iEnd[ic]; j++) {
    dx=x[j]-x[i2]; dy=y[j]-y[i2]; dz=z[j]-z[i2];
    x[j]=x[i2]+Cxx*dx+Cxy*dy+Cxz*dz; 
    y[j]=y[i2]+Cyx*dx+Cyy*dy+Cyz*dz; 
    z[j]=z[i2]+Czx*dx+Czy*dy+Czz*dz; 
  } 

  de_acc = -log(genrand()) / beta[ind];
  
  cellcalc("piv","after",0,i2 + 1,iEnd[ic],&de_exvol,&de_hbonds,&de_cc,&de_ca,e_max);  
  de += de_exvol + de_hbonds + de_cc + de_ca;

  if (de < e_max) {
    localsa("piv","after",ia,ia,&de_local);
    bias("piv","after",ia,ia,&de_bias);
    hp("piv","after",ia1,ia2,&de_hp);
    de += de_local + de_bias + de_hp;
  }
  
  if (de < de_acc) {  
    ph[iup] = phnew;
    for (i = i2 + 1; i <= iEnd[ic]; i++) {
      xo[i] = x[i]; 
      yo[i] = y[i]; 
      zo[i] = z[i];
    }
    E = ( (Ecc   += de_cc) +
	  (Ecp   += de_ca) +
	  (Eev   += de_exvol) + 
	  (Ehb   += de_hbonds) +
	  (Eloc  += de_local) + 
	  (Ebias += de_bias) +
	  (Ehp   += de_hp) );

    cellcalc("piv","acc",0,i2 + 1,iEnd[ic],&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    return 1;
  }
  else {
    for (i = i2 + 1; i <= iEnd[ic]; i++) {
      x[i] = xo[i]; 
      y[i] = yo[i]; 
      z[i] = zo[i];
    }
    cellcalc("piv","rej",0,i2 + 1,iEnd[ic],&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    return 0;
  }
}
/****************************************************************************/
int pivdw(int iflag,int ia)
{
  int i,j,i1,i2,ia1,ia2,iup,ic;
  double phnew;
  double de = 0,de_cc = 0,de_ca = 0,de_acc;
  double de_exvol = 0,de_hbonds = 0;
  double de_local = 0,de_bias = 0,de_hp = 0;
  double ex,ey,ez;
  double dx,dy,dz; 
  double Axx,Axy,Axz,Ayy,Ayz,Azz;
  double Bxx,Bxy,Bxz,Byy,Byz,Bzz;
  double Cxx,Cxy,Cxz,Cyx,Cyy,Cyz,Czx,Czy,Czz;
  double dph,cdph,sdph,renorm;
  
  if (NCH == 0) return 0;

  if (iflag<0) {
    return -1;
  }

  iup = (genrand() < 0.5 ? 3 * ia : 3 * ia + 1);
  ic = a2c[ia];
  ia1 = Nb[ic];

  //  printf("pivdw ia %i \n",ia);

  if (iup%3 == 0) {
    i1 = iN[ia]; i2 = iCa[ia]; ia2 = ia - 1;
  }
  else if (iup%3 == 1) {
    i1 = iCa[ia]; i2 = iC[ia]; ia2 = ia;
  }
  else {
    printf("%i -- error in pivot\n",iup); exit(-1);
  }

  /* energy calc */
  
  cellcalc("piv","before",0,iBeg[ic],i1 + 2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
  localsa("piv","before",ia,ia,&de_local);
  bias("piv","before",ia,ia,&de_bias);
  hp("piv","before",ia1,ia2,&de_hp);

  /* move */
  
  dph = PIVSTP * (genrand() - 0.5); 
  phnew = ph[iup] + dph;

  cdph = cos(dph); sdph = sin(dph);
  ex = x[i1] - x[i2]; 
  ey = y[i1] - y[i2]; 
  ez = z[i1] - z[i2];
  renorm=sqrt(ex*ex+ey*ey+ez*ez);
  ex/=renorm; ey/=renorm; ez/=renorm;

  Axx=ex*ex; Axy=ex*ey; Axz=ex*ez; Ayy=ey*ey; Ayz=ey*ez; Azz=ez*ez;
  Bxx=Byy=Bzz=cdph; Bxy=-ez*sdph; Byz=-ex*sdph; Bxz=ey*sdph; 
  Cxx=Bxx+(1-Bxx)*Axx-Bxy*Axy-Bxz*Axz;
  Cxy=Bxy+(1-Bxx)*Axy-Bxy*Ayy-Bxz*Ayz;
  Cxz=Bxz+(1-Bxx)*Axz-Bxy*Ayz-Bxz*Azz;
  Cyx=-Bxy+Bxy*Axx+(1-Byy)*Axy-Byz*Axz;
  Cyy=Byy+Bxy*Axy+(1-Byy)*Ayy-Byz*Ayz;
  Cyz=Byz+Bxy*Axz+(1-Byy)*Ayz-Byz*Azz;
  Czx=-Bxz+Bxz*Axx+Byz*Axy+(1-Bzz)*Axz;
  Czy=-Byz+Bxz*Axy+Byz*Ayy+(1-Bzz)*Ayz;
  Czz=Bzz+Bxz*Axz+Byz*Ayz+(1-Bzz)*Azz;

  for (j = i1 + 2; j >= iBeg[ic]; j--) {
    if (j == i1 || j == i2) continue;
    dx=x[j]-x[i1]; dy=y[j]-y[i1]; dz=z[j]-z[i1];
    x[j]=x[i1]+Cxx*dx+Cxy*dy+Cxz*dz; 
    y[j]=y[i1]+Cyx*dx+Cyy*dy+Cyz*dz; 
    z[j]=z[i1]+Czx*dx+Czy*dy+Czz*dz; 
  }

  de_acc = -log(genrand()) / beta[ind];
  
  cellcalc("piv","after",0,iBeg[ic],i1 + 2,&de_exvol,&de_hbonds,&de_cc,&de_ca,e_max);
  de += de_exvol + de_hbonds + de_cc + de_ca;

  if (de < e_max) {
    localsa("piv","after",ia,ia,&de_local);
    bias("piv","after",ia,ia,&de_bias);
    hp("piv","after",ia1,ia2,&de_hp);
    de += de_local + de_bias + de_hp;
  }
  
  if (de < de_acc) {
    ph[iup] = phnew;
    for (i = i1 + 2; i >= iBeg[ic]; i--) {
      if (i == i1 || i == i2) continue;
      xo[i] = x[i]; 
      yo[i] = y[i]; 
      zo[i] = z[i];
    }
    
    E = ( (Ecc   += de_cc) +
	  (Ecp   += de_ca) +
	  (Eev   += de_exvol) + 
	  (Ehb   += de_hbonds) +
	  (Eloc  += de_local) + 
	  (Ebias += de_bias) +
	  (Ehp   += de_hp) );

    cellcalc("piv","acc",0,iBeg[ic],i1 + 2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    return 1;
  }
  else {
    for (i = i1 + 2; i >= iBeg[ic]; i--) {
      if (i == i1 || i == i2) continue;
      x[i] = xo[i]; 
      y[i] = yo[i]; 
      z[i] = zo[i];
    }
    cellcalc("piv","rej",0,iBeg[ic],i1 + 2,&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    return 0;
  }
}
/****************************************************************************/
int trans(int iflag) 
{
  double dx,dy,dz;
  double de = 0,de_acc;
  double de_cc = 0,de_ca = 0;
  double de_exvol = 0,de_hbonds = 0;
  double de_hp = 0;
  int i,ic;
  static long int natt[NCH],nacc[NCH];
  FILE *fp;

  if (NCH == 0) return 0;

  if (iflag < 0) {
    printf("init trans\n");
    fp = fopen("_trans_acc","w");
    fclose(fp);
    for (i = 0; i < NCH; i++) natt[i] = nacc[i] = 0;
    printf("  .. DONE\n");    
    return 0;
  }
  
  if (iflag == 0) {
    ic = NCH * genrand();

    natt[ic]++;

    /* energy calc */
    
    cellcalc("piv","before",0,iBeg[ic],iEnd[ic],&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
    hp("piv","before",Nb[ic],Ne[ic],&de_hp);
    
    /* move */
    
    dx = TRLSTP * (genrand() - .5);
    dy = TRLSTP * (genrand() - .5);
    dz = TRLSTP * (genrand() - .5);
    
    //  printf("trans: ich %i dx %f dy %f dz %f\n",ich,dx,dy,dz);
    
    movech(dx,dy,dz,ic);
    
    de_acc = -log( genrand() ) / beta[ind];
    
    cellcalc("piv","after",0,iBeg[ic],iEnd[ic],&de_exvol,&de_hbonds,&de_cc,&de_ca,e_max);
    de += de_exvol + de_hbonds + de_cc + de_ca;
    
    hp("piv","after",Nb[ic],Ne[ic],&de_hp);
    de += de_hp;
    
    if (de < de_acc) {
      nacc[ic]++;
      for (i=iBeg[ic];i<=iEnd[ic];i++) {
	xo[i]=x[i]; 
	yo[i]=y[i]; 
	zo[i]=z[i];
      }
      
      E = ( (Eloc + Ebias) +
	    (Ecc   += de_cc) +
	    (Ecp   += de_ca) +
	    (Eev   += de_exvol) + 
	    (Ehb   += de_hbonds) +
	    (Ehp   += de_hp) );
      
      cellcalc("piv","acc",0,iBeg[ic],iEnd[ic],&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
      return 1;
    } else {
      for (i=iBeg[ic];i<=iEnd[ic];i++) {
	x[i]=xo[i]; 
	y[i]=yo[i]; 
	z[i]=zo[i];
      }
      cellcalc("piv","rej",0,iBeg[ic],iEnd[ic],&de_exvol,&de_hbonds,&de_cc,&de_ca,0);
      return 0;
    }
    return 0;
  }

  if (iflag > 0) {
    fp = fopen("_trans_acc","w");
    for (i = 0; i < NCH; i++) 
      fprintf(fp,"%i %2li %2li %2.8lf\n",i,natt[i],nacc[i],(double) nacc[i] / (double) natt[i]);
    fclose(fp);
    return 0;
  }

  return 0;  
}
/****************************************************************************/
void movech(double dx,double dy,double dz,int ic) 
{
  int i;
  for (i=iBeg[ic];i<=iEnd[ic];i++) {
    x[i]+=dx;
    y[i]+=dy;
    z[i]+=dz;
  }
}
/****************************************************************************/
int flip(double e)
{
  /* sampling joint distribution p(E,k) ~ exp(- beta_k * E - g_k) */
  /* marginal distribution p(k) flat if g_k = - beta F_k */

  int nind = (genrand() > 0.5) ? ind + 1 : ind - 1; 

  if (nind == NTMP || nind == -1)
    return 0;

  if ( genrand() < exp( -( (beta[nind] - beta[ind]) * e + g[nind] - g[ind] ) ) ) { 
    ind = nind; 
    return 1;
  } else {
    return 0;
  }
}
/****************************************************************************/
int get_cell(int a,int *pnt) {
  if (a >= 0) while ((a = pnt[a]) >= 0);
  return -a-2;  
}
/****************************************************************************/
void get_cell_xyz(short *ix,short *iy,short *iz,int ic,int ns) {
  /* ic = ix + ns * (iy + ns * iz); */
  int h2=ns*ns;
  (*iz) = ic / h2; ic -= (*iz) * h2;
  (*iy) = ic / ns; ic -= (*iy) * ns;
  (*ix) = ic;
  
}
/****************************************************************************/
void remove_from_list(int ic,int *pnt,short *cell,int aout) {
  int a = cell[ic],a0; 

  if (a == aout) {
    cell[ic] = max(pnt[a],-1);
    pnt[a] = -1;
    return;
  } 
  
  while ((a = pnt[(a0 = a)]) >= 0 && a != aout) ;
  if (a < 0) {printf("Error remove_from_list... exiting\n"); exit(-1);}
  pnt[a0] = pnt[a];
  pnt[a] = -1;
  return ;
}
/****************************************************************************/
void insert_into_list(int ic,int *pnt,short *cell,int ain) {
  int a = cell[ic],a0;

  if (a < 0 || a > ain) {
    cell[ic] = ain;
    pnt[ain] = (a < 0 ? -2-ic : a);
    return ;
  } 

  while ((a = pnt[(a0 = a)]) >= 0 && a < ain);
  if (a == ain) {printf("Error insert_into_list... exiting\n"); exit(-1);}
  pnt[ain] = a;
  pnt[a0] = ain;
}
/****************************************************************************/
double dist2_cell_crowd(short ix, short iy,short iz,int i,double cutg) {
  double ixc = ix;
  double iyc = iy;
  double izc = iz;
  double dx,dy,dz;

  dx = xcr[i] - (ixc + .5) * cutg; bc(&dx);
  dy = ycr[i] - (iyc + .5) * cutg; bc(&dy);
  dz = zcr[i] - (izc + .5) * cutg; bc(&dz); 
  
  dx = max(0,fabs(dx) - .5 * cutg);
  dy = max(0,fabs(dy) - .5 * cutg);
  dz = max(0,fabs(dz) - .5 * cutg);  

  return dx * dx + dy * dy + dz * dz;
}
/****************************************************************************/
void build_list_crowd(short ix,short iy,short iz,int *nlc,short *listc,
		      double cutg) {
  int i;

  for (i = 0; i < NCR; i++) {
    if (dist2_cell_crowd(ix,iy,iz,i,cutg) < cut_crowd_atom2) 
      listc[(*nlc)++] = i;
  }
}
/****************************************************************************/
void build_list_inc(int a,int *pnt,int *nl,short *list,int sel,int a1,int a2) {

  if (a < 0)
    return;

  if (sel < 0) { /* exclude a1..a2 */

    if (a < a1 || a > a2)
      list[(*nl)++] = a;

    while ((a = pnt[a]) >= 0)
      if (a < a1 || a > a2)
	list[(*nl)++] = a;

    return ;
  }

  if (sel > 0) { /* include a1..a2 */

    if (a >= a1 && a <= a2)
      list[(*nl)++] = a;

    while ((a = pnt[a]) >= 0)
      if (a >= a1 && a <= a2)
	list[(*nl)++] = a;

    return ;
  }
  
  /* include all */

  list[(*nl)++] = a;
  while ((a = pnt[a]) >= 0) list[(*nl)++] = a;
  
  return ;
}
/****************************************************************************/
void build_list_fwd(int ic,int *pnt,short *cell,
		    short ix,short iy,short iz,
		    int *nl,short *list,int ns,
		    int sel,int a1,int a2) { 
  int cn;
  int h2=ns*ns;
  int h3=h2*ns;

  cn=ic+1; 
  if (ix+1==ns) cn-=ns;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+ns; 
  if (iy+1==ns) cn-=h2;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+h2; 
  if (iz+1==ns) cn-=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+1+ns; 
  if (ix+1==ns) cn-=ns; if (iy+1==ns) cn-=h2;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+1+h2; 
  if (ix+1==ns) cn-=ns; if(iz+1==ns) cn-=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+ns+h2;
  if (iy+1==ns) cn-=h2; if (iz+1==ns) cn-=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+1-ns; 
  if (ix+1==ns) cn-=ns; if (iy==0) cn+=h2; 
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+1-h2; 
  if (ix+1==ns) cn-=ns; if (iz==0) cn+=h3; 
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+ns-h2;
  if (iy+1==ns) cn-=h2; if (iz==0) cn+=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+1+ns+h2;
  if (ix+1==ns) cn-=ns; if (iy+1==ns) cn-=h2; if (iz+1==ns) cn-=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+1+ns-h2;
  if (ix+1==ns) cn-=ns; if (iy+1==ns) cn-=h2; if (iz==0) cn+=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+1-ns+h2;
  if (ix+1==ns) cn-=ns; if (iy==0) cn+=h2; if (iz+1==ns) cn-=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-1+ns+h2;
  if (ix==0) cn+=ns; if (iy+1==ns) cn-=h2; if (iz+1==ns) cn-=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  return ;
}
/****************************************************************************/
void build_list_bwd(int ic,int *pnt,short *cell,
		    short ix,short iy,short iz,
		    int *nl,short *list,int ns,
		    int sel,int a1,int a2) {
  int cn;
  int h2=ns*ns;
  int h3=h2*ns;

  cn=ic-1; 
  if (ix==0) cn+=ns;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-ns; 
  if (iy==0) cn+=h2;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-h2; 
  if (iz==0) cn+=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-1-ns; 
  if (ix==0) cn+=ns; if (iy==0) cn+=h2;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-1-h2; 
  if (ix==0) cn+=ns; if(iz==0) cn+=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-ns-h2;
  if (iy==0) cn+=h2; if (iz==0) cn+=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-1+ns; 
  if (ix==0) cn+=ns; if (iy+1==ns) cn-=h2; 
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-1+h2; 
  if (ix==0) cn+=ns; if (iz+1==ns) cn-=h3; 
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-ns+h2;
  if (iy==0) cn+=h2; if (iz+1==ns) cn-=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-1-ns-h2;
  if (ix==0) cn+=ns; if (iy==0) cn+=h2; if (iz==0) cn+=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-1-ns+h2;
  if (ix==0) cn+=ns; if (iy==0) cn+=h2; if (iz+1==ns) cn-=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic-1+ns-h2;
  if (ix==0) cn+=ns; if (iy+1==ns) cn-=h2; if (iz==0) cn+=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  cn=ic+1-ns-h2;
  if (ix+1==ns) cn-=ns; if (iy==0) cn+=h2; if (iz==0) cn+=h3;
  build_list_inc(cell[cn],pnt,nl,list,sel,a1,a2);

  return ;
}
/****************************************************************************/
int crowd_clash(){
  int i,j,n = 0;
  double dx,dy,dz,r2;
  double r2min_crcr = 2 * (rcrowd - sigcr) * 2 * (rcrowd - sigcr);
  double r2min_crat = (rcrowd - sigcr) * (rcrowd - sigcr);
  
  for (i=0;i<NCR;i++)  {
    for (j=0;j<i;j++)  {
      dx = xcr[i] - xcr[j]; bc(&dx);
      dy = ycr[i] - ycr[j]; bc(&dy);
      dz = zcr[i] - zcr[j]; bc(&dz);
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < r2min_crcr) {
	printf("crowd-crowd clash: %i %i %lf %lf\n",i,j,sqrt(r2min_crcr),sqrt(r2));
	n++;
      }
    }
  }

  for (i=0;i<NCR;i++)  {
    for (j=0;j<NTO;j++)  {
      dx = xcr[i] - x[j]; bc(&dx);
      dy = ycr[i] - y[j]; bc(&dy);
      dz = zcr[i] - z[j]; bc(&dz);
      r2 = dx*dx + dy*dy + dz*dz;
      if (r2 < r2min_crat) {
	printf("crowd-atom clash: %i %i %lf %lf\n",i,j,sqrt(r2min_crat),sqrt(r2));
	n++;
      }
    }
  }

  return n;
}
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/


# if 0


/****************************************************************************/
int bgsdw_wrong(char *move, int ia2)
{
  int i,j,k,ia1,ic,i2,iph[8],nph=0;
  double de=0;
  //  double e1=Ebias,e2=Ehb,e3=Ehp,e4=Eloc;
  //  double de = 0,de_crowd = 0,de_exvol = 0;
  double ix[8],iy[8],iz[8],bx[8],by[8],bz[8],ab[8];
  double rx[3],ry[3],rz[3];
  double dx[3][8],dy[3][8],dz[3][8];
  double A[8][8],p[8],psi[8],dph[8],sum,r1,r2;
  double wfw,wbw;

  ia1=ia2-3;
  ic=a2c[ia1];

  //  i1=iN[ia1]; 
  i2=iC[ia2]-1; 
  // printf("dw ia1 %i ia2 %i\n",ia1,ia2);

  /* energy calc */

  //  if (!strcmp(move,"mc")) {
  //    de_crowd -= crowders("chain","before",0,iBeg[ic],i2,0);
  //    de_exvol -= exvol("bgs","before",iBeg[ic],i2,0);
  //  }
  
  /** forward move **/

  for (i=0;i<4;i++){
    ix[nph]=x[iCa[ia2-i]];
    iy[nph]=y[iCa[ia2-i]];
    iz[nph]=z[iCa[ia2-i]];
    bx[nph]=x[iCa[ia2-i]]-x[iC[ia2-i]];
    by[nph]=y[iCa[ia2-i]]-y[iC[ia2-i]];
    bz[nph]=z[iCa[ia2-i]]-z[iC[ia2-i]];
    ab[nph]=b[2];
    iph[nph]=3*(ia2-i)+1;
    nph++;
    if (seq[ia2-i]!='P') {
      ix[nph]=x[iN[ia2-i]];
      iy[nph]=y[iN[ia2-i]];
      iz[nph]=z[iN[ia2-i]];
      bx[nph]=x[iN[ia2-i]]-x[iCa[ia2-i]];
      by[nph]=y[iN[ia2-i]]-y[iCa[ia2-i]];
      bz[nph]=z[iN[ia2-i]]-z[iCa[ia2-i]];
      ab[nph]=b[1];
      iph[nph]=3*(ia2-i);
      nph++;
    }
  } 
  i=iCa[ia1];    rx[0]=x[i]; ry[0]=y[i]; rz[0]=z[i];
  i=iN[ia1];     rx[1]=x[i]; ry[1]=y[i]; rz[1]=z[i];
  i=iN[ia1]+1;   rx[2]=x[i]; ry[2]=y[i]; rz[2]=z[i]; // H
  for (i=0;i<3;i++) {
    for (j=0;j<nph;j++) {
      dx[i][j]=(by[j]*(rz[i]-iz[j])-bz[j]*(ry[i]-iy[j]))/ab[j];
      dy[i][j]=(bz[j]*(rx[i]-ix[j])-bx[j]*(rz[i]-iz[j]))/ab[j];
      dz[i][j]=(bx[j]*(ry[i]-iy[j])-by[j]*(rx[i]-ix[j]))/ab[j];
    }
  }

  /* Matrix A */

  for (i=0;i<nph;i++) {
    for (j=i;j<nph;j++) {
      A[i][j]=0;
      for (k=0;k<3;k++) {
	A[i][j]+=dx[k][i]*dx[k][j]+dy[k][i]*dy[k][j]+dz[k][i]*dz[k][j];
      }
      A[i][j]*=bbgs;
      if (i==j) A[i][j]+=1;
      A[i][j]*=abgs/2;
    }
  }

  /* Cholesky decomp of A=LL^T (put in lower triangle A and p[]) */

  for (i=0;i<nph;i++) {
    for (j=i;j<nph;j++) {
      for (sum=A[i][j],k=i-1;k>=0;k--) sum-=A[i][k]*A[j][k];
      if (i==j)
	p[i]=sqrt(sum);
      else
	A[j][i]=sum/p[i];
    }
  }

  /* Gaussian random numbers psi */

  for (i=0;i<8;i+=2) {
    r1=sqrt(-log( genrand() ));
    r2=genrand();
    psi[i]=r1*cos(pi2*r2);
    psi[i+1]=r1*sin(pi2*r2);
  }

  /* solve L^T * dph = psi */

  for (i=nph-1;i>=0;i--) {
    for (sum=psi[i],k=i+1;k<nph;k++) if (k<nph) sum-=A[k][i]*dph[k];
    dph[i]=sum/p[i];
  }

  /* detailed balance */

  sum=0; for (i=0;i<nph;i++) sum+=psi[i]*psi[i];
  wfw=exp(-sum); 
  for (i=0;i<nph;i++) wfw*=p[i]; 
  
  /** backward probability **/

  for (i=0;i<nph;i++) ph[iph[i]]+=dph[i];
  int2cart(ia2,-1);

  if (!strcmp(move,"move"))
    return 1;
  
  nph=0;
  for (i=0;i<4;i++){
    ix[nph]=x[iCa[ia2-i]];
    iy[nph]=y[iCa[ia2-i]];
    iz[nph]=z[iCa[ia2-i]];
    bx[nph]=x[iCa[ia2-i]]-x[iC[ia2-i]];
    by[nph]=y[iCa[ia2-i]]-y[iC[ia2-i]];
    bz[nph]=z[iCa[ia2-i]]-z[iC[ia2-i]];
    ab[nph]=b[2];
    iph[nph]=3*(ia2-i)+1;
    nph++;
    if (seq[ia2-i]!='P') {
      ix[nph]=x[iN[ia2-i]];
      iy[nph]=y[iN[ia2-i]];
      iz[nph]=z[iN[ia2-i]];
      bx[nph]=x[iN[ia2-i]]-x[iCa[ia2-i]];
      by[nph]=y[iN[ia2-i]]-y[iCa[ia2-i]];
      bz[nph]=z[iN[ia2-i]]-z[iCa[ia2-i]];
      ab[nph]=b[1];
      iph[nph]=3*(ia2-i);
      nph++;
    }
  } 
  i=iCa[ia1];    rx[0]=x[i]; ry[0]=y[i]; rz[0]=z[i];
  i=iN[ia1];     rx[1]=x[i]; ry[1]=y[i]; rz[1]=z[i];
  i=iN[ia1]+1;   rx[2]=x[i]; ry[2]=y[i]; rz[2]=z[i]; 
  for (i=0;i<3;i++) {
    for (j=0;j<nph;j++) {
      dx[i][j] = ( by[j] * (rz[i]-iz[j]) - bz[j] * (ry[i]-iy[j]) ) / ab[j];
      dy[i][j] = ( bz[j] * (rx[i]-ix[j]) - bx[j] * (rz[i]-iz[j]) ) / ab[j];
      dz[i][j] = ( bx[j] * (ry[i]-iy[j]) - by[j] * (rx[i]-ix[j]) ) / ab[j];
    }
  }
  for (i=0;i<nph;i++) {
    for (j=i;j<nph;j++) {
      A[i][j]=0;
      for (k=0;k<3;k++) {
	A[i][j]+=dx[k][i]*dx[k][j]+dy[k][i]*dy[k][j]+dz[k][i]*dz[k][j];
      }
      A[i][j]*=bbgs;
      if (i==j) A[i][j]+=1;
      A[i][j]*=abgs/2;
    }
  }
  for (i=0;i<nph;i++) {
    for (j=i;j<nph;j++) {
      for (sum=A[i][j],k=i-1;k>=0;k--) sum-=A[i][k]*A[j][k];
      if (i==j)
	p[i]=sqrt(sum);
      else
	A[j][i]=sum/p[i];
    }
  }
  for (i=0;i<nph;i++) {
    psi[i]=p[i]*dph[i]; for (j=i+1;j<nph;j++) if (j<nph) psi[i]+=A[j][i]*dph[j];  
  }
  sum=0; for (i=0;i<nph;i++) sum+=psi[i]*psi[i];
  wbw=exp(-sum); 
  for (i=0;i<nph;i++) wbw*=p[i]; 			  
    
  if (de < - log( (wfw / wbw) * genrand() ) / beta[ind] ) {
    for (i=i2;i>=iBeg[ic];i--) {
      xo[i]=x[i]; 
      yo[i]=y[i]; 
      zo[i]=z[i];
    }
    return 1;
  }
  else {
    for (i=0;i<nph;i++) ph[iph[i]]-=dph[i];
    for (i=i2;i>=iBeg[ic];i--) {
      x[i]=xo[i]; 
      y[i]=yo[i]; 
      z[i]=zo[i];
    }

    return 0;
  }
}

#endif
