# include <stdio.h>
# include <math.h>
# include <string.h>
# include <stdlib.h>
# include "sys.h"
# include "defs.h"
# include "global.h"

/****************************************************************************/
# define BB "N  ","CA ","C  ","O  "
const char atomname[20][15][4]={
{BB},                                                                /* gly */
{BB,"CB "},                                                          /* ala */
{BB,"CB ","CG1","CG2"},                                              /* val */
{BB,"CB ","CG ","CD1","CD2"},                                        /* leu */
{BB,"CB ","CG1","CG2","CD1"},                                        /* ile */
{BB,"CB ","OG "},                                                    /* ser */
{BB,"CB ","OG1","CG2"},                                              /* thr */
{BB,"CB ","SG "},                                                    /* cys */
{BB,"CB ","CG ","SD ","CE "},                                        /* met */
{BB,"CB ","CG ","CD "},                                              /* pro */
{BB,"CB ","CG ","OD1","OD2"},                                        /* asp */
{BB,"CB ","CG ","OD1","ND2"},                                        /* asn */
{BB,"CB ","CG ","CD ","OE1","OE2"},                                  /* glu */
{BB,"CB ","CG ","CD ","OE1","NE2"},                                  /* gln */
{BB,"CB ","CG ","CD ","CE ","NZ "},                                  /* lys */
{BB,"CB ","CG ","CD ","NE ","CZ ","NH1","NH2"},                      /* arg */
{BB,"CB ","CG ","ND1","CD2","CE1","NE2"},                            /* his */
{BB,"CB ","CG ","CD1","CD2","CE1","CE2","CZ "},                      /* phe */
{BB,"CB ","CG ","CD1","CD2","CE1","CE2","CZ ","OH "},                /* tyr */
{BB,"CB ","CG ","CD1","CD2","NE1","CE2","CE3","CZ2","CZ3","CH2"}     /* trp */
};
const char aminoname[20][4]={"GLY","ALA","VAL","LEU","ILE","SER","THR",
		             "CYS","MET","PRO","ASP","ASN","GLU","GLN",
			     "LYS","ARG","HIS","PHE","TYR","TRP"};
const char aminolett[20]={'G','A','V','L','I','S','T','C','M','P',
			  'D','N','E','Q','K','R','H','F','Y','W'};
const int aminosize[20]={4,5,7,8,8,6,7,6,8,7,8,8,9,9,9,11,10,11,12,14};
const int nrotamers[20]={0,1,3,4,4,2,3,2,4,0,2,3,3,4,5,4,2,2,3,2};
#undef BB 
/****************************************************************************/
int get_atm(int iaa,char *name) 
{
  int i=0;
  while (i<20 && strcmp(atomname[iaa][i],name)!=0) i++;
  return i;
}
/****************************************************************************/
int get_aatyp1(char c) 
{
  int i=0;
  while (i<20 && aminolett[i]!=c) i++;
  return i;
}
/****************************************************************************/
int get_aatyp3(char *name) 
{
  int i=0;
  while (i<20 && strcmp(aminoname[i],name)!=0) i++;
  return i;
}
/****************************************************************************/
void get_aminoname(int ia,char *name) {
  strcpy(name,aminoname[get_aatyp1(seq[ia])]);
  return ;
}
/****************************************************************************/
void get_atomname(int iatm,char *name) 
{
  int j,k,ia,aminotyp;

  ia=i2a[iatm];
  aminotyp=get_aatyp1(seq[ia]);
  
  if (iatm==iN[ia]-1 && seq[ia]!='P') {strcpy(name,"H  "); return;}
  if (iatm==iN[ia]) {strcpy(name,"N  "); return;}
  if (iatm==iCa[ia]) {strcpy(name,"CA "); return;}
  if (iatm==iCa[ia]+1) {strcpy(name,"HA "); return;}
  if (iatm==iCa[ia]+2 && seq[ia]=='G') {strcpy(name,"HA2"); return;}
  if (iatm==iCa[ia]+2 && seq[ia]!='G') {strcpy(name,"CB "); return;}
  if (iatm==iC[ia]) {strcpy(name,"C  "); return;}
  if (iatm==iC[ia]+1) {strcpy(name,"O  "); return;}

  k=5;
  for (j=iCa[ia]+3;j<iC[ia];j++) {
    if (atyp[j]==3 || atyp[j]==6 || atyp[j]<0) continue;
    if (iatm==j) {strcpy(name,atomname[aminotyp][k]); return;}
    k++;
  }
  strcpy(name,"HX "); return ;
}

/****************************************************************************/
/**** PDB INPUT/OUTPUT ******************************************************/
/****************************************************************************/

# define KMAX 10000 

int readpdb(int flag,char *fn,double x[],double y[],double z[],int n1,int n2)
/****************************************************************************/
/* flag  = 1, Ca atoms                                                      */
/* flag  = 2, N,Ca,C atoms                                                  */
/* flag >= 3, all atoms (N,Ca,C,O + non-H sidechain atoms)                  */
/****************************************************************************/
{
  int i,k=0,iaa,naa=0,natm=0,err=0,verb=1;
  FILE *fp; 

  extern int readaa(int ,FILE *,int *,double *,double *,double *,int *) ;

  if (flag<0) {flag=abs(flag); verb=0;};

  if (verb) printf("  Reading... %s\n",fn);

  fp=fopen(fn,"r");

  natm=readaa(-1,fp,&iaa,x+k,y+k,z+k,&err);

  k=0;
  for (i=0;i<=n2;i++) {
    natm=readaa(flag,fp,&iaa,x+k,y+k,z+k,&err);
    if ((natm!=1 && flag==1) || (natm!=3 && flag==2) || 
	(natm!=aminosize[iaa] && flag>=3)) {
      printf("missing atoms %s: ",fn);
      printf("%s %i, %i read\n",aminoname[iaa],naa+1,natm);
      natm=aminosize[iaa];
      if (flag==3) {printf("exiting...\n");exit(-1);}
    }
    //    printf("%f\n",x[k]);
    if (i>=n1) {
      k+=natm;
      naa++;
    }
  }

  if (verb) printf("  # amino acids %i\n",naa);
  if (verb) printf("  # atoms %i (flag %i)\n",k,flag);

  fclose(fp);
  return k;
}
/****************************************************************************/
int readaa(int flag,FILE *fp,int *aa,double *x,double *y,double *z,int *err) 
{
  static char pdbline[100];
  char str_beg[7];     /* pdb entry type            col 1-6                 */
  char str_atm[4];     /* atom type                 col 14-16               */
  char str_aa[4];      /* amino acid type           col 18-20               */
  char str_no[6];      /* amino acid no             col 22-27               */
  char str_x[9];       /* x                         col 32-38               */
  char str_y[9];       /* y                         col 40-46               */
  char str_z[9];       /* z                         col 48-54               */

  /* misc */
  int k=0,iatm,iaa,fend=1;

  extern int fgetline(char *line,int max,FILE *fp);
  extern void substr(char *sub,char *str,int pos, int length);

  /* find next amino acid */
  
  if (flag<0) {
    strcpy(pdbline,"");
    return 0;
  }

  if (flag==1) fend=fgetline(pdbline,100,fp);
  
  do {
    substr(str_beg,pdbline,0,6);
    substr(str_atm,pdbline,13,3);
    substr(str_aa,pdbline,17,3);
    substr(str_no,pdbline,21,5);
    substr(str_x,pdbline,31,8);
    substr(str_y,pdbline,39,8);
    substr(str_z,pdbline,47,8);

    if (strcmp(str_beg,"ATOM  ")!=0) {
      fend=fgetline(pdbline,100,fp);
      continue; 
    } 

    /* find amino acid type */
    for (iaa=0; iaa<20 && strcmp(aminoname[iaa],str_aa)!=0; iaa++);
    if (iaa==20) {
      fend=fgetline(pdbline,100,fp);
      continue; 
    } 
    
    /* find atom type */
    for (iatm=0; iatm<aminosize[iaa] && 
	   strcmp(atomname[iaa][iatm],str_atm)!=0; iatm++);
    if (iatm==aminosize[iaa]) (*err)++;

    if (iatm==0 && flag>1) break;
    if (iatm==1 && flag==1) break;

    fend=fgetline(pdbline,100,fp);
  } while (fend!=0);

  if (fend==0) return -1;

  (*aa)=iaa;
  /* store information */
  if (flag==1) {
    x[0]=atof(str_x);
    y[0]=atof(str_y);
    z[0]=atof(str_z);
    return 1;
  }
  
  do {
    /* store information */
    if ((iatm<3 && flag==2) || (iatm<aminosize[(*aa)] && flag>=3)) {
      x[iatm]=atof(str_x);
      y[iatm]=atof(str_y);
      z[iatm]=atof(str_z);
      //            printf("%s %i ***%s*** %f %f %f\n",aminoname[(*aa)],iatm,
      //      	     str_atm,x[iatm],y[iatm],z[iatm]);
      k++;
    }
    
    if (k==3 && flag==2) break;
    if (k==aminosize[(*aa)] && flag>=3) break;

    fend=fgetline(pdbline,100,fp);
    substr(str_beg,pdbline,0,6);
    substr(str_atm,pdbline,13,3);
    substr(str_aa,pdbline,17,3);
    substr(str_no,pdbline,21,5);
    substr(str_x,pdbline,31,8);
    substr(str_y,pdbline,39,8);
    substr(str_z,pdbline,47,8);

    if (strcmp(str_beg,"ATOM  ")!=0) continue;

    /* find amino acid type */
    for (iaa=0; iaa<20 && strcmp(aminoname[iaa],str_aa)!=0; iaa++);

    /* if wrong type return */
    if (iaa!=(*aa)) return k;
    
    /* find atom type */
    for (iatm=0; iatm<aminosize[iaa] && 
	   strcmp(atomname[iaa][iatm],str_atm)!=0; iatm++);

    /* if unknown type continue */
    if (iatm==aminosize[iaa]) (*err)++;    
    
    /* if next amino acid found return */
    if ((iatm==1 && flag==1) || (iatm==0 && flag>1)) return k;
    
  } while (fend!=0);

  if (fend==0) return -1;
  return k;
}
/****************************************************************************/
void substr(char *sub,char *str,int pos, int n){
  strncpy(sub,str+pos,n);
  *(sub+n)='\0';
}
/****************************************************************************/
int fgetline(char *line,int max,FILE *fp){
  if (fgets(line,max,fp)==NULL && feof(fp)>0)
    return 0;
  else
    return strlen(line);
}
/****************************************************************************/
void printatom(int iatm) {
  int ia=i2a[iatm],aai;
  char atomi[20];

  get_atomname(iatm,atomi);
  aai=get_aatyp1(seq[ia]);
  printf("%3i %s %3i %s  ",ia,aminoname[aai],iatm,atomi);
}
/****************************************************************************/
void printatom2(int iatm,int jatm) {
  int ia=i2a[iatm],aai;
  int ja=i2a[jatm],aaj;
  char atomi[20],atomj[20];

  get_atomname(jatm,atomj);
  get_atomname(iatm,atomi);
  aai=get_aatyp1(seq[ia]);
  aaj=get_aatyp1(seq[ja]);
  printf("%3i %s %3i %s -- %3i %s %3i %s  ",
	 ia,aminoname[aai],ja,aminoname[aaj],iatm,atomi,jatm,atomj);
}
/****************************************************************************/
/******** RMSD CALCULATIONS *************************************************/
/****************************************************************************/
double rmsd2(double x1[],double y1[],double z1[],
	     double x2[],double y2[],double z2[],int n){
  double gyr1,gyr2;
  double xcm1,ycm1,zcm1;
  double xcm2,ycm2,zcm2;

  gyr1=center_of_mass(x1,y1,z1,&xcm1,&ycm1,&zcm1,n);
  gyr2=center_of_mass(x2,y2,z2,&xcm2,&ycm2,&zcm2,n);
  return max(0,gyr1+gyr2-correlation(x1,y1,z1,x2,y2,z2,n));
}
/****************************************************************************/
double rmsd2_nopt(double x1[],double y1[],double z1[],
		  double x2[],double y2[],double z2[],int n){
  int i;
  double r2=0;

  for (i=0;i<n;i++) r2+=((x1[i]-x2[i])*(x1[i]-x2[i])+
			 (y1[i]-y2[i])*(y1[i]-y2[i])+
			 (z1[i]-z2[i])*(z1[i]-z2[i]));
  return r2/n;
}
/****************************************************************************/
void dumppdb(char *fn,int mc,int nobs,double o[],char format[100])
{
  const char atom[7][10] =      {"N  ","H  ","CA ","HA1","CB ","C  ","O  "};
  const char atom_gly[7][10] =  {"N  ","H  ","CA ","HA1","HA2","C  ","O  "};
  const char amino[20][10] =    {"GLY","ALA","VAL","LEU","ILE","SER","THR",
				 "CYS","MET","PRO","ASP","ASN","GLU","GLN",
				 "LYS","ARG","HIS","PHE","TYR","TRP"};
  int i,j,k,ia,n=0;
  FILE *fp;

  fp=fopen(fn,"w");

  ch2box();
  cr2box();

  fprintf(fp,"REMARK MC step %i\n",mc);
  fprintf(fp,"REMARK temp ind %i\n",ind);  
  for (j=1;j<nobs;j++) fprintf(fp,"REMARK %i %f\n",j,o[j]);

  /* protein chains */
  
  for (j=0;j<NCH;j++) {
    for (i=Nb[j];i<=Ne[j];i++) {
      switch (seq[i]) {  
	case 'G' : {ia=0; break;} 
	case 'L' : {ia=3; break;}
	case 'S' : {ia=5; break;}
	default : {printf("%i -- aa unknown\n",i);exit(-1); 
	}
      }
      for (k=0;k<7;k++) {
	if (strncmp("CA",format,2) == 0 && strncmp("CA",atom[k],2) != 0) continue; 
	if (seq[i]=='G') {
	  fprintf(fp,"ATOM   %4u  %s %s  %4u    %8.3f%8.3f%8.3f  1.00\n",
		  ++n,atom_gly[k],amino[ia],i+1,x[iN[i]+k],y[iN[i]+k],
		  z[iN[i]+k]);
	} else {
	  fprintf(fp,"ATOM   %4u  %s %s  %4u    %8.3f%8.3f%8.3f  1.00\n",
		  ++n,atom[k],amino[ia],i+1,x[iN[i]+k],y[iN[i]+k],
		  z[iN[i]+k]);
	}
      }
    }
  }

  /* crowders */
  
  for (i=n=0;i<NCR;i++) {
    fprintf(fp,"HETATM %4u  %s %s  %4u    %8.3f%8.3f%8.3f  1.00\n",
	    ++n," O ","HOH",i+1,xcr[i],ycr[i],zcr[i]);
  }
  
  fclose(fp);
}
/****************************************************************************/
double center_of_mass(double x[],double y[],double z[],
		      double *xcm,double *ycm,double *zcm,int n)
{
  int j;
  double gyr2=0;

  (*xcm)=(*ycm)=(*zcm)=0;

  for (j=0;j<n;j++) {
    (*xcm)+=x[j]; 
    (*ycm)+=y[j]; 
    (*zcm)+=z[j];
  }
  (*xcm)*=1.0/n; (*ycm)*=1.0/n; (*zcm)*=1.0/n;
  
  for (j=0;j<n;j++){
    x[j]-=(*xcm); 
    y[j]-=(*ycm); 
    z[j]-=(*zcm);
    gyr2+=x[j]*x[j]+y[j]*y[j]+z[j]*z[j];
  }
  return gyr2/n;
}
/****************************************************************************/
double correlation(double x1[],double y1[],double z1[],
		   double x2[],double y2[],double z2[],int n)
{
  double R[3][3],RtR[3][3];
  double a,b,c,q,r,theta,detR,w[3];
  double pi2=2*acos(-1.);
  int i,j,im;
 
  for (i=0;i<3;i++) for (j=0;j<3;j++) R[i][j]=0;
  for (i=0;i<n;i++) {
    R[0][0]+=x1[i]*x2[i]; R[0][1]+=x1[i]*y2[i]; R[0][2]+=x1[i]*z2[i];
    R[1][0]+=y1[i]*x2[i]; R[1][1]+=y1[i]*y2[i]; R[1][2]+=y1[i]*z2[i];
    R[2][0]+=z1[i]*x2[i]; R[2][1]+=z1[i]*y2[i]; R[2][2]+=z1[i]*z2[i];
  }
  detR=+R[0][0]*R[1][1]*R[2][2]
       +R[0][1]*R[1][2]*R[2][0]
       +R[0][2]*R[1][0]*R[2][1]
       -R[2][0]*R[1][1]*R[0][2]
       -R[2][1]*R[1][2]*R[0][0]
       -R[2][2]*R[1][0]*R[0][1];

  if (detR==0) printf("****** WARNING! *******\ndetR=0\n");

  for (i=0;i<3;i++) {
    for (j=0;j<3;j++) {
      RtR[i][j]=R[0][i]*R[0][j]+R[1][i]*R[1][j]+R[2][i]*R[2][j];
    }
  }

  a=-(RtR[0][0]+RtR[1][1]+RtR[2][2]);
  b=RtR[0][0]*RtR[1][1]+RtR[0][0]*RtR[2][2]+RtR[1][1]*RtR[2][2]-
    RtR[0][1]*RtR[0][1]-RtR[0][2]*RtR[0][2]-RtR[1][2]*RtR[1][2];
  c=-detR*detR;

  q=(a*a-3*b)/9;
  r=(2*a*a*a-9*a*b+27*c)/54;
  if (r*r/(q*q*q)>=1) {
    printf("error rmsd calc r^2=%e q^3=%e %e\n",r*r,q*q*q,r*r/(q*q*q));
    return 0;
  }
  q=sqrt(q);
  theta=acos(r/(q*q*q));
  w[0]=sqrt(fabs(-2*q*cos(theta/3)-a/3));
  w[1]=sqrt(fabs(-2*q*cos((theta+pi2)/3)-a/3));
  w[2]=sqrt(fabs(-2*q*cos((theta-pi2)/3)-a/3));
  
  if (detR<0) {
    im = w[1] > w[0] ? 1 : 0;
    return 2.*(w[im]+fabs(w[(im+1)%3]-w[(im+2)%3]))/n;
  } else
    return 2.*(w[0]+w[1]+w[2])/n;
}
/****************************************************************************/
double correlation2(double x1[],double y1[],double z1[],
		    double x2[],double y2[],double z2[],
		    double A[3][3],int n)
{
  double R[3][3],RRt[3][3],RtR[3][3];
  double w[3],u[3][3],v[3][3];
  double a,b,c,q,r,theta,detR;
  double ex,ey,e,tmp;
  double pi2=2*acos(-1.);
  int i,j,k,im;
  int s[3];
  
  for (i=0;i<3;i++) 
    for (j=0;j<3;j++) 
      R[i][j]=RtR[i][j]=RRt[i][j]=A[i][j]=0;

  for (i=0;i<n;i++) {
    R[0][0]+=x2[i]*x1[i]; R[0][1]+=x2[i]*y1[i]; R[0][2]+=x2[i]*z1[i];
    R[1][0]+=y2[i]*x1[i]; R[1][1]+=y2[i]*y1[i]; R[1][2]+=y2[i]*z1[i];
    R[2][0]+=z2[i]*x1[i]; R[2][1]+=z2[i]*y1[i]; R[2][2]+=z2[i]*z1[i];
  }
  detR=+R[0][0]*R[1][1]*R[2][2]
       +R[0][1]*R[1][2]*R[2][0]
       +R[0][2]*R[1][0]*R[2][1]
       -R[2][0]*R[1][1]*R[0][2]
       -R[2][1]*R[1][2]*R[0][0]
       -R[2][2]*R[1][0]*R[0][1];
  if (detR==0) printf("*** detR=0 ***\n");

  for (i=0;i<3;i++)
    for (j=0;j<3;j++)
      for (k=0;k<3;k++) {
	RtR[i][j]+=R[k][i]*R[k][j];
	RRt[i][j]+=R[i][k]*R[j][k];
      }

  /* eigenvalues */
  a=-(RtR[0][0]+RtR[1][1]+RtR[2][2]);
  b=+RtR[0][0]*RtR[1][1]+RtR[0][0]*RtR[2][2]+RtR[1][1]*RtR[2][2]
    -RtR[0][1]*RtR[0][1]-RtR[0][2]*RtR[0][2]-RtR[1][2]*RtR[1][2];
  c=-detR*detR;

  q=(a*a-3*b)/9;
  r=(2*a*a*a-9*a*b+27*c)/54;
  if (r*r/(q*q*q)>=1) {
    printf("error r^2=%e q^3=%e %e\n",r*r,q*q*q,r*r/(q*q*q));
    exit(-1);
  }
  q=sqrt(q);
  theta=acos(r/(q*q*q));
  
  w[0]=-2*q*cos(theta/3)-a/3;
  w[1]=-2*q*cos((theta+pi2)/3)-a/3;
  w[2]=-2*q*cos((theta-pi2)/3)-a/3;

  /* eigenvectors */
  for (i=0;i<3;i++) {
    for (k=0;k<3;k++) {RtR[k][k]-=w[i]; RRt[k][k]-=w[i];}
    ey=RtR[0][2]*RtR[1][2]-RtR[0][1]*RtR[2][2];
    ey/=RtR[0][1]*RtR[1][2]-RtR[0][2]*RtR[1][1];
    ex=-(RtR[0][1]*ey+RtR[0][2])/RtR[0][0];
    e=sqrt(ex*ex+ey*ey+1.);
    v[0][i]=ex/e;
    v[1][i]=ey/e;
    v[2][i]=1./e;
    ey=RRt[0][2]*RRt[1][2]-RRt[0][1]*RRt[2][2];
    ey/=RRt[0][1]*RRt[1][2]-RRt[0][2]*RRt[1][1];
    ex=-(RRt[0][1]*ey+RRt[0][2])/RRt[0][0];
    e=sqrt(ex*ex+ey*ey+1.);
    u[0][i]=ex/e;
    u[1][i]=ey/e;
    u[2][i]=1./e;
    for (k=0;k<3;k++) {RtR[k][k]+=w[i]; RRt[k][k]+=w[i];}
  }

  /* orientation */
  for (i=0;i<3;i++) {
    tmp=R[i][0]*v[0][i]+R[i][1]*v[1][i]+R[i][2]*v[2][i];
    s[i]=tmp/u[i][i] > 0 ? 1 : -1;
  }
  if (w[1] < w[2]) im=1; else im=2;
  if (w[0] < w[im]) im=0;
  if (detR<0) s[im]*=-1;

  /* rotation matrix */
  for(i=0;i<3;i++)
    for(j=0;j<3;j++)
      for (k=0;k<3;k++) 
	A[i][j]+=v[i][k]*s[k]*u[j][k];

  w[0]=sqrt(w[0]);
  w[1]=sqrt(w[1]);
  w[2]=sqrt(w[2]);
  if (detR<0) 
    return 2.*(w[(im+1)%3]+w[(im+2)%3]-w[im])/n;
  else
    return 2.*(w[0]+w[1]+w[2])/n;
}
/****************************************************************************/
double genrand(void) {

  return ran3n(&seed);
   //  return genrand64_real3(void);
}
/****************************************************************************/
/******** RAN3N RANDOM NUMBER GENERATION ************************************/
/****************************************************************************/

#define MBIG 1000000000
#define MSEED 161803398
#define MZ 0
#define FAC (1.0/MBIG)

double ran3n(long *idum)
/*******************************************************/
/*  returns uniform deviate r, 0<r<1                   */
/*  NR sec 7.1,  changed 1) float->double 2) r>0       */
/*******************************************************/
{
  static int inext,inextp;
  static long ma[56];
  static int iff=0;
  long mj,mk;
  int i,ii,k;
  double ret_val;
  
  if (*idum < 0 || iff == 0) {
    iff=1;
    mj=MSEED-(*idum < 0 ? -*idum : *idum);
    mj %= MBIG;
    ma[55]=mj;
    mk=1;
    for (i=1;i<=54;i++) {
      ii=(21*i) % 55;
      ma[ii]=mk;
      mk=mj-mk;
      if (mk < MZ) mk += MBIG;
      mj=ma[ii];
    }
    for (k=1;k<=4;k++)
      for (i=1;i<=55;i++) {
	ma[i] -= ma[1+(i+30) % 55];
	if (ma[i] < MZ) ma[i] += MBIG;
      }
    inext=0;
    inextp=31;
    *idum=1;
  }
  else {(*idum)++;}
  if (++inext == 56) inext=1;
  if (++inextp == 56) inextp=1;
  mj=ma[inext]-ma[inextp];
  if (mj < MZ) mj += MBIG;
  ma[inext]=mj;
  ret_val = mj*FAC;
  if (mj == 0) ret_val = FAC;
  if (mj == MBIG) ret_val = 1. - FAC;
  return ret_val;
}
#undef MBIG
#undef MSEED
#undef MZ
#undef FAC

/****************************************************************************/
/* (C) Copr. 1994 Feb 4th Jan 2:27pm  Potti-Soft  registered Trademark      */
/****************************************************************************/

/****************************************************************************/
/******** RAN3N RANDOM NUMBER GENERATION ************************************/
/****************************************************************************/
/* 
   A C-program for MT19937-64 (2014/2/23 version).
   Coded by Takuji Nishimura and Makoto Matsumoto.

   This is a 64-bit version of Mersenne Twister pseudorandom number
   generator.

   Before using, initialize the state by using init_genrand64(seed)  
   or init_by_array64(init_key, key_length).

   Copyright (C) 2004, 2014, Makoto Matsumoto and Takuji Nishimura,
   All rights reserved.                          

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions
   are met:

     1. Redistributions of source code must retain the above copyright
        notice, this list of conditions and the following disclaimer.

     2. Redistributions in binary form must reproduce the above copyright
        notice, this list of conditions and the following disclaimer in the
        documentation and/or other materials provided with the distribution.

     3. The names of its contributors may not be used to endorse or promote 
        products derived from this software without specific prior written 
        permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
   "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
   LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
   A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
   CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
   PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
   PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
   LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
   SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

   References:
   T. Nishimura, ``Tables of 64-bit Mersenne Twisters''
     ACM Transactions on Modeling and 
     Computer Simulation 10. (2000) 348--357.
   M. Matsumoto and T. Nishimura,
     ``Mersenne Twister: a 623-dimensionally equidistributed
       uniform pseudorandom number generator''
     ACM Transactions on Modeling and 
     Computer Simulation 8. (Jan. 1998) 3--30.

   Any feedback is very welcome.
   http://www.math.hiroshima-u.ac.jp/~m-mat/MT/emt.html
   email: m-mat @ math.sci.hiroshima-u.ac.jp (remove spaces)
*/


// #include <stdio.h>
#include "mt64.h"

#define NN 312
#define MM 156
#define MATRIX_A UINT64_C(0xB5026F5AA96619E9)
#define UM UINT64_C(0xFFFFFFFF80000000) /* Most significant 33 bits */
#define LM UINT64_C(0x7FFFFFFF) /* Least significant 31 bits */


/* The array for the state vector */
static uint64_t mt[NN]; 
/* mti==NN+1 means mt[NN] is not initialized */
static int mti=NN+1; 

/* initializes mt[NN] with a seed */
void init_genrand64(uint64_t seed)
{
    mt[0] = seed;
    for (mti=1; mti<NN; mti++) 
        mt[mti] =  (UINT64_C(6364136223846793005) * (mt[mti-1] ^ (mt[mti-1] >> 62)) + mti);
}

/* initialize by an array with array-length */
/* init_key is the array for initializing keys */
/* key_length is its length */
void init_by_array64(uint64_t init_key[],
		     uint64_t key_length)
{
    unsigned int i, j;
    uint64_t k;
    init_genrand64(UINT64_C(19650218));
    i=1; j=0;
    k = (NN>key_length ? NN : key_length);
    for (; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * UINT64_C(3935559000370003845)))
          + init_key[j] + j; /* non linear */
        i++; j++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
        if (j>=key_length) j=0;
    }
    for (k=NN-1; k; k--) {
        mt[i] = (mt[i] ^ ((mt[i-1] ^ (mt[i-1] >> 62)) * UINT64_C(2862933555777941757)))
          - i; /* non linear */
        i++;
        if (i>=NN) { mt[0] = mt[NN-1]; i=1; }
    }

    mt[0] = UINT64_C(1) << 63; /* MSB is 1; assuring non-zero initial array */ 
}

/* generates a random number on [0, 2^64-1]-interval */
uint64_t genrand64_int64(void)
{
    int i;
    uint64_t x;
    static uint64_t mag01[2]={UINT64_C(0), MATRIX_A};

    if (mti >= NN) { /* generate NN words at one time */

        /* if init_genrand64() has not been called, */
        /* a default initial seed is used     */
        if (mti == NN+1) 
            init_genrand64(UINT64_C(5489)); 

        for (i=0;i<NN-MM;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+MM] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
        }
        for (;i<NN-1;i++) {
            x = (mt[i]&UM)|(mt[i+1]&LM);
            mt[i] = mt[i+(MM-NN)] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];
        }
        x = (mt[NN-1]&UM)|(mt[0]&LM);
        mt[NN-1] = mt[MM-1] ^ (x>>1) ^ mag01[(int)(x&UINT64_C(1))];

        mti = 0;
    }
  
    x = mt[mti++];

    x ^= (x >> 29) & UINT64_C(0x5555555555555555);
    x ^= (x << 17) & UINT64_C(0x71D67FFFEDA60000);
    x ^= (x << 37) & UINT64_C(0xFFF7EEE000000000);
    x ^= (x >> 43);

    return x;
}

/* generates a random number on [0, 2^63-1]-interval */
int64_t genrand64_int63(void)
{
    return (int64_t)(genrand64_int64() >> 1);
}

/* generates a random number on [0,1]-real-interval */
double genrand64_real1(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740991.0);
}

/* generates a random number on [0,1)-real-interval */
double genrand64_real2(void)
{
    return (genrand64_int64() >> 11) * (1.0/9007199254740992.0);
}

/* generates a random number on (0,1)-real-interval */
double genrand64_real3(void)
{
    return ((genrand64_int64() >> 12) + 0.5) * (1.0/4503599627370496.0);
}


void save_twister(char *fn) {
  FILE *fp;
  
  if ( (fp = fopen(fn,"w")) != NULL) {
    fwrite(&mti,sizeof(int),1,fp);
    fwrite(mt,sizeof(uint64_t),NN,fp);
    fclose(fp);
  }

}


void restore_twister(char *fn) {
  FILE *fp;
  
  if ( (fp = fopen(fn,"r")) != NULL) {
    fread(&mti,sizeof(int),1,fp);
    fread(mt,sizeof(uint64_t),NN,fp);
    fclose(fp);
  }

}
