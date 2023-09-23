/*******************************************************************/
/* calculates constants needed by pdaf.c                            */
/*******************************************************************/
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

int seq[500];                      /* current amino acid sequence           */
int des[500];                      /* variable amino acid positions         */
int Nb[500],Ne[500];               /* first, last aa in chain               */

int main(int argc,char *argv[])
{
  int i,j,n=0,nto=0,ndf=0,nrt=0;
  int c,NCH=1,NCR=0;
  FILE *fp;
  
  if (argc == 1) {
    printf("Usage: Give input filename, #chains (NCH) and #crowders (NCR).\n");
    printf("Output: File sys.h with constants N, NTO, NDF, etc\n");
    exit(-1);
  }

  printf("argc %i\n",argc);
  
  if (argc >= 3) NCH = atoi(argv[2]);
  if (argc == 4) NCR = atoi(argv[3]);
 
  fp=fopen("input","r");   
  for (i=j=0;i<NCH;i++) {
    Nb[i]=j;                   /* initial amino acid sequence */
    while ((c=getc(fp))!='\n') seq[j++]=c; 
    Ne[i]=j-1;
  }
  fclose(fp);

  for (i=0;i<NCH;i++) {
    printf("\nChain   First      Last     Sequence\n");
    printf("%3i %9i %9i     ",i,Nb[i],Ne[i]); 
    for (j=0;j<=Ne[i]-Nb[i];j++) {
      printf("%c",seq[j+Nb[i]]); 
      if ((j+1)%5==0) putchar(' ');
    }
    printf("\n");
  }
  printf("\n");

  for (i=0;i<NCH;i++) {
    for (j=Nb[i];j<=Ne[i];j++) {
      switch (seq[j]) {
      case 'S' : {
	nto+=7;
	ndf+=2;
	nrt++;
	break;
      }
      case 'L' : {
	nto+=7;
	ndf+=2;
	nrt++;
	break;
      }
      default : {
	nto+=7;
	ndf+=2;
	nrt++;
	break;
      }
      }
      n++;
    }
  }
  
  printf("\n");

  printf("NCH=%i N=%i NTO=%i NDF=%i NRT=%i NCR=%i\n",NCH,n,nto,ndf,nrt,NCR);

  fp=fopen("sys.h","w");

  fprintf(fp,"/************* geometry ************/\n");
  fprintf(fp,"# define N %i     /* # amino acids */\n",n);
  fprintf(fp,"# define NCH %i   /* # chains */\n",NCH);
  fprintf(fp,"# define NTO %i   /* # atoms (total) */\n",nto);
  fprintf(fp,"# define NDF %i   /* # backbone dofs */\n",ndf);
  fprintf(fp,"# define NRT %i   /* # rotamer dofs */\n",nrt);
  fprintf(fp,"# define NCR %i   /* # crowders */\n",NCR);
  
  fclose(fp);

  return 0;
}
/**************************************************************/ 
