/************* MC parameters ************************************************/
double beta[NTMP];                 /* inverse temperatures                  */
double g[NTMP];                    /* sim. temp. parameters                 */
int ind;                           /* index of dynamical parameter          */
const double abgs=300.0;           /* bgs parameter                         */
const double bbgs=10.0;            /* bgs parameter                         */
/************* interactions *************************************************/
const double epshb=3.22;           /* depth H bond                          */
//const double epshb=0;           /* depth H bond                          */
const double sighb=2.0;            /* min H bond                            */
const double cuthb=4.5;            /* cutoff H bond                         */
const double sigsa[NA]={1.75,      /* C                                     */
	                1.55,      /* N                                     */
                        1.42,      /* O                                     */
                        1.00,      /* H                                     */
                        1.00,      /* amide H                               */
			1.75};     /* Cb                                    */
const double ksa=0.10;             /* strength sa                           */
const double cut=6.6;              /* cutoff sa                             */
const double epshp=0.805;          /* strength hp --- ORIG VALUE (0.805)    */
//const double epshp=0;          /* strength hp --- ORIG VALUE (0.805)    */
const double sighp=5.0;            /* min hp interaction                    */
const double dchp=4.0;             /* cutoff hp = sighp+dchp                */
const double kbias=50.0;           /* bias strength                         */
const double eclash=1e6;           /* crowder clash penalty                 */
const double e_max=1e5;            /* crowder clash penalty                 */
double gam[N][N];                  /* scale factors hb energy               */
int fixed[NCH];                    /* >0 no updates                         */
/****************************************************************************/
double E,Ebias,Ehb,Eev,Eloc,Ehp,Ecc,Ecp; /* energies                        */
double Eo,Ebiaso,Ehbo,Eevo,Eloco,Ehpo,Ecco,Ecpo; /* energies                        */
int atyp[NTO];                     /* atom type                             */
short aa[NTO][NTO];                /* type of atom-atom pair                */
/************* crowder interactions *****************************************/
const double rcrowd = 12.0;        /* Radius of crowders                    */
const double sigcr = 3.0;          /* Softness crowders                     */
const double epsilonrep = 1.0;     /* Strength repulsion                    */
const double epsilonatt = 1.5;    /* Strength attraction                   */
const double cut_crowd_crowd  = 2 * rcrowd + sigcr;               /* cutoff */
const double cut_crowd_crowd2 = cut_crowd_crowd * cut_crowd_crowd;
const double cut_crowd_atom   = rcrowd + 2 + sigcr;               /* cutoff */
const double cut_crowd_atom2  = cut_crowd_atom * cut_crowd_atom;
/************* chain geometry ***********************************************/
const double b[3]={1.33,1.46,1.52};/* backbone CN,NCa,CaC bond lengths      */
const double bNH=1.00;             /* NH bond length                        */
const double bCaHa=1.08;           /* CaHa bond length                      */
const double bCaCb=1.53;           /* CaCb bond length                      */
const double bCO=1.23;             /* CO bond length                        */
/****************************************************************************/
int seq[N];                        /* sequence                              */
double ph[NBA];                    /* all backbone torsional angles         */
double x[NTO],y[NTO],z[NTO];       /* atom coordinates                      */
double xb[NTO],yb[NTO],zb[NTO];    /* within boundary box [0,BOX]           */
double xo[NTO],yo[NTO],zo[NTO];    /* backup                                */
double xcr[NCR],ycr[NCR],zcr[NCR]; /* crowder coordinates                   */
double xcrb[NCR],ycrb[NCR],zcrb[NCR]; /* within boundary box [0,BOX]        */
double xcro[NCR],ycro[NCR],zcro[NCR]; /* backup                             */
double xnat[N],ynat[N],znat[N]; 
double xnat2[N],ynat2[N],znat2[N]; 
/************ indices *******************************************************/
int iN[N];                         /* N atoms                               */
int iCa[N];                        /* Ca atoms                              */
int iC[N];                         /* C atoms                               */
int Nb[NCH],Ne[NCH];               /* first, last aa in chain               */
int iBeg[NCH],iEnd[NCH];           /* first, last atom in chain             */
int a2c[N];                        /* aa to chain                           */
int i2a[NTO];                      /* atom to aa                            */
int idonacc[NTO];                  /* donors and acceptors                  */
/************* misc *********************************************************/
const double boxhf=BOX/2.;         /* half box length                       */
const double vbox=BOX*BOX*BOX;     /* box volume                            */
const double hb_lim_per_aa = 3*epshb;
double pi,pi2,pid2;                /* pi,2*pi,pi/2                          */
long seed=-15309167;               /* random number seed                    */
uint64_t seed_mt;                  /* random number seed, Marsenne Twister  */
long orig_seed;
int relax=0;
long int imc0;
double so[NTMP][NOBS],son[NTMP][NOBS],sou[NTMP][NOBS],sot[NTMP][NOBS];
double minobs[NOBS];
double Phic;
/************* obs **********************************************************/
int contact_pair_i[MAXCONT];	
int contact_pair_j[MAXCONT];
int boxnum[MAXCONT];
int ncontacts;
/****************************************************************************/
