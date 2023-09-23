# include <inttypes.h>
/************* interactions *************************************************/
extern double E,Ebias,Ehb,Eev,Eloc,Ehp,Ecc,Ecp;
extern double Eo,Ebiaso,Ehbo,Eevo,Eloco,Ehpo,Ecco,Ecpo; 
extern double Mhp[][NCH];   
extern double Mhb[][NCH];   
extern const double epshb;     
extern const double sighb;         
extern const double cuthb;         
extern const double sigsa[];
extern const double ksa;
extern const double cut;
extern const double epshp;
extern const double sighp;
extern const double dchp;
extern const double kbias;
extern const double kcr;
extern const double Phi;             /*  Volume fraction of crowders */
extern const double rcrowd;          /*  Radius of crowders  */
extern const double sigcr;           /*  Coefficient of softness */
extern const double epsilonrep;      /*  Strength of crowding interactions repulsion */
extern const double epsilonatt;      /*  Strength of crowding interactions attraction */
extern const double cut_crowd_crowd,cut_crowd_crowd2;
extern const double cut_crowd_atom, cut_crowd_atom2;
extern double gam[][N];
extern int atyp[];      
extern short aa[][NTO];
extern const double eclash,e_max;
extern int relax;
/************* chain geometry ***********************************************/
extern int seq[];       
extern int iN[];        
extern int iCa[];       
extern int iC[];        
extern const double b[];
extern const double bNH;
extern const double bCaHa;
extern const double bCaCb;
extern const double bCO;
extern double x[],y[],z[]; 
extern double xb[],yb[],zb[]; 
extern double xo[],yo[],zo[]; 
extern double xcr[],ycr[],zcr[]; 
extern double xcrb[],ycrb[],zcrb[]; 
extern double xcro[],ycro[],zcro[]; 
extern double xnat[],ynat[],znat[]; 
extern double xnat2[],ynat2[],znat2[]; 
extern double xpdb[],ypdb[],zpdb[]; 
extern double ph[]; 
/************* monte carlo **************************************************/
extern double beta[];     
extern double g[];        
extern int ind;           
extern const double abgs;
extern const double bbgs;
extern const double boxhf;
extern const double vbox;
/************* crowders *********************************************************/
extern double xcr[],ycr[],zcr[];
/************* misc *********************************************************/
extern double pi,pi2,pid2;
extern long seed,orig_seed;
extern uint64_t seed_mt;        
extern int a2c[],i2a[];
extern int idonacc[];
extern int Nb[],Ne[];  
extern int iBeg[],iEnd[];
extern long int imc0;
extern double so[][NOBS];
extern double minobs[NOBS];
extern double Phic;
extern const double hb_lim_per_aa;
/************** obs *********************************************************/
extern int contact_pair_i[];	
extern int contact_pair_j[];
extern int boxnum[];
extern int ncontacts;
/****************************************************************************/
/****************************************************************************/
/****************************************************************************/
/**** utils.c ****/
double genrand(void);
void save_twister(char *fn);
void restore_twister(char *fn);
int readpdb(int flag,char *fn,double x[],double y[],double z[],int n1,int n2);
void dumppdb(char *fn,int mc,int nobs,double o[],char format[100]);
double rmsd2(double x1[],double y1[], double z1[],
	     double x2[],double y2[], double z2[],int n);
double rmsd2_nopt(double x1[],double y1[],double z1[],
		  double x2[],double y2[],double z2[],int n);
double ran3n(long *seed);
int get_atm(int iaa,char *name);
void get_atomname(int a,char *name);
void get_aminoname(int ia,char *name);
void printatom(int iatm);
void printatom2(int iatm,int jatm);
void superimpose(double *xt,double *yt,double *zt,int a1,int a2);
double center_of_mass(double x[],double y[],double z[],
		      double *xcm,double *ycm,double *zcm,int n);
double correlation(double x1[],double y1[],double z1[],
		   double x2[],double y2[],double z2[],int n);
double correlation2(double x1[],double y1[],double z1[],
		    double x2[],double y2[],double z2[],
		    double A[3][3],int n);
/**** geometry.c ****/
void int2cart(int ia,int dir); 
void int2cart_all(int dir); 
double rgyr(int ic);
double cm(double *xcm,double *ycm,double *zcm,int ic,int offset);
void in2box(double *xb,double *yb,double *zb,double *x,double *y,double *z,
	    int n1,int n2);
void movech(double dx,double dy,double dz,int ic);
void bc(double *x);
double dist(int i1,int i2);
double dist_atom_crowd(int i1,int i2);
double dist_cr(int ic1,int ic2);
void backup_coord();
/**** energy.c ****/
void ecalc();
double crowd_crowd_ecalc(int iflag,double r2);
double crowd_atom_ecalc(int iflag,int atmtyp,double r2);
double crowd_atom_attr_ecalc(int iflag,int atmtyp,char aatyp,double r2);
double hp_ecalc(int iflag,int l,int m);
double hbonds_ecalc(int iflag,int i,int j,double r2);
double exvol_ecalc(int iflag,int a,double r2);
void localsa(char *move,char *mc,int a1,int a2,double *esa);
void bias(char *move,char *mc,int a1,int a2,double *ebi);
void hp(char *move,char *mc,int a1,int a2,double *ehp);
void cellcalc(char *move,char *mc,int icr,int a1,int a2,double *eev,double *ehb,
	      double *ecc,double *ecp,double emax);
void listcalc_crowd(int iflag,int inc,short list[],int nlc,short listc[],
		    double mcp[],double *ecp,double emax);
void listcalc(int iflag,int in_cell,int nl,short list[],
	      double *eev,double *ehb);
int sel_aa(char *mov);
int sel_aa2(int iflag,int a1,int a2);
void pivot(int iflag,long *nup,double *accup,long *ndw,double *accdw);
int pivup(int iflag,int ia);
int pivdw(int iflag,int ia);
void bgs(int iflag,long *nup,double *accup,long *ndw,double *accdw);
int bgsup(char *move, int ia1);
int bgsdw(char *move, int ia2);
int bgsdw_alt(char *move,int ia1);
void qsi(long *nup,double *accup,long *ndw,double *accdw);
int qsiup(int iflag, int ia1);
int qsidw(int iflag, int ia2);
int trans(int iflag);
int rotate(int iflag);
int rotate2(int iflag);
int crowd_translate();
int flip(double e);
int get_cell(int a,int *pnt);
void get_cell_xyz(short *ix,short *iy,short *iz,int ic,int ns);
void remove_from_list(int ic,int *pnt,short *cell,int aout);
void remove_from_list_range(int ic,int *pnt,short *cell,int a1,int a2);
void insert_into_list(int ic,int *pnt,short *cell,int ain);
double dist2_cell_crowd(short ix, short iy,short iz,int icr,double cutg);
void build_list_crowd(short ix,short iy,short iz,int *nlc,short *listc,double cutg);
void build_list_inc(int a,int *pnt,int *nl,short *list,
		    int sel,int a1,int a2);
void build_list_fwd(int ic,int *pnt,short *cell,
		    short ix,short iy,short iz,
		    int *nl,short *list,int ns,
		    int sel,int a1,int a2);
void build_list_bwd(int ic,int *pnt,short *cell,
		    short ix,short iy,short iz,
		    int *nl,short *list,int ns,
		    int sel,int a1,int a2);
int crowd_clash();
/**** obs.c ****/
//double hbonds_intra_chain(int ic);
//double hp_intra_chain(int iflag,int ic);
double mindist_cr_at();
double mindist_cr_cr();
void crowd_stat(int iflag);
double rgyr_from_coord(double *xt,double *yt,double *zt,int a1,int a2);
int depth(int clr_no,int icur,int *alist);
int clustering();
void hist_sq(int iflag,int ind,int icr);
void pair_correlation_crcr(int iflag,int ind);
void pair_correlation_crat(int iflag,int iatom,int ind);
void pair_correlation_cm(int iflag,int ich,int ind); 
double ree(int ich);
int dist2(int x1,int y1,int z1,int x2,int y2,int z2);
int cont_made(int iflag,int i, int j);
int contacts();
int contacts_all();
int contacts_in_box(int b);
int contacts_res(int iaa);
int hhcont(int c1,int c2);
double calcrmsd_ca(int iflag,double *xt,double *yt,double *zt,int a1,int a2);
double cmdist(int iflag,int ic1,int ic2);
int helix(int ic);
int sheet(int ic);
int helix_aa(int i);
void histo_cont(int iflag,int ind,int n);
void histo_cont_all(int iflag,int ind,int n);
void histo_cont_rg(int iflag,int itmp,int cont,double x);
void histo_contmap(int iflag,int ind);
void helixprof(int iflag,int itmp);
void histoe(int iflag,int itmp,double x);
void histod1(int iflag,int itmp,double x);
void histod2(int iflag,int itmp,double x);
void histod3(int iflag,int itmp,double x);
void histod4(int iflag,int itmp,double x);
void histoed1(int iflag,int itmp,double x1,double x2);
void histoed2(int iflag,int itmp,double x1,double x2);
void histoed3(int iflag,int itmp,double x1,double x2);
void histodd1(int iflag,int itmp,double x1,double x2);
void histodd2(int iflag,int itmp,double x1,double x2);
void histonn1(int iflag,int itmp,int m,int n);
double interhp();
double interhb();
double gyros(void);
double delta2(int flag);
void rama();
/**** misc.c ****/
double settings_read(char str_target[]);
void printheading(double frcrd[],double frpiv[],double frbgs[],double frtrl[]);
void write_conf_ver2(int iflag,char *fn,char *fmode);
void write_conf(char *fn,char *fmode);
void write_checkpnt(long *imc,double so[][NOBS],double mino[]);
void read_checkpnt(long *imc,double so[][NOBS],double mino[]);
int read_conf_ver2(int iflag,char *fn,int m);
int read_conf(int iflag,char *fn,int m);
void plot_rama(int ia,int ch);
void runtime(long it,double o[],int n);
void runtime_allcols(long it,double o[],int n);
void ch2box();
void cr2box();
void init(void);
void write_averages(char *fn,int ntmp,int nobs,double so[][NOBS]);
int acceptance(long *nflip,double *accflip,
	       long *ncrd,double *acccrd,
	       long *npivup,double *accpivup,
	       long *npivdw,double *accpivdw,
	       long *nbgsup,double *accbgsup,
	       long *nbgsdw,double *accbgsdw,
	       long *ntrl,double *acctrl);
/****************************************************************************/



