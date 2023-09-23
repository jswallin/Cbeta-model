/************* chain geometry ***********************************************/
# define NBA (3*N)                 /* # backbone atoms                      */
# define NA 6                      /* # atom types                          */
# define BOX 100                   /* box length of periodic boundary       */
/************* MC parameters ************************************************/
# define MCCYC 500000000           /* # cycles                              */
# define MCCON 100                 /* # steps per cycle                     */
# define NTMP 1                   /* # temperatures                        */
# define TMAX 0.53                 /* max temperature                       */
# define TMIN 0.53                 /* min temperature                       */
# define NTHERM 10000            /* # discarded steps                     */
# define INORM 10                 /* normalize                             */
# define IRT 1000                  /* run time                              */
# define IRAMA 100000              /* Ramachandran data                     */
# define ICONF 1000                /* configurations                        */
# define ICHECK 10000              /* write checkpoint data                 */
# define ISTART 0                  /* 0 hot, 1 helix, 2 conf, 3 checkpoint  */
# define IREADG 1                  /* 1 read g, else g=0                    */
# define ISEED 1                   /* 1 randomize seed (/dev/urandom)       */
# define PIVSTP 6.28               /* step size pivot                       */
# define TRLSTP 2.00               /* step size translation                 */
# define RGRSTP 0.628              /* step size rigid rotation              */
# define CRSTEP 5.00               /* step size crowder translation         */
/************* interactions *************************************************/
# define POWA 0.5                  /* powers in hbonds                      */
# define POWB 0.5                  /* powers in hbonds                      */
# define ALPHA1 0.75               /* scale factor nonlocal pairs           */
# define ALPHA2 1.25               /* scale factor nonlocal pairs           */
# define ISOFT 0                   /* 1 soft attraction on, 0 off           */
/************* rmsd calculation *********************************************/
# define RMSD 1                    /* >0 on; 0 off                          */
# define NR1 0                     /* read NATIVE1/NATIVE2 from NR1 to NR2  */
# define NR2 34                    /* store in xnat1[], xnat2[], etc        */
/************* misc *********************************************************/
# define NOBS 13                   /* # observables                         */
# define MAXCELL 100000            /* maximum # cells                       */
# define MAXLOC (100*N)            /* maximum # local pairs                 */
# define NBIN 100                  /* # bins in histograms                  */
# define NBINL 500                 /* # bins in histograms (large)          */
# define NB 5000
# define MAXCONT 600
# define RAN_RAN3 1
# define RAN_MARSENNE_TWISTER 0
/************* INPUT ********************************************************/
# define NATIVE1 "./A35_emin.pdb"
# define NATIVE2 "./A35_emin.pdb"
# define STARTCONF "./2h_1h_mine.conf"
# define INPUT "input" 
# define INPUTG "inputg" 
# define CONTACTS "./native_contacts_A35_boxnum" 
/************* OUTPUT *******************************************************/
# define RT "tmp_rt" 
# define CONF "tmp_conf" 
# define RAMA "tmp_rama" 
# define AVERAGES "averages"
# define OUTPDB "_pdb"
# define OUTCONF "_conf"
# define OUTPUTG "outputg" 
# define CHECKCONF "./checkpnt/_conf"
# define CHECKDATA "./checkpnt/_data"
# define NATIVE "./native.pdb"
# define ACCSTAT "./acceptance"
/************* functions ****************************************************/
# define max(A,B) ((A)>(B)?(A):(B)) 
# define min(A,B) ((A)<(B)?(A):(B))
/****************************************************************************/
