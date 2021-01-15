/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   https://lammps.sandia.gov/, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Carsten Svaneborg (SDU)
------------------------------------------------------------------------- */

/* 8
   01-09-21
   coeff does not set masses
   exterminate coeff :(
   unexterminate coeff and debug :|
   remove PairZero::coeff text
   now using cout to print
   added <iostream>
   manually setting masses because I do not fear a deity
   FINISH 8, MASSES ASSIGNED, ATOMS JUMPING OUT OF BOX
   9
   01-09-21
   set cutoff manually
   FUNCTIONS
   10
   01-09-21
   try to implement box2lattice
*/

extern "C" {
  #include "mlib.h"
}
#include "pair_maise.h"

#include "atom.h"
#include "comm.h"
#include "memory.h"
#include "error.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "lattice.h"
#include "domain.h"
#include "compute.h"
#include "compute_ke.h"
#include "modify.h"

#include <cstring>
#include <iostream>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMaise::PairMaise(LAMMPS *lmp) : Pair(lmp) {
  coeffflag=1;
  writedata=1;
  single_enable=1;
  respa_enable=1;
  scale = nullptr;

  j = 0;
  
  comm_forward = 38;
  comm_reverse = 30;
}

/* ---------------------------------------------------------------------- */

PairMaise::~PairMaise()
{
  if (allocated) {
    memory->destroy(setflag);
    memory->destroy(cutsq);
    memory->destroy(cut);
    memory->destroy(scale);
  }
  delete [] map;
}

/* ---------------------------------------------------------------------- */

void PairMaise::compute(int eflag, int vflag)
{
  int i;
  int q;
  
  ev_init(eflag,vflag);

  double **f = atom->f;
  ANN mR;
  PRS mP;
  PRS mW[9];
  LNK mL;
  Cell mC;
  int mCODE;
  int mN;
  int mNM;
  int mND;
  int mNP;
  int mXT;
  int mATMN[4];
  double mLAT[3][3];
  double mPOS[4][3];
  double **x;
  double **v;

  double mH;
  double mFRC[4][3];
  double mSTR[6];
  double mLATf[3*3];
  double mPOSf[4*3];
  double mFRCf[4*3];
  mCODE = 0;
  mN    = 4;
  mNM   = 500;
  mND   = 0;
  mNP   = 4;
  mXT   = 1;
  mC.p  = 0;
  x = atom->x;
  v = atom->v;

  double mx;
  double my;
  double mz;


  for(i=0;i<3;i++){
    if (i==0){
      mx = domain->boxhi[0] - domain->boxlo[0];
      my = 0;
      mz = 0;
    }
    if (i==1){
      mx = domain->xy;
      my = domain->boxhi[1] - domain->boxlo[1];
      mz = 0;
    }
    if (i==2){
      mx = domain->xz;
      my = domain->yz;
      mz = domain->boxhi[2] - domain->boxlo[2];
    }
    /*domain->lattice->box2lattice(mx,my,mz);*/
    mLAT[i][0] = mx;
    mLAT[i][1] = my;
    mLAT[i][2] = mz;
  }  

  double ke;
  mPOS[ 1][ 0] = -2.0;
  mPOS[ 3][ 1] = -2.0;

  mATMN[0] = 0; mATMN[1] = 0; mATMN[2] = 0; mATMN[3] = 0;
  mL.B = 0;
  for(i=0;i<3;i++)
    for(q=0;q<3;q++)
      mLATf[3*i+q] = mLAT[i][q];
  for(i=0;i<mN;i++)
    for(q=0;q<3;q++){
      if(x[i][q] >= 10)
	mPOSf[3*i+q] = x[i][q] - 20;
      else
	mPOSf[3*i+q] = x[i][q]; 
    }
  mH = CALL_MAISE(&mR, &mP, mW, &mL, &mC, mCODE, mN, mNM, mND, mNP, mXT, mATMN, mLATf, mPOSf, mFRCf, mSTR);
  for (int i = 0; i < mN; i++){
    ke +=0.5* mass[0]*(v[i][0]*v[i][0] + v[i][1]*v[i][1] + v[i][2]*v[i][2]);
  }
  eng_coul = ke-mC.E;
  if (vflag_fdotr) virial_fdotr_compute();
  for(i=0;i<mN;i++)
    for(q=0;q<3;q++)
      f[i][q] = mFRCf[3*i+q];
  if(j<=2)
  for(i=0;i<mN;i++)
    for(q=0;q<3;q++)
      printf("%d %d force: % lf % lf % lf\n",i,q,f[i][q],mPOSf[3*i+q],mC.E);
  j++;

}

/* ---------------------------------------------------------------------- */

void PairMaise::compute_outer(int eflag, int vflag)
{
 ev_init(eflag,vflag);
}

/* ----------------------------------------------------------------------
   allocate all arrays
------------------------------------------------------------------------- */

void PairMaise::allocate()
{
  allocated = 1;
  int n = atom->ntypes;
  map = new int[n+1];

  memory->create(setflag,n+1,n+1,"pair:setflag");
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  memory->create(cutsq,n+1,n+1,"pair:cutsq");
  memory->create(cut,n+1,n+1,"pair:cut");
  memory->create(scale,n+1,n+1,"pair:scale");
}

/* ----------------------------------------------------------------------
   global settings
------------------------------------------------------------------------- */

void PairMaise::settings(int narg, char **arg)
{
  if ((narg != 1) && (narg != 2))
    error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
  if (narg == 2) {
    if (strcmp("nocoeff",arg[1]) == 0) coeffflag=0;
    else error->all(FLERR,"Illegal pair_style command");
  }

  // reset cutoffs that have been explicitly set

  if (allocated) {
    int i,j;
    for (i = 1; i <= atom->ntypes; i++)
      for (j = i+1; j <= atom->ntypes; j++)
        cut[i][j] = cut_global;
  }
}

/* ----------------------------------------------------------------------
   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMaise::coeff(int narg, char **arg)
{
    int m,n;
    if ((narg < 2) || (coeffflag && narg > 4))
    error->all(FLERR,"Incorrect args for pair coefficients");

  if (!allocated) allocate();
  
  if (strcmp(arg[0], "*") != 0 || strcmp(arg[1], "*") != 0)
    error->all(FLERR, "Incorrect args for pair coefficients");

  if (nelements) {
    for (int i=0; i < nelements; i++) delete [] elements[i];
    delete [] elements;
    delete [] mass;
  }

  nelements = 2;
  
  elements = new char*[nelements];
  mass = new double[nelements];

  for (int i = 0; i < nelements; i++) {
    n = strlen(arg[i+2]) + 1;
    elements[i] = new char[n];
    strcpy(elements[i], arg[i+2]);
  }

  for (int i = 4 + nelements; i < narg; i++){
    m = i - (4+nelements) + 1;
    int j;
    for (j = 0; j < nelements; j++)
      if (strcmp(arg[i],elements[j]) == 0) break;
    if (j < nelements) map[m] = j;
    else if (strcmp(arg[i],"NULL") == 0) map[m] = -1;
    else error->all(FLERR,"Incorrect args for pair coefficients");
  }

  n = atom->ntypes;
  for (int i = 1; i <= n; i++)
    for (int j = i; j <= n; j++)
      setflag[i][j] = 0;

  mass[0] = 63.546;

    for (int i = 0; i < nelements; i++)
    if (mass[i] != 0)
      std::cout << "MASS: " << i << " " << mass[i] << "\n" << std::endl;

  int count = 0;
   for (int i = 1; i <= n; i++) {
    for (int j = i; j <= n; j++) {
      if (map[i] >= 0 && map[j] >= 0) {
        setflag[i][j] = 1;
        if (i == j){
	  atom->set_mass(FLERR,i,mass[map[i]]);
        }
	  count++;
      }
      scale[i][j] = 1.0;
    }

   }
   int ilo,ihi,jlo,jhi;
   utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
   utils::bounds(FLERR,arg[1],1,atom->ntypes,jlo,jhi,error);

   double cut_one = cut_global;
   if (coeffflag && (narg == 3)) cut_one = utils::numeric(FLERR,arg[2],false,lmp);
   cut[0][0] = 20;
   cut[0][1] = 20;
   cut[1][1] = 20;
}
/* ----------------------------------------------------------------------
   init for one type pair i,j and corresponding j,i
------------------------------------------------------------------------- */

double PairMaise::init_one(int i, int j)
{
  if (setflag[i][j] == 0) {
    cut[i][j] = mix_distance(cut[i][i],cut[j][j]);
  }

  return cut[i][j];
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMaise::write_restart(FILE *fp)
{
  write_restart_settings(fp);

  int i,j;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      fwrite(&setflag[i][j],sizeof(int),1,fp);
      if (setflag[i][j]) {
        fwrite(&cut[i][j],sizeof(double),1,fp);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMaise::read_restart(FILE *fp)
{
  read_restart_settings(fp);
  allocate();

  int i,j;
  int me = comm->me;
  for (i = 1; i <= atom->ntypes; i++)
    for (j = i; j <= atom->ntypes; j++) {
      if (me == 0) utils::sfread(FLERR,&setflag[i][j],sizeof(int),1,fp,nullptr,error);
      MPI_Bcast(&setflag[i][j],1,MPI_INT,0,world);
      if (setflag[i][j]) {
        if (me == 0) {
          utils::sfread(FLERR,&cut[i][j],sizeof(double),1,fp,nullptr,error);
        }
        MPI_Bcast(&cut[i][j],1,MPI_DOUBLE,0,world);
      }
    }
}

/* ----------------------------------------------------------------------
   proc 0 writes to restart file
------------------------------------------------------------------------- */

void PairMaise::write_restart_settings(FILE *fp)
{
  fwrite(&cut_global,sizeof(double),1,fp);
  fwrite(&coeffflag,sizeof(int),1,fp);
}

/* ----------------------------------------------------------------------
   proc 0 reads from restart file, bcasts
------------------------------------------------------------------------- */

void PairMaise::read_restart_settings(FILE *fp)
{
  int me = comm->me;
  if (me == 0) {
    utils::sfread(FLERR,&cut_global,sizeof(double),1,fp,nullptr,error);
    utils::sfread(FLERR,&coeffflag,sizeof(int),1,fp,nullptr,error);
  }
  MPI_Bcast(&cut_global,1,MPI_DOUBLE,0,world);
  MPI_Bcast(&coeffflag,1,MPI_INT,0,world);
}

/* ----------------------------------------------------------------------
   proc 0 writes to data file
------------------------------------------------------------------------- */

void PairMaise::write_data(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    fprintf(fp,"%d\n",i);
}

/* ----------------------------------------------------------------------
   proc 0 writes all pairs to data file
------------------------------------------------------------------------- */

void PairMaise::write_data_all(FILE *fp)
{
  for (int i = 1; i <= atom->ntypes; i++)
    for (int j = i; j <= atom->ntypes; j++)
      fprintf(fp,"%d %d %g\n",i,j,cut[i][j]);
}

/* ---------------------------------------------------------------------- */

double PairMaise::single(int /*i*/, int /*j*/, int /* itype */, int /* jtype */,
                        double /* rsq */, double /*factor_coul*/,
                        double /* factor_lj */, double &fforce)
{
  fforce = 0.0;
  return 0.0;
}

