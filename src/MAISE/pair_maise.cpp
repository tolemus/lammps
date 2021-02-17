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
#include <fstream>

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairMaise::PairMaise(LAMMPS *lmp) : Pair(lmp) {
  coeffflag=1;
  writedata=1;
  single_enable=1;
  respa_enable=1;
  scale = nullptr;

  j = 0;
  loflag = 0;
  
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

// default settings

  double **f = atom->f;
  int mN = atom->natoms;
  int NT = atom->ntypes;
  ANN mR;
  PRS mP;
  PRS mW[9];
  LNK mL;
  Cell mC;
  int mCODE;
  int mNM;
  int mND;
  int mNP;
  int mXT;
  int mATMN[mN];
  double mLAT[3][3];
  double mPOS[mN][3];
  double **x;
  double **v;

  double mH;
  double mFRC[mN][3];
  double mSTR[6];
  double mLATf[3*3];
  double mPOSf[mN*3];
  double mFRCf[mN*3];
  mCODE = 0;
  mNM   = 500;
  mND   = 0;
  mNP   = 4;
  mXT   = 1;
  mC.p  = 0;
  x = atom->x;
  v = atom->v;
  mP.EFS = 1;

  double mx;
  double my;
  double mz;

// set maise lattice
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

// set ND based on boundary properties
  for (i = 0; i < 3; i++)
    for (q = 0; q < 2; q++)
      if (domain->boundary[i][q] == 0){
	mND = 3;
	break;
      }
  mL.B = 0; 

// set atom types
  
  for (i=0;i<mN;i++)
    mATMN[i] = atom->type[i] - 1;

// flatten 2d arrays into 1d arrays

  for(i=0;i<3;i++)
    for(q=0;q<3;q++)
      mLATf[3*i+q] = mLAT[i][q];
  for(i=0;i<mN;i++)
    for(q=0;q<3;q++){
	mPOSf[3*i+q] = x[i][q]; 
    }
  
// call maise
/*  for (i=0;i<NT;i++)
    mC.spcz[i] = spcz[i];
*/
  mH = CALL_MAISE(&mR, &mP, mW, &mL, &mC, mCODE, mN, mNM, mND, mNP, mXT, mATMN, mLATf, mPOSf, mFRCf, mSTR);

// set energy

  eng_coul = mC.E;
  if (vflag_fdotr) virial_fdotr_compute();

// set lammps forces

  for(i=0;i<mN;i++)
    for(q=0;q<3;q++)
      f[i][q] = mFRCf[3*i+q];

  if(j<1 && loflag==1){
    LOUTCAR(&mC);
    j++;
  }

// print lammps forces three times
/*
  if(j<=2)
  for(i=0;i<mN;i++)
    for(q=0;q<3;q++)
      printf("%d %d forces, pos, type = % lf, % lf, % d\n",i, q, f[i][q], x[i][q], atom->type[i]);
  j++;
*/
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
  int i,q;  

  if (!allocated) allocate();

  if (narg > 2)
    error->all(FLERR,"Illegal pair_style command");

  cut_global = utils::numeric(FLERR,arg[0],false,lmp);
  
  loflag = utils::inumeric(FLERR,arg[1],false,lmp);

  if (loflag != 1 && loflag != 0)
    error->all(FLERR,"Illegal pair_style command");

  // reset cutoffs that have been explicitly set
  if (allocated) {
    for (i = 1; i <= atom->ntypes; i++)
      for (q = 1; q <= atom->ntypes; q++){
        cut[i][q] = cut_global;
	setflag[i][q] = 1;
      }
  }
}

/* ----------------------------------------------------------------------

   set coeffs for one or more type pairs
------------------------------------------------------------------------- */

void PairMaise::coeff(int narg, char **arg)
{
/*    int i,q;
    if ((narg == 0))
    error->all(FLERR,"Incorrect args for pair coefficients");

    if (narg != atom->ntypes)
      error->all(FLERR,"Incorrect args for pair coefficients");

    if (!allocated) allocate(); 

    for (i=0;i<narg;i++)
      spcz[i] = utils::numeric(FLERR,arg[i],false,lmp);



// set bounds

    int ilo,ihi,qlo,qhi;
    utils::bounds(FLERR,arg[0],1,atom->ntypes,ilo,ihi,error);
    utils::bounds(FLERR,arg[1],1,atom->ntypes,ilo,ihi,error);
    double cut_one = cut_global;

// required for setting flags setflag
  int count = 0;
  for (i=ilo;i<=ihi;i++)
    for (q=MAX(qlo,i);q<=qhi;q++){
      cut[i][q] = cut_one;
      setflag[i][q] = 1;
      count++;
    }
*/   
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

void PairMaise::LOUTCAR(Cell *mC)
{

  FILE *file;
  std::string l;
  int NLin = 0;
  int LNLin = 0;
  char s[200];

  printf("Hello LOUTCAR is working\n");

  if (file = fopen("OUTCAR", "r")) {
    system("cp OUTCAR OUTCAR.0");
    system("rm OUTCAR");
    CELL_OUT(mC);
    system("cp OUTCAR LOUTCAR");
    system("rm OUTCAR");
    system("cp OUTCAR.0 OUTCAR");
    system("rm OUTCAR.0");
    
    std::ifstream ifs("OUTCAR");
    std::ifstream Lifs("LOUTCAR");
//    while(std::getline(ifs, l))
//      NLin++;
//    while(std::getline(Lifs, l))
//      LNLin++;
    sprintf(s,"head -n %d OUTCAR | tail -n %d >> temp", mC->N+23, mC->N);
    system(s);
    sprintf(s,"head -n %d LOUTCAR | tail -n %d >> ltemp",mC->N+9, mC->N);
    system(s);
    system("diff temp ltemp >> comp");
    std::ifstream Fifs("comp");
    while(std::getline(Fifs, l))
      NLin++;
    system("rm temp ltemp");
    if (NLin==0){

      printf("FORCES MATCH\n");

    }
  } else {
   
    error->all(FLERR,"NO OUTCAR PLEASE CHANGE FLAG OR INSERT OUTCAR\n");
 
  }
}
