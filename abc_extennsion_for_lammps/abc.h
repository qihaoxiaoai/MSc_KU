/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under 
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: SUN,Lixin (MIT) 
------------------------------------------------------------------------- */

#ifdef COMMAND_CLASS

CommandStyle(abc,ABC)

#else

#ifndef LMP_ABC_H
#define LMP_ABC_H

#include "pointers.h"

namespace LAMMPS_NS {

class ABC : protected Pointers {
 public:
  ABC(class LAMMPS *);  
  ~ABC();
  void command(int, char **);
  void init();
  void getarg(int, char **);
  bool goback(bigint);
  void run(int);
  void minimize();
  void freeArg();
  //void min(int, char **);
 private:
  void addPenalty();
  void addState();
  void addBlock();

  char **penaltyarg, **plusEnergy; //parameter for savconf
  char **blockarg,  **savconfarg;
  char **runarg, **minarg, **min_stylearg;
  int penaltyargN;

  int MaxPenaltyN, MaxStateN;
  int penaltyN, stateN, blockN;

  double minima_tolerance,minima_drop;
  int tradition_flag;
  
  int originStep;
  
  class FixPenalty **penaltyList, **blockList;
  class FixABCstat *statistic;
  class FixSavConf **stateList;
};

}

#endif
#endif
