/* ----------------------------------------------------------------------
   Contributing author: SUN,Lixin (MIT) 
------------------------------------------------------------------------- */

#include	"abc.h"
#include	"atom.h"
#include	"comm.h"
#include	"compute.h"
#include	"displace_atoms.h"
#include	"domain.h"
#include	"error.h"
#include	"finish.h"
#include	"fix_abcstat.h"
#include	"fix_penalty.h"
#include	"fix_savconf.h"
#include	"fix.h"
#include	"force.h"
#include	"group.h"
#include	"math.h"
#include	"min.h"
#include	"minimize.h"
#include	"modify.h"
#include	"neighbor.h"
#include	"output.h"
#include	"run.h"
#include	"stdio.h"
#include	"stdlib.h"
#include	"string.h"
#include	"timer.h"
#include	"update.h"
#include	<time.h>

using namespace LAMMPS_NS;
#define avgTEMP 10
#define NAME 0
/* ---------------------------------------------------------------------- */

ABC::ABC(LAMMPS *lmp) : Pointers(lmp) {
  
  penaltyN = 0;
  stateN = 0;
  blockN = 0;
  min_stylearg = NULL;
  minarg = NULL;
  penaltyarg = NULL;
  blockarg = NULL;
  savconfarg = NULL;

  //default setting: 
  tradition_flag = 1; // traditional abc
  penaltyargN = -1;

  //default minimize style: sd
  if (min_stylearg == NULL){
    min_stylearg = new char*[1];
    min_stylearg[0] = new char[20];
    sprintf(min_stylearg[0],"sd");}
  error->warning(FLERR,"Change the minimize style to sd");
  update->create_minimize(1,min_stylearg);

}
/* ---------------------------------------------------------------------- */

void ABC::command(int narg, char **arg){

  if((comm->me ==0) && screen) fprintf(screen,"start abc\n");
  if((comm->me ==0) && logfile)  fprintf(logfile,"start abc\n");

  if (domain->box_exist == 0) error->all(FLERR,"ABC command before simulation box is defined");

  getarg(narg,arg); //read all the parameters
  int MaxBlockN = avgTEMP*MaxStateN+1;


  init(); //allocate memory and etc
 
  addState(); //add original state as state 0

  char bstring[100];
  bool lastMinima=true,visited=false;
  int state;
  int counter_attemp=1;
  double previous_pe=-1e20,drop;

  //record the number of penalty function and dump the configuration
  statistic->get_penalty(penaltyN);
  statistic->dump((char*)"penalty");
  statistic->dump((char*)"newconf");
  statistic ->print();
  bigint ntimestep;

  /* Begin the ABC sampling: at the beginning of the itenary, 
   * the activated group will be applied on a new group of
   * function, which is parameterized by the current coordinates,
   * and the penalty function exists until the end of abc sampling.
   *
   * ABCe: the system will jump back to the initial configuration 
   * after it finds a new basin.
   */
  while (( penaltyN < MaxPenaltyN ) && (stateN < MaxStateN) && (blockN < MaxBlockN) ){
    /* in ABCe, if the last state was a minimum, prepare a container
     * of 'block' penalty function for the pathway to be explored later
     * and also record such action in the output document */
    if (lastMinima){
      if(!tradition_flag) {
        //create the container
        addBlock();
        //count the times of attemping to find a new basin
        statistic->get_attemp(counter_attemp);
        statistic->breakpoint((char*)"new attemp !");}}
  
    //record the current potential energy;
    previous_pe = statistic->compute_vector(0);

    // record current coordinates and apply penalty function
    addPenalty(); 
    run(0);
    //optimization, the minimization style is defined ahead.
    minimize();

    statistic->dump((char*)"penalty");
    statistic ->print();
    drop = previous_pe - statistic->compute_vector(0);
    if ((drop > minima_drop) 
            && (statistic->compute_vector(2) < minima_tolerance)){
        //if current state is local minima
  
      lastMinima = true;
      counter_attemp++;

      if (!tradition_flag){
        blockList[blockN-1] ->enable(); //block the current pathway

        //check whether the state is visited before
        visited = false;
        for(state = 0;state < stateN;state++){
          if (stateList[state]->same()) { 
            visited = true;
            break;}
        }
        // if it's a new state, record it, type it and dump it.
        if (!visited) {
          addState();
          sprintf(bstring,"newstate %d %d", stateN,counter_attemp);
          statistic->breakpoint(bstring);
          statistic->dump((char*)"newconf");
        }
        // if the simulation is not finished
        // go back the most initial state
        if (( penaltyN < MaxPenaltyN ) && (stateN < MaxStateN)){
          if (!goback(ntimestep))  { 
            error->all(FLERR,"cannot goback in abc multiple");}
    	  else {
            sprintf(bstring,"goback");
            statistic->breakpoint(bstring);
            statistic->dump((char*)"penalty");
            statistic ->print(1);}
        }
      }
      else {
        //traditional ABC, record it, dummp it
        stateN++;
        sprintf(bstring,"newstate %d %d", stateN,counter_attemp);
        statistic->dump((char*)"newconf");
        statistic->breakpoint(bstring);
      }
    }//if local minima
    else   lastMinima = false;
  }//while 

  if((comm->me ==0) && screen) fprintf(screen,"finish abc\n");
  if((comm->me ==0) && logfile)  fprintf(logfile,"finish abc\n");
}

/* ---------------------------------------------------------------------- */

ABC::~ABC(){
  if ((comm->me == 0) && screen) fprintf(screen,"ABC::~ABC() %ld\n",update->ntimestep);
  char Temp[100];
  int i;
  error->warning(FLERR,"all the penalty functions and statistic fix commands are deleted");
  for(i = 0;i < penaltyN;i++) {
    sprintf(Temp,"penalty%d",i);
    modify->delete_fix(Temp);}
  for(i = 0;i < stateN;i++) {
    sprintf(Temp,"state%d",i);
    modify->delete_fix(Temp);}
  for(i = 0;i < blockN;i++) {
    sprintf(Temp,"block%d",i);
    modify->delete_fix(Temp);}
  modify->delete_fix("abcstat");
  
  error->warning(FLERR,"asdfadf");
  delete [] penaltyList;
  delete [] stateList;
  delete [] blockList;

  error->warning(FLERR,"asdfadf");
  delete [] minarg;
  delete [] penaltyarg;
  if (!tradition_flag) delete [] blockarg;
  delete [] savconfarg;
  delete [] plusEnergy;
  delete [] runarg;

  error->warning(FLERR,"asdfadf");
  update->restrict_output =0;
  error->warning(FLERR,"all dump and restart functions are recovered");
}

void ABC::getarg(int narg, char **arg){
  char str[248];
  if((comm->me ==0) && screen) fprintf(screen,"ABC::getarg\n");
  if((comm->me ==0) && logfile)  fprintf(logfile,"ABC::getarg\n");
/*
 group
 step, upper limit of steps
 state, upper limit of states
 min_style, sd
 min_arg 1-4 arg
 perturb yes
 penalty xyz W sigma2 cut_off perturb
 extend xyz W_block sigma2_block cut_off
 minima_cri, delta,drop ,single diff ,avg diff
*/
  
  if ((narg < 0)||(narg > 30)) error->all(FLERR,"Illegal ABC command");
  int iarg = 1;
  int i,n;
  while (iarg < narg) {
    if((comm->me ==0) && screen) fprintf(screen,"ABC::getarg while %s\n",arg[iarg]);
    if (strcmp(arg[iarg],"step") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal ABC command");
      MaxPenaltyN = atoi(arg[iarg+1]);
      if (MaxPenaltyN <= 0)  error->all(FLERR,"Illegal ABC command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"state") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal ABC command");
      MaxStateN = atoi(arg[iarg+1]);
      if (MaxStateN <= 0)  error->all(FLERR,"Illegal ABC command");
      iarg += 2;

    } else if (strcmp(arg[iarg],"min_style") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal ABC command");
      sprintf(min_stylearg[0],"%s",arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"min_arg") == 0) {
      if (iarg+5 > narg) error->all(FLERR,"Illegal ABC command");
      minarg = new char*[4];
      for (i=0;i<4;i++){
        n = strlen(arg[iarg+i+1]) + 1;
        minarg[i] = new char [n];
        strcpy(minarg[i], arg[iarg+1+i]);}
      iarg += 5;

    } else if (strcmp(arg[iarg],"perturb") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal ABC command");
      if (strcmp(arg[iarg+1],"yes") == 0) penaltyargN = 8;
      else if (strcmp(arg[iarg+1],"no") == 0 )penaltyargN = 7;
      else error->all(FLERR,"Illegal ABC command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"penalty") == 0) {
      if ( penaltyargN == -1 ) error->all(FLERR,"Illegal ABC command");
      if ((iarg+penaltyargN-2) > narg) error->all(FLERR,"Illegal ABC command");
      penaltyarg = new char*[penaltyargN];
      penaltyarg[0] = (char *) "fixPenalty";
      n = strlen(arg[0]) + 1;
      penaltyarg[1] = new char [n];
      strcpy(penaltyarg[1], arg[0]);
      penaltyarg[2] = (char *) "penalty";
      for (int i=3;i<(penaltyargN);i++){
        penaltyarg[i] = new char [strlen(arg[iarg+i-2]) + 1];
        strcpy(penaltyarg[i], arg[iarg+i-2]);}
      iarg += penaltyargN-2;
    } else if (strcmp(arg[iarg],"extend") == 0) {
      if ((iarg+5) > narg) error->all(FLERR,"Illegal ABC command");
      tradition_flag = 0;
      blockarg = new char*[7];
      blockarg[0] = (char *) "fixBlock";
      n = strlen(arg[0]) + 1;
      blockarg[1] = new char [n];
      strcpy(blockarg[1], arg[0]);
      blockarg[2] = (char *) "penalty";
      for (int i=3;i<7;i++){
        blockarg[i] = new char [strlen(arg[iarg+i-2]) + 1];
        strcpy(blockarg[i], arg[iarg+i-2]);}
      iarg += 5;
    } else if (strcmp(arg[iarg],"minima_cri") == 0) {
      if ((iarg+5) > narg) error->all(FLERR,"Illegal ABC command");
      savconfarg = new char*[5];
      savconfarg[0] = (char *) "savConf";
      savconfarg[1] =(char*) "all";
      savconfarg[2] =(char*) "savconf";
      savconfarg[3] = new char [strlen(arg[iarg+1]) + 1];
      strcpy(savconfarg[3], arg[iarg+1]);
      savconfarg[4] = new char [strlen(arg[iarg+2]) + 1];
      strcpy(savconfarg[4], arg[iarg+2]);
      minima_tolerance = atof(arg[iarg+3]);
      minima_drop = atof(arg[iarg+4]);
      iarg += 5;
    } else{ 
      sprintf(str,"Illegal ABC command %s ", arg[iarg]);
      error->all(FLERR,str);}
  }
  if ((minarg == NULL) || (penaltyarg == NULL) || ((!tradition_flag) && (blockarg == NULL)) || (savconfarg == NULL))
      error->all(FLERR,"Illegal ABC command");
  delete [] min_stylearg[0];
  delete [] min_stylearg;
}

void ABC::init(){

  if ((comm->me ==0)&& logfile) fprintf(logfile,"ABC::init\n");
  if ((comm->me ==0)&& screen) fprintf(screen,"ABC::init\n");

  /* memory for penalty function, block function 
   * and configurations of different stages */
  penaltyList = new class FixPenalty*[MaxPenaltyN+1];
  stateList = new class FixSavConf*[MaxStateN+1];
  blockList = new class FixPenalty*[avgTEMP*MaxStateN+1];

  // ban all dumps
  update->restrict_output =1;
  error->warning(FLERR,"all dump and restart are banned"); 

  // create fix for statistic and output
  if(modify->find_fix("abcstat") == -1){
    char **abcsarg = new char*[4];
    abcsarg[0] = (char*)"abcstat";
    abcsarg[1] = (char*)"all";
    abcsarg[2] = (char*)"abcstat";
    abcsarg[3] = (char*)"1";
    modify->add_fix(4,abcsarg);
    statistic = (class FixABCstat*) modify->fix[modify->nfix-1];
    delete [] abcsarg;}
  else statistic = (class FixABCstat*) modify->fix[modify->find_fix("abcstat")];

  // store frequently used arguments 
  plusEnergy = new char*[3];
  plusEnergy[0] =(char*) "Fix_Modify";
  plusEnergy[1] =(char*) "energy";
  plusEnergy[2] =(char*) "yes";

  runarg = new char*[1];
  runarg[0] =(char*) "0";

  run(0);
}

void ABC::addPenalty(){
  if ((comm->me ==0)&& logfile) fprintf(logfile,"addPenalty\n");
  if ((comm->me ==0)&& screen) fprintf(screen,"addPenalty\n");
  char TempName[64];
  sprintf(TempName,"penalty%d",penaltyN);
  penaltyarg[NAME]=TempName;
  plusEnergy[NAME]=TempName;
  fprintf(screen,"\n");
  modify->add_fix(penaltyargN,penaltyarg);
  penaltyList[penaltyN] = (class FixPenalty*)modify->fix[modify->nfix-1];
  modify->modify_fix(3,plusEnergy);

  penaltyN++;
  statistic->get_penalty(penaltyN);
}

void ABC::addState(){

  if ((comm->me ==0)&& logfile) fprintf(logfile,"add state\n");
  if ((comm->me ==0)&& screen) fprintf(screen,"add state \n");
  char TempName[64];
  sprintf(TempName,"state%d",stateN);

  savconfarg[NAME]=TempName;
  modify->add_fix(5,savconfarg);
  stateList[stateN] = (class FixSavConf*) modify->fix[modify->nfix-1];

  stateN++;
}

void ABC::addBlock(){  
  if ((comm->me ==0)&& logfile) fprintf(logfile,"add block\n");
  if ((comm->me ==0)&& screen) fprintf(screen,"add block \n");
  char TempName[64];
  sprintf(TempName,"block%d",blockN);

  blockarg[NAME] = TempName;
  plusEnergy[NAME] = TempName;

  modify->add_fix(7,blockarg);
  modify->modify_fix(3,plusEnergy);
  blockList[blockN] = (class FixPenalty*)modify->fix[modify->nfix-1];

  blockList[blockN]->disable();

  statistic->get_block(blockList[blockN]);
  blockN++;
}

bool ABC::goback(bigint ntimestep){
  if ((comm->me ==0)&& logfile) fprintf(logfile,"goback\n");
  if ((comm->me ==0)&& screen) fprintf(screen,"goback \n");
  if ((penaltyN >= MaxPenaltyN) || (stateN >= MaxStateN)){
    return false;}
  fprintf(logfile,"goback()\n");
  stateList[0]->goback();
  //update->reset_timestep(update->ntimestep+1);
  run(0);
  return true;
}

void ABC::run(int n){
  if ((comm->me ==0)&& logfile) fprintf(logfile,"run %d\n",n);
  if ((comm->me ==0)&& screen) fprintf(screen,"run %d\n",n);

  class Run * run;
  char temp[100];
  sprintf(temp,"%d",n);
  run = new class Run(lmp);
  runarg [0] = temp;
  run->command(1,runarg);
  delete run;
}

void ABC::minimize(){
  if ((comm->me ==0)&& logfile) fprintf(logfile,"minimize\n");
  if ((comm->me ==0)&& screen) fprintf(screen,"minimize\n");

  class Minimize * minimize;
  minimize = new class Minimize(lmp);
  minimize->command(4,minarg);    
  delete minimize;
  //min(4,minarg);
}
