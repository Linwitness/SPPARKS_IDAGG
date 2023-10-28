/* -------------------------------------------------------
AppPotts_eng class source for abnormal grain growth
--
Read-Shockley implementation developed by Efrain Hernandez-Rivera (2017--2018)
US Army Research Laboratory
--
THIS SOFTWARE IS MADE AVAILABLE ON AN "AS IS" BASIS
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, NEITHER
EXPRESSED OR IMPLIED
------------------------------------------------------- */

#include "stdio.h"
#include "string.h"
#include "stdlib.h"
#include "domain.h"
#include "math.h"
#include "app.h"
#include "app_potts_eng.h"
#include "random_park.h"
#include "comm_lattice.h"
#include "error.h"
#include <fstream>
#include <iostream>
#include <type_traits>
#include <list>
using namespace SPPARKS_NS;

#define MY_PI 3.14159265358979323846 // pi
#define MY_2PI 6.28318530717958647692 // 2pi

template<typename T>
std::ostream& operator<< (std::ostream& out, const std::vector<T>& v) {
    out << "{";
    size_t last = v.size() - 1;
    for(size_t i = 0; i < v.size(); ++i) {
        out << v[i];
        if (i != last)
            out << ", ";
    }
    out << "}";
    return out;
}

/* ---------------------------------------------------- */

// I think this part is used for initialization
AppPotts_eng::AppPotts_eng(SPPARKS *spk, int narg, char **arg) :
AppPotts(spk,narg,arg)
{
  ninteger = 1;
  //1 double array per Euler angle
  ndouble = 3;
  // delpropensity = 6;
  // delevent = 0;
  // add the extra arrays
  recreate_arrays();
  // only error check for this class, not derived classes
  if (strcmp(arg[0],"potts/eng") == 0 && narg < 2)
  error->all(FLERR,"Illegal app_style command");
  //cutoff misorientation angle
  thetam=25.0/180.0*MY_PI;
  //interaction (interfacial) energy
  Jij=1.0;
  //Mobility parameters
  // nmob=4.0; bmob=5.0;
  nspins=atoi(arg[1]);
  //Symmetry operator
  Osym=24;
  if (narg == 3)
  Osym=atoi(arg[2]);
  //Inclination default parameters
  maxIncEnergy = 1; numBins = 5;
  triple_energy = "ave";
  // smoothing algorithm iteration times
  interval = 3;
  neighbor_length = interval+1;
  clip = 0;
  storage_flag = 0;
  reference_axis[0] = 1;
}

/* -------------------------------------------------------
Destructor
------------------------------------------------------- */
AppPotts_eng::~AppPotts_eng()
{
  //free up memory from quaternion symmetry operator
  for (int i = 0; i<Osym; i++)
  delete[] symquat[i];
  delete[] symquat;
}

/* -------------------------------------------------------
Initialize before each run
check validity of site values
------------------------------------------------------- */
void AppPotts_eng::init_app()
{

  // the coordination of site i: xyz[i][0] - x; xyz[i][1] - y; xyz[i][2] - z;
  // xyz = app->xyz;

  int flag = 0;
  //Check angles are within corresponding range
  for (int i = 0; i < nlocal; i++) {
    if (phi1[i] < 0 || phi1[i] >= MY_2PI){
      fprintf(screen, "phi1 %d\n",i);
      flag = 1;
    }
    if (phi2[i] < 0 || phi2[i] >= MY_2PI) {
      fprintf(screen, "phi2 %d\n",i);
      flag = 1;
    }
    if (Phi[i] < 0 || Phi[i] >= MY_PI) {
      fprintf(screen, "Phi %d\n",i);
      flag = 1;
    }
    if (spin[i] < 1 || spin[i] > nspins+1) {
      fprintf(screen, "Spin %d\n",i);
      flag = 1;
    }
  }

  //Initialize symmetry operator as quaternion vectors
  //Osym = 24 (cubic), 12 (hexagonal)
  symmat(&symquat);

  comm->all();   // what's your function?

  if (logfile)
    fprintf(logfile," Pairs misorientation map created\n");
  if (screen && me==0)
    fprintf(screen," Pairs misorientation map created\n");

  int flagall;
  MPI_Allreduce(&flag,&flagall,1,MPI_INT,MPI_SUM,world);
  if (flagall)
    error->all(FLERR,"One or more sites have invalid values");

  // initial the parameters and matrix we need
  nx = ceil(domain->boxxhi);
  ny = ceil(domain->boxyhi);
  nz = ceil(domain->boxzhi);
  dimension = domain->dimension;

  if (dimension == 2) {
    linear_vector_matrix_i.assign(2*interval+3, std::vector<double>(2*interval+3, 0));
    linear_vector_matrix_j.assign(2*interval+3, std::vector<double>(2*interval+3, 0));
    // create a 2d matrix to store the site ID;
    structure_matrix_global.assign(ny ,std::vector<int>(nx));
    for (int smg_i=0; smg_i < ny; smg_i++) {
      for (int smg_j=0; smg_j < nx; smg_j++) {
        structure_matrix_global[smg_i][smg_j] = smg_i*nx + smg_j;
      }
    }

    // New neareast neighbor functions
    numneigh_near = 8;
    neighbor_near.assign(numneigh_near, 0);
  }
  else if (dimension == 3) {
    linear_vector_matrix_i_3d.assign(2*interval+3, std::vector<std::vector<double>>(2*interval+3, std::vector<double>(2*interval+3, 0)));
    linear_vector_matrix_j_3d.assign(2*interval+3, std::vector<std::vector<double>>(2*interval+3, std::vector<double>(2*interval+3, 0)));
    linear_vector_matrix_k_3d.assign(2*interval+3, std::vector<std::vector<double>>(2*interval+3, std::vector<double>(2*interval+3, 0)));
    // create a 3d matrix to store the site ID; &&&
    structure_matrix_global_3d.assign(nz ,std::vector<std::vector<int>>(ny, std::vector<int>(nx)));
    for (int smg_i=0; smg_i < nz; smg_i++) {
      for (int smg_j=0; smg_j < ny; smg_j++) {
        for (int smg_k=0; smg_k < nx; smg_k++) {
          structure_matrix_global_3d[smg_i][smg_j][smg_k] = smg_i*nx*ny + smg_j*nx + smg_k;
        }
      }
    }

    // New neareast neighbor functions
    numneigh_near = 26;
    neighbor_near.assign(numneigh_near, 0);
  }

  delete [] sites; // include all object(like neighbors) the i site can flip
  delete [] unique;
  sites = new int[1 + numneigh_near]; // !!! Be careful !!!
  unique = new int[1 + numneigh_near]; // !!! Be careful !!!

  dt_sweep = 1.0/numneigh_near; // !!! Be careful !!!

  // Inclination vector in lab and self frame assignment
  vector_i_lab_frame.assign(nglobal, std::vector<double>(3,0));
  vector_nei_lab_frame.assign(nglobal, std::vector<double>(3,0));
  average_normal_i_lab_frame.assign(nglobal, std::vector<double>(3,0));
  average_normal_nei_lab_frame.assign(nglobal, std::vector<double>(3,0));
  average_normal_lab_frame_tri.assign(nglobal, std::vector<double>(3,0));
  inclination_i_frame.assign(nglobal, std::vector<double>(3,0));
  inclination_nei_frame.assign(nglobal, std::vector<double>(3,0));
  inclination_quaternion_i_frame.assign(nglobal, std::vector<double>(4,0));
  inclination_quaternion_nei_frame.assign(nglobal, std::vector<double>(4,0));
  inclination_quaternion_symmetry_i_frame.assign(nglobal, std::vector<double>(4,0));
  inclination_quaternion_symmetry_nei_frame.assign(nglobal, std::vector<double>(4,0));

  // storage for inclination
  inc_storage.assign(nglobal,std::unordered_map<int, std::vector<double>>());
  //Read in misorientations from a file if present, else calculate misorientations and write file
  std::ifstream misoFile("MisoEnergy.txt");
  if (misoFile.peek() == std::ifstream::traits_type::eof()) {
    // //Save misorientation list to file

  }
  //Read in misorientations and replace 0 grain boundary energy entries
  else {
    read_misorientation();
  }
  misoFile.close();

  std::ifstream vectorMatrixFile("vectorMatrix/vectorMatrix_" +std::to_string(interval)+ ".txt");
  if (vectorMatrixFile.peek() == std::ifstream::traits_type::eof()) {
    // //Save vectorMatrixFile list to file

  }
  //Read in vectorMatrixFile
  else {
    read_linear_vector_matrix();
  }
  vectorMatrixFile.close();

  //If there isnt an inclination lookup table file make one, else read it in
  incBins.assign(numBins+1,0.0); //Store inclination values that mark bin thresholds
  for (double i = 0; i<numBins+1; i++){
    incBins[i] = -1.0+(2.0/numBins*i)+0.01;
  }
  incTable.assign(numBins,std::vector<double>(numBins)); //Set size of inclination table
  std::ifstream incFile("incTable.txt");
  if (incFile.peek() == std::ifstream::traits_type::eof()) {
    init_inclination_table();
  }
  else{
    read_inclination_table();
  }
  incFile.close();

  //If there isnt an energy file make one, else read it in
  std::ifstream energyFile("Energy.txt");
  if (energyFile.peek() == std::ifstream::traits_type::eof()) {
    init_energy();
  }
  else{
    read_energy();
  }
  energyFile.close();

  // create a output file to store the vector information
  std::ifstream vectorFile("Vector.txt");
  // std::ofstream vectorFile("Vector.txt", std::ios::out | std::ios::app);
  // vectorFile << "i spin v0 v1 xyz0 xyz1\n";
  vectorFile.close();

}

/* -------------------------------------------------------
Read in misorientations and replace 0 grain boundary energy entries
------------------------------------------------------- */
void AppPotts_eng::read_misorientation()
{
  double num;
  std::ifstream misoFile("MisoEnergy.txt");
  for (int j = 2; j<nspins+1; j++){
    for (int i = 1; i<j; i++){
      std::pair <int, int> spins = std::make_pair(i,j);
      misoFile>>num;
      if (num==0){
        num=0.0001;
      }
      misos[spins]=num;
    }
  }
}
/* -------------------------------------------------------
Initialize default (all zeros) energy file
------------------------------------------------------- */
void AppPotts_eng::init_energy()
{
  double num;
  std::ofstream energyFile("Energy.txt");
  for (int j = 2; j<nspins+1; j++){
    for (int i = 1; i<j; i++){
      energyFile << 0 << "\n";
      std::pair <int, int> spins = std::make_pair(i,j);
      energy[spins]=0;
    }
  }
}
/* -------------------------------------------------------
Read in energy file
------------------------------------------------------- */
void AppPotts_eng::read_energy()
{
  double num;
  std::ifstream energyFile("Energy.txt");
  for (int j = 2; j<nspins+1; j++){
    for (int i = 1; i<j; i++){
      std::pair <int, int> spins = std::make_pair(i,j);
      energyFile>>num;
      energy[spins]=num;
    }
  }
}

/* -------------------------------------------------------
Create an inclination look up table based on user provided parameters
------------------------------------------------------- */
void AppPotts_eng::init_inclination_table()
{
  double xEn, yEn, incVal;
  std::ofstream incFile("incTable.txt");
  for (int i = 0; i<numBins; i++){ //X dimension bins
    for (int j = 0; j<numBins; j++){ //Y dimension bins
      xEn = fabs((incBins[i]+incBins[i+1])/2);
      yEn = fabs((incBins[j]+incBins[j+1])/2);
      incVal = maxIncEnergy/2*(xEn+yEn);
      if (incVal < 0.0001){
        incVal = 0;
      }
      incTable[i][j] = incVal;
      incFile << incTable[i][j] << "\n";
    }
  }
}
/* -------------------------------------------------------
read in inclination look up table
------------------------------------------------------- */
void AppPotts_eng::read_inclination_table()
{
  double num;
  std::ifstream incFile("incTable.txt");
  for (int i = 0; i<numBins; i++){ //X dimension bins
    for (int j = 0; j<numBins; j++){ //Y dimension bins
      incFile >> incTable[i][j];
    }
  }
}

void AppPotts_eng::read_linear_vector_matrix() {
  int ma_len = 2*interval+3;
  std::ifstream vectorMatrixFile("vectorMatrix/vectorMatrix_" +std::to_string(interval)+ ".txt");
  if (dimension==2) {
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        vectorMatrixFile >> linear_vector_matrix_i[i][j];
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        vectorMatrixFile >> linear_vector_matrix_j[i][j];
  }
  else if (dimension==3) {
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        for (int k=0; k<ma_len; k++)
          vectorMatrixFile >> linear_vector_matrix_i_3d[i][j][k];
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        for (int k=0; k<ma_len; k++)
          vectorMatrixFile >> linear_vector_matrix_j_3d[i][j][k];
    for (int i=0; i<ma_len; i++)
      for (int j=0; j<ma_len; j++)
        for (int k=0; k<ma_len; k++)
          vectorMatrixFile >> linear_vector_matrix_k_3d[i][j][k];
  }
}

/* -------------------------------------------------------
Set site value ptrs each time iarray/darray are
reallocated
------------------------------------------------------- */
void AppPotts_eng::grow_app()
{
  // set pointers
  // to define these, use command
  // create_sites box iN and set iN
  spin = iarray[0];
  phi1 = darray[0];
  Phi = darray[1];
  phi2 = darray[2];
}

/* -------------------------------------------------------
User defined optional parameters
------------------------------------------------------- */
void AppPotts_eng::input_app(char *command, int narg, char **arg)
{
  //Redefine mobility parameters (n,b)
  if (strcmp(command,"mobility") == 0) {

  }
  //Cutoff angle for Read-Shockley
  else if (strcmp(command,"cutoff") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal cutoff angle command\n");
    thetam=fabs(atof(arg[0]))/180.0*MY_PI;
    if (thetam>MY_2PI)
      error->all(FLERR,"Cutoff angle must be defined in "
        "terms of degrees (0,360)\n");

    if (logfile)
      fprintf(logfile," Low-to-high angle cutoff reset "
        "to %s deg\n",arg[0]);
    if (screen && me==0)
      fprintf(screen," Low-to-high angle cutoff reset "
        "to %s deg\n",arg[0]);
  }
  //Potts interfacial energy scaler
  else if (strcmp(command,"energy_scaling") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal scaling energy command\n");
    Jij=atof(arg[0]);
    if (Jij<0)
      error->all(FLERR,"Illegal energy value (>0)\n");
    if (logfile)
      fprintf(logfile," PMC energy scaling by %g.\n",Jij);
    if (screen && me==0)
      fprintf(screen," PMC energy scaling by %g.\n",Jij);
  }
  else if (strcmp(command,"incParams") == 0) {
    if (narg != 4){
      error->all(FLERR,"Illegal incParams flag: requires "
        "two arguments, type double maximum inclination energy and type int number of discretized energy bins\n "
        "(e.g. incParams 1 5)\n");
    }
    else{
      maxIncEnergy = atof(arg[0]);
      numBins = atof(arg[1]);
      E_delta = atof(arg[2]);
      E_m = atof(arg[3]);
    }
  }
  else if (strcmp(command,"interval") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal smoothing algorithm iteration command\n");
    else {
      interval=atof(arg[0]);
      neighbor_length = interval+1;
      smthAlgo = arg[1];
    }
    if (logfile)
      fprintf(logfile," We use smoothing algorithm %s with %d iteration.\n",smthAlgo.c_str(), interval);
    if (screen && me==0)
      fprintf(screen," We use smoothing algorithm %s with %d iteration.\n",smthAlgo.c_str(), interval);
  }
  else if (strcmp(command,"clip") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal margin clip command\n");
    clip=atof(arg[0]);
    if (clip>interval)
      error->all(FLERR,"Illegal clipped margin size should be smaller than interval");
    if (logfile)
      fprintf(logfile," The margin of length %d are clipped to increase efficiency.\n",clip);
    if (screen && me==0)
      fprintf(screen," The margin of length %d are clipped to increase efficiency.\n",clip);
  }
  else if (strcmp(command,"storage_flag") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal margin clip command\n");
    storage_flag=atof(arg[0]);
    if (logfile)
      fprintf(logfile," The inclination storage is opened.\n");
    if (screen && me==0)
      fprintf(screen," The inclination storage is opened.\n");
  }
  else if (strcmp(command,"reference_axis") == 0) {
    if (narg<3)
      error->all(FLERR,"Illegal reference axis\n");
    reference_axis[0]=atof(arg[0]);
    reference_axis[1]=atof(arg[1]);
    reference_axis[2]=atof(arg[2]);
    if (logfile)
      fprintf(logfile," The reference axis is {%f, %f, %f}.\n", reference_axis[0], reference_axis[1], reference_axis[2]);
    if (screen && me==0)
      fprintf(screen," The reference axis is {%f, %f, %f}.\n", reference_axis[0], reference_axis[1], reference_axis[2]);
  }
  else if (strcmp(command,"TJ_energy_type") == 0) {
    if (narg<1)
      error->all(FLERR,"Illegal triple junction energy type\n");
    triple_energy=arg[0];
    if (logfile)
      fprintf(logfile," The constant triple energy is %s.\n", triple_energy.c_str());
    if (screen && me==0)
      fprintf(screen," The constant triple energy is %s.\n", triple_energy.c_str());
  }
  else
    error->all(FLERR,"Input command not recognized by app\n");
}

void AppPotts_eng::init_vector(int i, int oldstate, double energy0, double e0_add, double energy1, double e1_add)
{
  std::ofstream vectorFile("Vector.txt", std::ios::out | std::ios::app);
  // vectorFile << i << " " << spin[i] << " " << vectori[0] << " " << vectori[1] << " " << vectori[2] << " " << xyz[i][0] << " " << xyz[i][1] <<" "<< xyz[i][2] << "\n";
  vectorFile << i << " " << oldstate << " vector: " << \
  vector_i_lab_frame[i][0] << " " << vector_i_lab_frame[i][1] << " " << vector_i_lab_frame[i][2] << " inc: " << \
  inclination_i_frame[i][0] << " " << inclination_i_frame[i][1] << " " << inclination_i_frame[i][2] << " coordinate: " << \
  xyz[i][0] << " " << xyz[i][1] <<" "<< xyz[i][2] << " energy: " << \
  energy0 << " add_init: " << e0_add << "\n";

  vectorFile << i << " " << spin[i] << " vector: " << \
  " inc: " << \
  inclination_nei_frame[i][0] << " " << inclination_nei_frame[i][1] << " " << inclination_nei_frame[i][2] << " coordinate: " << \
  xyz[i][0] << " " << xyz[i][1] <<" "<< xyz[i][2] << " energy: " << \
  energy1 << " add_init: " << e1_add << "\n";
  vectorFile.close();
  // vector_nei_lab_frame[i][0] << " " << vector_nei_lab_frame[i][1] << " " << vector_nei_lab_frame[i][2] <<
}

void AppPotts_eng::get_nearest_neighbor(int i, int neighbor_size) {

  if (dimension == 2) {
    neighbor_near[0] = neighbor[i][2*neighbor_size*(neighbor_size+2)];
    neighbor_near[1] = neighbor[i][2*neighbor_size*(neighbor_size+2)+1];
    neighbor_near[2] = neighbor[i][2*neighbor_size*(neighbor_size+2)+2];
    neighbor_near[3] = neighbor[i][2*neighbor_size*neighbor_size+6*neighbor_size+3];
    neighbor_near[4] = neighbor[i][2*neighbor_size*neighbor_size+6*neighbor_size+4];
    neighbor_near[5] = neighbor[i][2*neighbor_size*neighbor_size+8*neighbor_size+5];
    neighbor_near[6] = neighbor[i][2*neighbor_size*neighbor_size+8*neighbor_size+6];
    neighbor_near[7] = neighbor[i][2*neighbor_size*neighbor_size+8*neighbor_size+7];
  }
  else if (dimension == 3) {
    neighbor_near[0] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+13*neighbor_size];
    neighbor_near[1] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+13*neighbor_size+1];
    neighbor_near[2] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+13*neighbor_size+2];
    neighbor_near[3] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+15*neighbor_size+3];
    neighbor_near[4] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+15*neighbor_size+4];
    neighbor_near[5] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+15*neighbor_size+5];
    neighbor_near[6] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+17*neighbor_size+6];
    neighbor_near[7] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+17*neighbor_size+7];
    neighbor_near[8] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+14*neighbor_size*neighbor_size+17*neighbor_size+8];

    neighbor_near[9] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+25*neighbor_size+9];
    neighbor_near[10] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+25*neighbor_size+10];
    neighbor_near[11] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+25*neighbor_size+11];
    neighbor_near[12] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+27*neighbor_size+12];
    neighbor_near[13] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+27*neighbor_size+13];
    neighbor_near[14] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+29*neighbor_size+14];
    neighbor_near[15] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+29*neighbor_size+15];
    neighbor_near[16] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+18*neighbor_size*neighbor_size+29*neighbor_size+16];

    neighbor_near[17] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+37*neighbor_size+17];
    neighbor_near[18] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+37*neighbor_size+18];
    neighbor_near[19] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+37*neighbor_size+19];
    neighbor_near[20] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+39*neighbor_size+20];
    neighbor_near[21] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+39*neighbor_size+21];
    neighbor_near[22] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+39*neighbor_size+22];
    neighbor_near[23] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+41*neighbor_size+23];
    neighbor_near[24] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+41*neighbor_size+24];
    neighbor_near[25] = neighbor[i][4*neighbor_size*neighbor_size*neighbor_size+22*neighbor_size*neighbor_size+41*neighbor_size+25];
  }


}

/* -------------------------------------------------------
Compute Hamiltonian of site
------------------------------------------------------- */
double AppPotts_eng::site_energy(int i)
{
  if (spin[i] > nspins) return 0.0;
  int nei, nei_num = 0;
  double eng = 0.0, thetar;
  nei_set.clear();
  get_nearest_neighbor(i, interval);

  // test
  // std::ofstream vectorFile("Vector.txt", std::ios::out | std::ios::app);
  // if (xyz[i][0] == 256) {
  //   vectorFile << i << " " << spin[i] << " xyz: " << \
  //   xyz[i][0] << " " << xyz[i][1] <<" "<< xyz[i][2] << "\n" << "numneigh_near: " << \
  //   numneigh_near << " nei: " << neighbor_near << "\n";
  // }
  // std::cout << "i: " << i;

  // If the site is not on the boundary
  for (int j = 0; j < numneigh_near; j++) {
    nei=neighbor_near[j];
    // if (xyz[i][0]==7 && xyz[i][1]==7 && xyz[i][2]==7) std::cout<<" |"<<nei<<"| ";
    // Get neighboring grain ID
    if (spin[i] == spin[nei]) continue;
    nei_set.insert(spin[nei]);
    // Get num of neighboring site
    nei_num += 1;
  }
  // std::cout << " NotB0√ ";
  if (!nei_set.empty()) {
    // std::cout << " num_nei: " << nei_set.size();
    inclination_calculation(i, i, vector_i_lab_frame[i], i);
  }
  // std::cout << " NotB1("<< nei_num <<")√ ";

  // energy
  // Boundary energy
  if (nei_set.size() == 1)
    for (int j = 0; j < numneigh_near; j++) {
      double tmp_eng=0;
      nei=neighbor_near[j];

      if (spin[i] == spin[nei]) continue;
      // std::cout << " BE0("<< xyz[nei][0] << " "<< xyz[nei][1] << " "<< xyz[nei][2] <<")√";

      // add the inclination energy into total site energy
      // tmp_eng += energy_ave(i, nei);
      tmp_eng += energy_ave_split(i, nei);
      eng+=tmp_eng;
      // std::cout << " BE1√ ";
      // test
      // std::cout << xyz[nei][0] <<" "<<xyz[nei][1] <<" "<<xyz[nei][2] <<\
      // "\n   vector_i: " << vector_i_lab_frame[i] << " nei: " << vector_i_lab_frame[nei] << " ave: " << average_normal_i_lab_frame[nei] <<\
      // "\n   vector_nei: " << vector_nei_lab_frame[nei] << " nei: " << vector_nei_lab_frame[i] << " ave: " << average_normal_nei_lab_frame[nei] <<\
      // "\n   tmp_energy: " << tmp_eng << " numneighbor: " << numneigh[i] << numneigh[nei] << "\n";

    }
  else if (nei_set.size() > 1) {
    if (triple_energy == "min") {
      // Triple depend on lowest energy boundary site
      double min_eng = 8;
      for (int j = 0; j < numneigh_near; j++) {
        double tmp_eng = 0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // add the inclination energy into total site energy
        tmp_eng += energy_ave_split(i, nei);
        if (min_eng > tmp_eng) min_eng = tmp_eng;
      }
      eng += nei_num * min_eng;
    }
    else if (triple_energy == "max") {
      // Triple depend on highest energy boundary site
      double max_eng = 0;
      for (int j = 0; j < numneigh_near; j++) {
        double tmp_eng = 0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // add the inclination energy into total site energy
        tmp_eng += energy_ave_split(i, nei);
        if (max_eng < tmp_eng) max_eng = tmp_eng;
      }
      eng += nei_num * max_eng;
    }
    else if (triple_energy == "ave") {
      // Triple depend on average energy boundary site
      double sum_eng = 0;
      std::vector<double> other_grain_eng(nei_set.size(), 0); // store energy for other grains except for spin[i]
      std::vector<int> other_grain_num(nei_set.size(), 0); // store num of sites for other grains except for spin[i]
      for (int j = 0; j < numneigh_near; j++) {
        double tmp_eng = 0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // // Get the inclination and eng for nei
        // inclination_calculation(nei, nei, vector_nei_lab_frame[nei], i);
        // vector_to_quternion(vector_nei_lab_frame[nei], reference_axis, inclination_quaternion_nei_frame[nei]);
        // Store the energy and num of sites into corresponding grain sequence
        std::set<int>::iterator iter;
        int k = 0;
        for (iter = nei_set.begin(); iter != nei_set.end(); ++iter) {
          if (spin[nei] == *iter) {
            other_grain_eng[k] += energy_ave_split(i, nei); // compute_inclination_energy(inclination_quaternion_nei_frame[nei]);
            other_grain_num[k] += 1;
          }
          k++;
        }
      }
      // Add average energy from other grains except for spin[i]
      for (int k = 0; k < other_grain_eng.size(); k++) sum_eng += other_grain_eng[k] / other_grain_num[k];
      // // Add average energy from spin[i] grain
      // vector_to_quternion(vector_i_lab_frame[i], reference_axis, inclination_quaternion_i_frame[i]);
      // sum_eng += compute_inclination_energy(inclination_quaternion_i_frame[i]);
      eng += nei_num * sum_eng / (nei_set.size() );
      // std::cout << " TE("<< xyz[i][0] << " "<< xyz[i][1] << " "<< xyz[i][2] <<")√ ";
    }
    else if (triple_energy == "con") {
      // Triple depend on average boundary site (constant surface)
      double sum_eng = 0;
      std::vector<double> other_grain_eng(nei_set.size(), 0); // store energy for other grains except for spin[i]
      std::vector<int> other_grain_num(nei_set.size(), 0); // store num of sites for other grains except for spin[i]
      for (int j = 0; j < numneigh_near; j++) {
        double tmp_eng = 0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // Get the inclination and eng for nei
        inclination_calculation(nei, nei, vector_nei_lab_frame[nei], i);
        vector_to_quternion(vector_nei_lab_frame[nei], reference_axis, inclination_quaternion_nei_frame[nei]);
        // Store the energy and num of sites into corresponding grain sequence
        std::set<int>::iterator iter;
        int k = 0;
        for (iter = nei_set.begin(); iter != nei_set.end(); ++iter) {
          if (spin[nei] == *iter) {
            other_grain_eng[k] += compute_inclination_energy(inclination_quaternion_nei_frame[nei]);
            other_grain_num[k] += 1;
          }
          k++;
        }
      }
      // Add average energy from other grains except for spin[i]
      for (int k = 0; k < other_grain_eng.size(); k++) sum_eng += other_grain_eng[k] / other_grain_num[k];
      // Add average energy from spin[i] grain
      vector_to_quternion(vector_i_lab_frame[i], reference_axis, inclination_quaternion_i_frame[i]);
      sum_eng += compute_inclination_energy(inclination_quaternion_i_frame[i]);
      eng += 4.0 * sum_eng / (nei_set.size() + 1);
    }
    else if (triple_energy == "sum") {
      // Triple depend on average boundary site (higher energy)
      double sum_eng = 0;
      std::vector<double> other_grain_eng(nei_set.size(), 0); // store energy for other grains except for spin[i]
      std::vector<int> other_grain_num(nei_set.size(), 0); // store num of sites for other grains except for spin[i]
      for (int j = 0; j < numneigh_near; j++) {
        double tmp_eng = 0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // Get the inclination and eng for nei
        inclination_calculation(nei, nei, vector_nei_lab_frame[nei], i);
        vector_to_quternion(vector_nei_lab_frame[nei], reference_axis, inclination_quaternion_nei_frame[nei]);
        // Store the energy and num of sites into corresponding grain sequence
        std::set<int>::iterator iter;
        int k = 0;
        for (iter = nei_set.begin(); iter != nei_set.end(); ++iter) {
          if (spin[nei] == *iter) {
            other_grain_eng[k] += compute_inclination_energy(inclination_quaternion_nei_frame[nei]);
            other_grain_num[k] += 1;
          }
          k++;
        }
      }
      // Add average energy from other grains except for spin[i]
      for (int k = 0; k < other_grain_eng.size(); k++) sum_eng += other_grain_eng[k] / other_grain_num[k];
      // Add average energy from spin[i] grain
      vector_to_quternion(vector_i_lab_frame[i], reference_axis, inclination_quaternion_i_frame[i]);
      sum_eng += compute_inclination_energy(inclination_quaternion_i_frame[i]);
      eng += nei_num * sum_eng;
    }
    else if (triple_energy == "consMin") {
      // Triple depend on constant boundary site (minimal energy)
      double global_min_eng = 0;
      global_min_eng = 1 - E_delta;
      eng += nei_num * global_min_eng;
    }
    else if (triple_energy == "consMax") {
      // Triple depend on constant boundary site (minimal energy)
      double global_max_eng = 0;
      global_max_eng = 1 + E_delta;
      eng += nei_num * global_max_eng;
    }
    else if (triple_energy == "consTest") {
      // Triple depend on constant boundary site (minimal energy)
      double global_constest_eng = 0.093;
      eng += nei_num * global_constest_eng;
    }
    else if (triple_energy == "old") {
      // Triple depend on constant boundary site (minimal energy)

      for (int j = 0; j < numneigh_near; j++) {
        double old_TJ_eng = 0;
        nei = neighbor_near[j];
        if (spin[i] == spin[nei]) continue;
        // add the inclination energy into total site energy
        old_TJ_eng = energy_ave_split(i, nei);
        eng += old_TJ_eng;
      }

    }

  }

  // if (nei_num != 0) current_gamma = eng / nei_num;
  // else
  current_gamma = 1.0;

  // std::cout << " eng: " << eng <<std::endl;
  return Jij*eng;
}

double AppPotts_eng::energy_ave(int i, int nei) {
  inclination_calculation(i, nei, vector_i_lab_frame[nei], i);
  vector_average(vector_i_lab_frame[i], vector_i_lab_frame[nei], average_normal_i_lab_frame[nei]); // Get the average normals for sites in two sides with lab frame

  // Calculate the normal vector of nei's grain in lab frame
  inclination_calculation(nei, nei, vector_nei_lab_frame[nei], i); // ...
  inclination_calculation(nei, i, vector_nei_lab_frame[i], i); // ...
  vector_average(vector_nei_lab_frame[nei], vector_nei_lab_frame[i], average_normal_nei_lab_frame[i]); // ...

  vector_average(average_normal_i_lab_frame[nei], average_normal_nei_lab_frame[i], average_normal_lab_frame_tri[nei]);

  // calculate the roatation of quaternion version
  vector_to_quternion(average_normal_lab_frame_tri[nei], reference_axis, inclination_quaternion_i_frame[nei]);

  // add the inclination energy into total site energy
  return compute_inclination_energy(inclination_quaternion_i_frame[nei]);
}

double AppPotts_eng::energy_ave_split(int i, int nei) {
  inclination_calculation(i, nei, vector_i_lab_frame[nei], i);
  vector_average(vector_i_lab_frame[i], vector_i_lab_frame[nei], average_normal_i_lab_frame[nei]); // Get the average normals for sites in two sides with lab frame

  vector_to_quternion(average_normal_i_lab_frame[nei], reference_axis, inclination_quaternion_i_frame[nei]);

  // Calculate the normal vector of nei's grain in lab frame
  inclination_calculation(nei, nei, vector_nei_lab_frame[nei], i); // ...
  inclination_calculation(nei, i, vector_nei_lab_frame[i], i); // ...
  vector_average(vector_nei_lab_frame[nei], vector_nei_lab_frame[i], average_normal_nei_lab_frame[nei]); // ...

  // vector_average(average_normal_i_lab_frame[nei], average_normal_nei_lab_frame[i], average_normal_lab_frame_tri[nei]);

  // calculate the roatation of quaternion version
  vector_to_quternion(average_normal_nei_lab_frame[nei], reference_axis, inclination_quaternion_nei_frame[nei]);

  // add the inclination energy into total site energy
  return 0.5 * (compute_inclination_energy(inclination_quaternion_i_frame[nei]) + compute_inclination_energy(inclination_quaternion_nei_frame[nei]));
}

double AppPotts_eng::energy_diff(int i, int nei) {
  inclination_calculation(i, nei, vector_i_lab_frame[nei], i);
  vector_average(vector_i_lab_frame[i], vector_i_lab_frame[nei], average_normal_i_lab_frame[nei]); // Get the average normals for sites in two sides with lab frame

  // Calculate the normal vector of nei's grain in lab frame
  inclination_calculation(nei, nei, vector_nei_lab_frame[nei], i); // ...
  inclination_calculation(nei, i, vector_nei_lab_frame[i], i); // ...
  vector_average(vector_nei_lab_frame[nei], vector_nei_lab_frame[i], average_normal_nei_lab_frame[i]); // ...

  // calculate the roatation of quaternion version
  vector_to_quternion(average_normal_i_lab_frame[nei], reference_axis, inclination_quaternion_i_frame[nei]);

  // test
  // vectorFile << ">> nei: " << nei << " " << spin[nei] << " xyz: " << xyz[nei][0] <<" "<<xyz[nei][1]<<" "<<xyz[nei][2] <<\
  // "\n   vector_i: " << vector_i_lab_frame[i] << " vector_nei: " << vector_i_lab_frame[nei] << " ave: " << average_normal_i_lab_frame[nei] <<\
  // "\n   vector_nei: " << vector_nei_lab_frame[nei] << " vector_i: " << vector_nei_lab_frame[i] << " ave: " << average_normal_nei_lab_frame[i] <<\
  // "\n   quaternion rotation: " << inclination_quaternion_i_frame[nei] << "(" << 2*acos(inclination_quaternion_i_frame[nei][0])/MY_PI*180 << ")" << " && " << inclination_quaternion_nei_frame[i] <<\
  // "\n   rotation insymmetry: " << inclination_quaternion_symmetry_i_frame[nei] << " && " << inclination_quaternion_symmetry_nei_frame[i] <<\
  // "\n   tmp_energy: " << tmp_eng << "\n";

  // add the inclination energy into total site energy
  return compute_inclination_energy(inclination_quaternion_i_frame[nei]);
}

void AppPotts_eng::euler2quat(std::vector<double> eulerAngle, std::vector<double> & quternion_result)
{
    //Convert grain Euler angles to quaternion vector
    quternion_result[0]=cos(eulerAngle[1]/2.)*cos((eulerAngle[0]+eulerAngle[2])/2.);
    quternion_result[1]=sin(eulerAngle[1]/2.)*cos((eulerAngle[0]-eulerAngle[2])/2.);
    quternion_result[2]=sin(eulerAngle[1]/2.)*sin((eulerAngle[0]-eulerAngle[2])/2.);
    quternion_result[3]=cos(eulerAngle[1]/2.)*sin((eulerAngle[0]+eulerAngle[2])/2.);

}

/* -------------------------------------------------------
Convert symmetry matrix to quaternion form
------------------------------------------------------- */
void AppPotts_eng::mat2quat(const double O[3][3], double q[4])
{
    double q4 = 0;
    if ( (1 + O[0][0] + O[1][1] + O[2][2]) > 0) {
        q4 = sqrt(1 + O[0][0] + O[1][1] + O[2][2])/2;
        q[0] = q4;
        q[1] = (O[2][1] - O[1][2])/(4*q4);
        q[2] = (O[0][2] - O[2][0])/(4*q4);
        q[3] = (O[1][0] - O[0][1])/(4*q4);
    }
    else if ( (1 + O[0][0] - O[1][1] - O[2][2]) > 0) {
        q4 = sqrt(1 + O[0][0] - O[1][1] - O[2][2])/2;
        q[0] = (O[2][1] - O[1][2])/(4*q4);
        q[1] = q4;
        q[2] = (O[1][0] + O[0][1])/(4*q4);
        q[3] = (O[0][2] + O[2][0])/(4*q4);
    }
    else if ( (1 - O[0][0] + O[1][1] - O[2][2]) > 0) {
        q4 = sqrt(1 - O[0][0] + O[1][1] - O[2][2])/2;
        q[0] = (O[0][2] - O[2][0])/(4*q4);
        q[1] = (O[1][0] + O[0][1])/(4*q4);
        q[2] = q4;
        q[3] = (O[2][1] + O[1][2])/(4*q4);
    }
    else if ( (1 - O[0][0] - O[1][1] + O[2][2]) > 0) {
        q4 = sqrt(1 - O[0][0] - O[1][1] + O[2][2])/2;
        q[0] = (O[1][0] - O[0][1])/(4*q4);
        q[1] = (O[0][2] + O[2][0])/(4*q4);
        q[2] = (O[2][1] + O[1][2])/(4*q4);
        q[3] = q4;
    }
}

/* -------------------------------------------------------
Define the symmetry operator
------------------------------------------------------- */
void AppPotts_eng::symmat(double ***sym)
{
    //grow by number of symmetric operators
    (*sym) = new double*[Osym];

    //grow for symmetry quaternion vectors
    for (int o=0; o<Osym; o++)
        (*sym)[o] = new double[4];

    //buffer for quaternion
    double q[4];

    if (Osym == 24) {
        //cubic symmetry
        double SYM[24][3][3] =
                { {{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}},
                  {{ 1, 0, 0}, { 0,-1, 0}, { 0, 0,-1}},
                  {{ 1, 0, 0}, { 0, 0,-1}, { 0, 1, 0}},
                  {{ 1, 0, 0}, { 0, 0, 1}, { 0,-1, 0}},
                  {{-1, 0, 0}, { 0, 1, 0}, { 0, 0,-1}},
                  {{-1, 0, 0}, { 0,-1, 0}, { 0, 0, 1}},
                  {{-1, 0, 0}, { 0, 0,-1}, { 0,-1, 0}},
                  {{-1, 0, 0}, { 0, 0, 1}, { 0, 1, 0}},
                  {{ 0, 1, 0}, {-1, 0, 0}, { 0, 0, 1}},
                  {{ 0, 1, 0}, { 0, 0,-1}, {-1, 0, 0}},
                  {{ 0, 1, 0}, { 1, 0, 0}, { 0, 0,-1}},
                  {{ 0, 1, 0}, { 0, 0, 1}, { 1, 0, 0}},
                  {{ 0,-1, 0}, { 1, 0, 0}, { 0, 0, 1}},
                  {{ 0,-1, 0}, { 0, 0,-1}, { 1, 0, 0}},
                  {{ 0,-1, 0}, {-1, 0, 0}, { 0, 0,-1}},
                  {{ 0,-1, 0}, { 0, 0, 1}, {-1, 0, 0}},
                  {{ 0, 0, 1}, { 0, 1, 0}, {-1, 0, 0}},
                  {{ 0, 0, 1}, { 1, 0, 0}, { 0, 1, 0}},
                  {{ 0, 0, 1}, { 0,-1, 0}, { 1, 0, 0}},
                  {{ 0, 0, 1}, {-1, 0, 0}, { 0,-1, 0}},
                  {{ 0, 0,-1}, { 0, 1, 0}, { 1, 0, 0}},
                  {{ 0, 0,-1}, {-1, 0, 0}, { 0, 1, 0}},
                  {{ 0, 0,-1}, { 0,-1, 0}, {-1, 0, 0}},
                  {{ 0, 0,-1}, { 1, 0, 0}, { 0,-1, 0}} };

        //initialize global operator
        for (int o=0; o<Osym; o++) {
            mat2quat(SYM[o],q);
            for (int i=0; i<4; i++)
                (*sym)[o][i]=q[i];
        }
    }
    else if (Osym == 12) {
        //hexagonal symmetry
        double a = sqrt(3)/2;
        double SYM[12][3][3] =
                { {{ 1, 0, 0}, { 0, 1, 0}, { 0, 0, 1}},
                  {{-0.5, a, 0}, { -a,-0.5, 0}, { 0, 0, 1}},
                  {{-0.5, -a, 0}, { a,-0.5, 0}, { 0, 0, 1}},
                  {{ 0.5, a, 0}, { -a, 0.5, 0}, { 0, 0, 1}},
                  {{ -1, 0, 0}, { 0, -1, 0}, { 0, 0, 1}},
                  {{ 0.5, -a, 0}, { a, 0.5, 0}, { 0, 0, 1}},
                  {{-0.5, -a, 0}, { -a, 0.5, 0}, { 0, 0, -1}},
                  {{ 1, 0, 0}, { 0, -1, 0}, { 0, 0, -1}},
                  {{-0.5, a, 0}, { a, 0.5, 0}, { 0, 0, -1}},
                  {{ 0.5, a, 0}, { a,-0.5, 0}, { 0, 0, -1}},
                  {{ -1, 0, 0}, { 0, 1, 0}, { 0, 0, -1}},
                  {{ 0.5, -a, 0}, { -a,-0.5, 0}, { 0, 0, -1}} };

        //initialize global operator
        for (int o=0; o<Osym; o++) {
            mat2quat(SYM[o],q);
            for (int i=0; i<4; i++)
                (*sym)[o][i]=q[i];
        }
    }
}

void AppPotts_eng::quat_mult(const double qi[4], const double qj[4], double q[4])
{
  //Hamilton multiplication/product
  //multiplying quaternions and update
  q[0] = qi[0]*qj[0] - qi[1]*qj[1] - qi[2]*qj[2] - qi[3]*qj[3];
  q[1] = qi[0]*qj[1] + qi[1]*qj[0] + qi[2]*qj[3] - qi[3]*qj[2];
  q[2] = qi[0]*qj[2] - qi[1]*qj[3] + qi[2]*qj[0] + qi[3]*qj[1];
  q[3] = qi[0]*qj[3] + qi[1]*qj[2] - qi[2]*qj[1] + qi[3]*qj[0];
}

void AppPotts_eng::quaternions(const std::vector<double> qi, const std::vector<double> qj, std::vector<double> & q_result)
{
  double miso0, misom=MY_2PI;
  q_result.assign(4, 0);

  double q[4], qib[4], qjb[4], qmin[4]={0,0,0,0};
  double qi_array[4], qj_array[4];
  qi_array[0]=qi[0]; qi_array[1]=qi[1]; qi_array[2]=qi[2]; qi_array[3]=qi[3];
  qj_array[0]=qj[0]; qj_array[1]=qj[1]; qj_array[2]=qj[2]; qj_array[3]=qj[3];
  for (int o1=0; o1<Osym; o1++) {
    for (int o2=0; o2<Osym; o2++) {
      quat_mult(symquat[o1],qi_array,qib);
      quat_mult(symquat[o2],qj_array,qjb);

      //j-grain conjugate quaternion
      qjb[1]=-qjb[1];
      qjb[2]=-qjb[2];
      qjb[3]=-qjb[3];
      quat_mult(qib,qjb,q);
      miso0 = 2*acos(round(q[0]*1e5)/1e5);

      if (miso0 > MY_PI)
        miso0 = miso0-MY_2PI;
      if (fabs(miso0) < misom) {
        misom=fabs(miso0);
        qmin[0]=q[0]; qmin[1]=q[1]; qmin[2]=q[2]; qmin[3]=q[3];
        // std::ofstream vectorFile("Vector.txt", std::ios::out | std::ios::app);
        // vectorFile << ">>>> The qmin is " << qmin[0] << " " << qmin[1] << " " << qmin[2] << " " << qmin[3] <<"\n";
        // vectorFile.close();
      }
    }
  }

  miso0=2*acos(round(qmin[0]*1e5)/1e5);
  if (miso0 > MY_PI)
    miso0=miso0-MY_2PI;
  q_result[0] = qmin[0];
  q_result[1]=qmin[1];
  q_result[2]=qmin[2];
  q_result[3]=qmin[3];

//   cout<<fabs(miso0)<<endl;
  return;
}

/* -------------------------------------------------------
Calculate the rotation of quaternion version
between two vectors
------------------------------------------------------- */
void AppPotts_eng::vector_to_quternion(std::vector<double> & normals, std::vector<double> & reference_normals, std::vector<double> & q)
{
  // Initialization
  q.assign(4, 0);
  double cos_theta;

  // normals cross product reference_normals to get the rotation axis
  q[1] = normals[1]*reference_normals[2] - normals[2]*reference_normals[1];
  q[2] = normals[2]*reference_normals[0] - normals[0]*reference_normals[2];
  q[3] = normals[0]*reference_normals[1] - normals[1]*reference_normals[0];



  // normals with reference_normals to get the rotation angle
  cos_theta = round((normals[0]*reference_normals[0] + normals[1]*reference_normals[1] + normals[2]*reference_normals[2])*1e5)/1e5;
  q[0] = sqrt((cos_theta+1)/2);

  // get the quaternion by the rotation angle and axis
  q[1] = q[1]*sqrt((1-cos_theta)/2);
  q[2] = q[2]*sqrt((1-cos_theta)/2);
  q[3] = q[3]*sqrt((1-cos_theta)/2);

  if (q[1]==0 && q[2]==0 && q[3]==0){
    if (normals[2]==0 && reference_normals[2]==0){
      q[3] = 1;
    }else{
      q[1] = 1;
    }
  }

  // Normalized it
  double h = sqrt(q[1]*q[1] + q[2]*q[2] + q[3]*q[3]);
  if (h < 1){
    q[1] = q[1]/h;
    q[2] = q[2]/h;
    q[3] = q[3]/h;
  }
}


/* -------------------------------------------------------
Calculate the middle point vector by
averaging the two site vectors
------------------------------------------------------- */
void AppPotts_eng::vector_average(std::vector<double> vector_lab_frame_1, std::vector<double> vector_lab_frame_2, std::vector<double> & average_normal_lab_frame)
{
  average_normal_lab_frame.assign(3,0);
  average_normal_lab_frame[0] = vector_lab_frame_1[0] - vector_lab_frame_2[0];
  average_normal_lab_frame[1] = vector_lab_frame_1[1] - vector_lab_frame_2[1];
  average_normal_lab_frame[2] = vector_lab_frame_1[2] - vector_lab_frame_2[2];

  double normal_len = sqrt(average_normal_lab_frame[0]*average_normal_lab_frame[0]+average_normal_lab_frame[1]*average_normal_lab_frame[1]+average_normal_lab_frame[2]*average_normal_lab_frame[2]);

  if (normal_len < 1e-6) {
    average_normal_lab_frame[0] = vector_lab_frame_1[0];
    average_normal_lab_frame[1] = vector_lab_frame_1[1];
    average_normal_lab_frame[2] = vector_lab_frame_1[2];
    return;
  }
  average_normal_lab_frame[0] /= normal_len;
  average_normal_lab_frame[1] /= normal_len;
  average_normal_lab_frame[2] /= normal_len;

  return;
}

/* -------------------------------------------------------
Calculate the inclination on si by the crystal frame of i
------------------------------------------------------- */
void AppPotts_eng::inclination_calculation(int i, int si, std::vector<double> & vector_lab_frame, int center)
{
  // output the inclination vector by matrix method
  vector_lab_frame.assign(3,0);

  switch (dimension) {
    case 2:
      compute_normal_vector_2DLinear_matrix(i, si, vector_lab_frame, center);
      break;
    case 3:
      compute_normal_vector_3DLinear_matrix(i, si, vector_lab_frame, center);
      break;
  }
}

/* -------------------------------------------------------
Calculate the inclination from the vector
------------------------------------------------------- */
void AppPotts_eng::convert_inclination_refercrystal(int i, int i0, std::vector<double> & average_normal_lab_frame, std::vector<double> & inclination)
{
  inclination.assign(3,0);
  inclination[0] = 1.0 * (cos(phi1[i])*cos(phi2[i])-sin(phi1[i])*sin(phi2[i])*cos(Phi[i])) * average_normal_lab_frame[0] +\
                              (-cos(phi1[i])*sin(phi2[i])-sin(phi1[i])*cos(phi2[i])*cos(Phi[i])) * average_normal_lab_frame[1] +\
                              (sin(phi1[i])*sin(Phi[i])) * average_normal_lab_frame[2];
  inclination[1] = 1.0 * (sin(phi1[i])*cos(phi2[i])+cos(phi1[i])*sin(phi2[i])*cos(Phi[i])) * average_normal_lab_frame[0] +\
                              (-sin(phi1[i])*sin(phi2[i])+cos(phi1[i])*cos(phi2[i])*cos(Phi[i])) * average_normal_lab_frame[1] +\
                              (-cos(phi1[i])*sin(Phi[i])) * average_normal_lab_frame[2];
  inclination[2] = 1.0 * (sin(phi2[i])*sin(Phi[i])) * average_normal_lab_frame[0] +\
                              (cos(phi2[i])*sin(Phi[i])) * average_normal_lab_frame[1] +\
                              (cos(Phi[i])) * average_normal_lab_frame[2];
  return;
}

/* -------------------------------------------------------
Calculate the inclination energy at a voxel
------------------------------------------------------- */
double AppPotts_eng::compute_inclination_energy(std::vector<double> & input_i_frame)
{

  double inc_energy;

  // Tart calculation
  // inc_energy = 2*abs(acos(abs(input_i_frame[0]))-acos(abs(input_nei_frame[0])))/(MY_PI/2.85);

  double theta = 2*acos(abs(input_i_frame[0]));
  if (theta > MY_PI / 2) theta = MY_PI - theta;
  // inc_energy = 1.0 / (1 + E_delta) * (1 + E_delta * cos(E_m * theta));
  inc_energy = 1 + E_delta * cos(E_m * theta);

  // flat energy function
  // if (theta > 85.0 / 180 * MY_PI) inc_energy = 0.05;
  // else inc_energy = 1.0;

  return inc_energy;
}

void AppPotts_eng::get_boarder_place(int center, int & i_dis, int & j_dis, int & k_dis) {
  double x_min, y_min, z_min;
  double x_ind, y_ind, z_ind;
  x_min = 1.0*(domain->boxxhi - domain->boxxlo) / domain->procgrid[0];
  y_min = 1.0*(domain->boxyhi - domain->boxylo) / domain->procgrid[1];
  z_min = 1.0*(domain->boxzhi - domain->boxzlo) / domain->procgrid[2];
  i_dis = 0;
  j_dis = 0;
  k_dis = 0;

  if (dimension == 2) {
    x_ind = xyz[center][0];
    y_ind = xyz[center][1];
    while (x_ind>=x_min) x_ind = x_ind - x_min;
    while (y_ind>=y_min) y_ind = y_ind - y_min;

    if (x_ind < 1) j_dis = -1;
    else if (x_min - x_ind <= 1) j_dis = 1;
    if (y_ind < 1) i_dis = -1;
    else if (y_min - y_ind <= 1) i_dis = 1;

    // if ((xyz[center][0] - domain->subxlo) < 1 && (xyz[center][0] - domain->subxlo) >= 0) j_dis = -1;
    // else if (domain->subxhi - xyz[center][0] <= 1 && domain->subxhi - xyz[center][0] > 0) j_dis = 1;

    // if (xyz[center][1] - domain->subylo < 1 && xyz[center][1] - domain->subylo >= 0) i_dis = -1;
    // else if (domain->subyhi - xyz[center][1] <= 1 && domain->subyhi - xyz[center][1] > 0) i_dis = 1;

    if (i_dis==0 && j_dis==0) std::cout << "*** SPPARKS cannot find the correct boarder place 1 *** " << xyz[center][0] << " " << xyz[center][1] << std::endl;

  }
  else if (dimension == 3) {
    x_ind = xyz[center][0];
    y_ind = xyz[center][1];
    z_ind = xyz[center][2];
    while (x_ind>=x_min) x_ind = x_ind - x_min;
    while (y_ind>=y_min) y_ind = y_ind - y_min;
    while (z_ind>=z_min) z_ind = z_ind - z_min;

    if (x_ind < 1) j_dis = -1;
    else if (x_min - x_ind <= 1) j_dis = 1;
    if (y_ind < 1) i_dis = -1;
    else if (y_min - y_ind <= 1) i_dis = 1;
    if (z_ind < 1) k_dis = -1;
    else if (z_min - z_ind <= 1) k_dis = 1;

    // if (xyz[center][0] - domain->subxlo < 1 && xyz[center][0] - domain->subxlo >= 0) j_dis = -1;
    // else if (domain->subxhi - xyz[center][0] <= 1 && domain->subxhi - xyz[center][0] > 0) j_dis = 1;

    // if (xyz[center][1] - domain->subylo < 1 && xyz[center][1] - domain->subylo >= 0) i_dis = -1;
    // else if (domain->subyhi - xyz[center][1] <= 1 && domain->subyhi - xyz[center][1] > 0) i_dis = 1;

    // if (xyz[center][2] - domain->subzlo < 1 && xyz[center][2] - domain->subzlo >= 0) k_dis = -1;
    // else if (domain->subzhi - xyz[center][2] <= 1 && domain->subzhi - xyz[center][2] > 0) k_dis = 1;

    if (i_dis==0 && j_dis==0 && k_dis==0) std::cout << "*** SPPARKS cannot find the correct boarder place *** " << xyz[center][0] << " " << xyz[center][1] << " " << xyz[center][2] << std::endl;

  }

  return;
}

void AppPotts_eng::compute_normal_vector_2DLinear_matrix(int i, int si, std::vector<double> & vector_lab_frame, int center)
{
  // Get the normals at si, all spin == spin[i] is 1

  int sval, local_x, local_y;
  int j_dis = 0, i_dis = 0, k_dis = 0;
  double inc_x=0, inc_y=0, vec_len;

  if (numneigh[si] == (2*interval+3)*(2*interval+3)-1){
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j;
      if (ni < ((2*interval+3)*(2*interval+3)-1)/2) {
        norm_i = ni / (2*interval+3);
        norm_j = ni - (2*interval+3)*norm_i;
      }
      else {
        norm_i = (ni+1) / (2*interval+3);
        norm_j = (ni+1) - (2*interval+3)*norm_i;
      }

      // vectorFile << " sval: " << sval << " factor_i: " << factor_i << " matrix_i: " << linear_vector_matrix_i[norm_i][norm_j] << " matrix_j: " << linear_vector_matrix_j[norm_i][norm_j];

      inc_x += factor_i*linear_vector_matrix_i[norm_i][norm_j];
      inc_y += factor_i*linear_vector_matrix_j[norm_i][norm_j];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
  } else if (numneigh[si] == (2*interval+3)*(2*interval+2)-1) {
    // std::cout << "numneigh["<<xyz[i][0]<<"]["<<xyz[i][1]<<"]: " << numneigh[i] << std::endl;
    get_boarder_place(center, i_dis, j_dis, k_dis);
    if (i_dis*j_dis!=0) {
      if (xyz[si][1]==xyz[center][1]) i_dis = 0;
      if (xyz[si][0]==xyz[center][0]) j_dis = 0;
    }
    if (i_dis*j_dis!=0) {
      i_dis = 0; j_dis = 0; k_dis = 0;
      get_boarder_place(si, i_dis, j_dis, k_dis);
      i_dis = -i_dis; j_dis = -j_dis; k_dis = -k_dis;
    }
    if (i_dis*j_dis!=0) std::cout <<"Unexpected bool issue " << "ij_dis: " << i_dis<<" "<<j_dis << " si: " << xyz[si][0] << " " << xyz[si][1] << " center: " << xyz[center][0] << " " << xyz[center][1] << std::endl;
    // j_dis = xyz[i][0]-xyz[center][0];
    // i_dis = xyz[i][1]-xyz[center][1];
    // std::cout << "j_dis: " << j_dis << " i_dis: " << i_dis << std::endl;
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j;
      if (ni < (interval+1-(i_dis<0?1:0))*(2*interval+3-abs(j_dis))+interval+1-(j_dis<0?1:0)) {
        norm_i = ni / (2*interval+3-abs(j_dis)) + (i_dis<0? 1:0);
        norm_j = ni - (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0? 1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_i = (ni+1) / (2*interval+3-abs(j_dis)) + (i_dis<0? 1:0);
        norm_j = (ni+1) - (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0? 1:0)) + (j_dis<0?1:0);
      }

      // vectorFile << " sval: " << sval << " factor_i: " << factor_i << " matrix_i: " << linear_vector_matrix_i[norm_i][norm_j] << " matrix_j: " << linear_vector_matrix_j[norm_i][norm_j];
      // std::cout << "ni: " << ni << " norm_i: " << norm_i << " norm_j: " << norm_j << std::endl;
      inc_x += factor_i*linear_vector_matrix_i[norm_i][norm_j];
      inc_y += factor_i*linear_vector_matrix_j[norm_i][norm_j];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
      // std::cout << "norm_i: " << norm_i << " norm_j: " << norm_j << std::endl;
    }
    // std::cout << "Im done !" << std::endl;
  } else if (numneigh[si] == (2*interval+2)*(2*interval+2)-1) {
    // j_dis = xyz[si][0]-xyz[center][0];
    // i_dis = xyz[si][1]-xyz[center][1];
    get_boarder_place(si, i_dis, j_dis, k_dis);
    // std::cout << "j_dis: " << j_dis << " i_dis: " << i_dis << std::endl;
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j;
      if (ni < (interval+1-(i_dis<0?1:0))*(2*interval+3-abs(j_dis))+interval+1-(j_dis<0?1:0)) {
        norm_i = ni / (2*interval+3-abs(j_dis)) + (i_dis<0? 1:0);
        norm_j = ni - (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0? 1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_i = (ni+1) / (2*interval+3-abs(j_dis)) + (i_dis<0? 1:0);
        norm_j = (ni+1) - (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0? 1:0)) + (j_dis<0?1:0);
      }

      inc_x += factor_i*linear_vector_matrix_i[norm_i][norm_j];
      inc_y += factor_i*linear_vector_matrix_j[norm_i][norm_j];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
  }

  vec_len = sqrt(inc_x*inc_x+inc_y*inc_y);
  if (vec_len > 1e-5) {
    vector_lab_frame[0] = -inc_y/vec_len;
    vector_lab_frame[1] = -inc_x/vec_len;
    vector_lab_frame[2] = 0;
  }


  // Test
  // std::cout<<"vec_len: "<<vec_len <<" vector_i_lab_frame[0]:"<<vector_i_lab_frame[si][0]<< " local x y "<<local_x <<" "<<local_y << std::endl;

  return;
}

void AppPotts_eng::compute_normal_vector_3DLinear_matrix(int i, int si, std::vector<double> & vector_lab_frame, int center)
{
  // Get the normals at si, all spin == spin[i] is 1

  int sval, local_x, local_y, local_z;;
  int j_dis = 0, i_dis = 0, k_dis = 0;
  int j_dis_si = 0, i_dis_si = 0, k_dis_si = 0;
  double inc_x=0, inc_y=0, inc_z=0, vec_len;

  if (numneigh[si] == (2*interval+3)*(2*interval+3)*(2*interval+3)-1){
    // std::cout << " 3DL0 ";
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j, norm_k;
      if (ni < ((2*interval+3)*(2*interval+3)*(2*interval+3)-1)/2) {
        norm_k = ni / ((2*interval+3)*(2*interval+3));
        norm_i = (ni - ((2*interval+3)*(2*interval+3))*norm_k) / (2*interval+3);
        norm_j = ni - ((2*interval+3)*(2*interval+3))*norm_k - (2*interval+3) * norm_i;
      }
      else {
        norm_k = (ni+1) / ((2*interval+3)*(2*interval+3));
        norm_i = ((ni+1) - ((2*interval+3)*(2*interval+3))*norm_k) / (2*interval+3);
        norm_j = (ni+1) - ((2*interval+3)*(2*interval+3))*norm_k - (2*interval+3) * norm_i;
      }

      // vectorFile << " sval: " << sval << " factor_i: " << factor_i << " matrix_i: " << linear_vector_matrix_i[norm_i][norm_j] << " matrix_j: " << linear_vector_matrix_j[norm_i][norm_j];

      inc_x += factor_i*linear_vector_matrix_i_3d[norm_i][norm_j][norm_k];
      inc_y += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      inc_z += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
    // std::cout << "√ ";
  } else if (numneigh[si] == (2*interval+3)*(2*interval+3)*(2*interval+2)-1) {
    // std::cout << " 3DL1 ";
    get_boarder_place(center, i_dis, j_dis, k_dis);
    if (!(bool(i_dis)^bool(j_dis)^bool(k_dis)) || bool(i_dis*j_dis*k_dis)) {
      if (xyz[si][1]==xyz[center][1]) i_dis = 0;
      if (xyz[si][0]==xyz[center][0]) j_dis = 0;
      if (xyz[si][2]==xyz[center][2]) k_dis = 0;
    }
    if (!(bool(i_dis)^bool(j_dis)^bool(k_dis)) || bool(i_dis*j_dis*k_dis)) {
      i_dis_si = 0; j_dis_si = 0; k_dis_si = 0;
      get_boarder_place(si, i_dis_si, j_dis_si, k_dis_si);
      i_dis_si = -i_dis_si; j_dis_si = -j_dis_si; k_dis_si = -k_dis_si;
      i_dis = i_dis * (i_dis == i_dis_si);
      j_dis = j_dis * (j_dis == j_dis_si);
      k_dis = k_dis * (k_dis == k_dis_si);
    }
    if (!(bool(i_dis)^bool(j_dis)^bool(k_dis)) || bool(i_dis*j_dis*k_dis)) std::cout << "Unexpected bool issue!!" << "jik_dis: " << j_dis<<" "<<i_dis <<" "<<k_dis<<" jik_dis_si: "<<j_dis_si<<" "<<i_dis_si<<" "<<k_dis_si<< " si: " << xyz[si][0] << " " << xyz[si][1] <<" "<<xyz[si][2] << " center: " << xyz[center][0] << " " << xyz[center][1]<<" "<<xyz[center][2] << " expected: "<<numneigh[si] << std::endl;
    // j_dis = xyz[i][0]-xyz[center][0];
    // i_dis = xyz[i][1]-xyz[center][1];
    // std::cout << "j_dis: " << j_dis << " i_dis: " << i_dis << std::endl;
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j, norm_k;
      if (ni < (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(interval+1-(k_dis<0?1:0))+
      (2*interval+3-abs(j_dis))*(interval+1-(i_dis<0?1:0))+interval+1-(j_dis<0?1:0)) {
        norm_k = ni / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = (ni -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = ni - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_k = (ni+1) / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = ((ni+1) -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = (ni+1) - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }

      // vectorFile << " sval: " << sval << " factor_i: " << factor_i << " matrix_i: " << linear_vector_matrix_i[norm_i][norm_j] << " matrix_j: " << linear_vector_matrix_j[norm_i][norm_j];
      // std::cout << "ni: " << ni << " norm_i: " << norm_i << " norm_j: " << norm_j << std::endl;
      inc_x += factor_i*linear_vector_matrix_i_3d[norm_i][norm_j][norm_k];
      inc_y += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      inc_z += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
    // std::cout << "√ ";
    // std::cout << "Im done !" << std::endl;
  } else if (numneigh[si] == (2*interval+3)*(2*interval+2)*(2*interval+2)-1) {
    // std::cout << " 3DL2 ";
    get_boarder_place(center, i_dis, j_dis, k_dis);
    if (!((bool(i_dis) && bool(j_dis)) ^ (bool(j_dis) && bool(k_dis)) ^ (bool(i_dis) && bool(k_dis))) || bool(i_dis*j_dis*k_dis)) {
      if (xyz[si][1]==xyz[center][1]) i_dis = 0;
      if (xyz[si][0]==xyz[center][0]) j_dis = 0;
      if (xyz[si][2]==xyz[center][2]) k_dis = 0;
    }
    if (!((bool(i_dis) && bool(j_dis)) ^ (bool(j_dis) && bool(k_dis)) ^ (bool(i_dis) && bool(k_dis))) || bool(i_dis*j_dis*k_dis)) {
      i_dis_si = 0; j_dis_si = 0; k_dis_si = 0;
      get_boarder_place(si, i_dis_si, j_dis_si, k_dis_si);
      i_dis_si = -i_dis_si; j_dis_si = -j_dis_si; k_dis_si = -k_dis_si;
      i_dis = i_dis * (i_dis == i_dis_si);
      j_dis = j_dis * (j_dis == j_dis_si);
      k_dis = k_dis * (k_dis == k_dis_si);
    }
    if (!((bool(i_dis) && bool(j_dis)) ^ (bool(j_dis) && bool(k_dis)) ^ (bool(i_dis) && bool(k_dis))) || bool(i_dis*j_dis*k_dis)) std::cout << "Unexpected bool issue!!" << "jik_dis: " << j_dis<<" "<<i_dis <<" "<<k_dis<<" jik_dis_si: "<<j_dis_si<<" "<<i_dis_si<<" "<<k_dis_si<< " si: " << xyz[si][0] << " " << xyz[si][1] <<" "<<xyz[si][2] << " center: " << xyz[center][0] << " " << xyz[center][1]<<" "<<xyz[center][2] << " expected: "<<numneigh[si] << std::endl;
    // j_dis = xyz[i][0]-xyz[center][0];
    // i_dis = xyz[i][1]-xyz[center][1];
    // std::cout << "j_dis: " << j_dis << " i_dis: " << i_dis << std::endl;
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j, norm_k;
      if (ni < (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(interval+1-(k_dis<0?1:0))+
      (2*interval+3-abs(j_dis))*(interval+1-(i_dis<0?1:0))+interval+1-(j_dis<0?1:0)) {
        norm_k = ni / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = (ni -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = ni - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_k = (ni+1) / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = ((ni+1) -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = (ni+1) - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }

      // vectorFile << " sval: " << sval << " factor_i: " << factor_i << " matrix_i: " << linear_vector_matrix_i[norm_i][norm_j] << " matrix_j: " << linear_vector_matrix_j[norm_i][norm_j];
      // std::cout << "ni: " << ni << " norm_i: " << norm_i << " norm_j: " << norm_j << std::endl;
      inc_x += factor_i*linear_vector_matrix_i_3d[norm_i][norm_j][norm_k];
      inc_y += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      inc_z += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
      // std::cout << "norm_i: " << norm_i << " norm_j: " << norm_j << std::endl;
    }
    // std::cout << "√ ";
    // std::cout << "Im done !" << std::endl;
  } else if (numneigh[si] == (2*interval+2)*(2*interval+2)*(2*interval+2)-1) {
    // std::cout << " 3DL3 ";
    // j_dis = xyz[si][0]-xyz[center][0];
    // i_dis = xyz[si][1]-xyz[center][1];
    // k_dis = xyz[si][2]-xyz[center][2];
    get_boarder_place(center, i_dis, j_dis, k_dis);
    // std::cout << " jd: " << j_dis << " id: " << i_dis << " kd: " << k_dis << std::endl;
    for (int ni = 0; ni < numneigh[si]; ni++) {
      sval = neighbor[si][ni];
      int factor_i;
      if (spin[i]==spin[si]) factor_i=int(spin[sval]==spin[i]);
      else factor_i=int(spin[sval]!=spin[i]);

      int norm_i, norm_j, norm_k;
      if (ni < (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(interval+1-(k_dis<0?1:0))+
      (2*interval+3-abs(j_dis))*(interval+1-(i_dis<0?1:0))+interval+1-(j_dis<0?1:0)) {
        norm_k = ni / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = (ni -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = ni - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }
      else {
        norm_k = (ni+1) / ((2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))) + (k_dis<0? 1:0);
        norm_i = ((ni+1) -(2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0))) / (2*interval+3-abs(j_dis)) +
                 (i_dis<0?1:0);
        norm_j = (ni+1) - (2*interval+3-abs(i_dis))*(2*interval+3-abs(j_dis))*(norm_k-(k_dis<0? 1:0)) -
                      (2*interval+3-abs(j_dis))*(norm_i-(i_dis<0?1:0)) + (j_dis<0?1:0);
      }

      // std::cout << " ni: " << ni << " norm_i: " << norm_i << " norm_j: " << norm_j << " norm_k: " << norm_k << std::endl;
      inc_x += factor_i*linear_vector_matrix_i_3d[norm_i][norm_j][norm_k];
      inc_y += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      inc_z += factor_i*linear_vector_matrix_j_3d[norm_i][norm_j][norm_k];
      // We don't need to count on the central site because linear_vector_matrix_i/j = 0 at central site
    }
    // std::cout << "√ ";
  }


  vec_len = sqrt(inc_x*inc_x+inc_y*inc_y+inc_z*inc_z);
  if (vec_len > 1e-5) {
    double ff1 = 1/vec_len;
    vector_lab_frame[0] = -inc_z*ff1;
    vector_lab_frame[1] = -inc_y*ff1;
    vector_lab_frame[2] = -inc_x*ff1;
  }

  return;
}

int AppPotts_eng::int_mod(int base, int divide) {
  return int(base - floor(double(base)/divide)*divide);
}

/* -------------------------------------------------------
rKMC method
perform a site event with no null bin rejection
flip to random neighbor spin without null bin
------------------------------------------------------- */
void AppPotts_eng::site_event_rejection(int i, RandomPark *random)
{
  if (storage_flag!=0 && ((nsweeps)%int(maxneigh*storage_flag) == 0) && (cleanstep != nsweeps)) {
    std::cout << "naccept = " << naccept << std::endl;
    std::cout << "boundary sites: " << 1.0*(num_read+num_calc)/(nx*ny*nz*maxneigh)*100 << "%, p(calc_sites): " << 1.0*num_calc/(num_read+num_calc)*100 << "%, e/s is " << test_error_incl/(num_read+num_calc) << std::endl;
    cleanstep = nsweeps;
    num_read=0;
    num_calc=0;
    test_error_incl=0.0;
    for ( int cleani=0; cleani<inc_storage.size(); cleani++) {
      inc_storage[cleani].clear();
    }
  }
  // no events for a pinned site
  if (spin[i] > nspins) return;

  int oldstate=spin[i];
  double iphi[3]={phi1[i],Phi[i],phi2[i]};

  // events = spin flips to neighboring site different than self

  int j,nei;
  int nevent = 0;
  float tmp_error =0.0;
  get_nearest_neighbor(i, interval);

  //Nearest-neighbor sampling
  for (j = 0; j < numneigh_near; j++) {
    nei=neighbor_near[j];
    if (spin[i]==spin[nei])
      continue;
    sites[nevent++]=nei;
  }

  if (nevent == 0) return;
  int iran = (int) (nevent*random->uniform()); // what is this mean??  ->iniform()??
  if (iran >= nevent) iran = nevent-1;

  double einitial = site_energy(i); //, qold[4];  // qold?? how to neutralize it?
  // double einitial_add = compute_inclination_energy(inclination_i_frame[i]);
  // einitial += einitial_add;

  spin[i] = spin[sites[iran]];
  phi1[i] = phi1[sites[iran]];
  phi2[i] = phi2[sites[iran]];
  Phi[i] = Phi[sites[iran]];

  double efinal = site_energy(i);
  // double efinal_add = compute_inclination_energy(inclination_nei_frame[i]);
  // efinal += efinal_add;

  // Test output
  // init_vector(i, oldstate, einitial-einitial_add, einitial_add, efinal-efinal_add, efinal_add);

  //Determing misorientation between ij states to
  //calculate mobility
  double thetar;
  int smin = MIN(oldstate,spin[i]);
  int smax = MAX(oldstate,spin[i]);

  std::pair <int, int> spins = std::make_pair(smin,smax);

  thetar=misos[spins];

  // now the p0 is 1 in aniso
  double p0=1;
  // double p0=current_gamma;

  //Check for isotropic case
  // if (thetam<1e-8) p0=1;

  // accept or reject via Boltzmann criterion
  if (efinal <= einitial) {
    if ((thetar < 1e-8) || (random->uniform() < p0)) {
    }
    else {
      spin[i] = oldstate;
      phi1[i] = iphi[0];
      phi2[i] = iphi[2];
      Phi[i] = iphi[1];
    }
  }
  else if (temperature == 0.0) {
    spin[i] = oldstate;
    phi1[i] = iphi[0];
    phi2[i] = iphi[2];
    Phi[i] = iphi[1];
  }
  else if (random->uniform() > p0*exp((einitial-efinal)*t_inverse/current_gamma)) {
    spin[i] = oldstate;
    phi1[i] = iphi[0];
    phi2[i] = iphi[2];
    Phi[i] = iphi[1];
  }
  else {
  }

  if (spin[i] != oldstate) {
    naccept++;
  }
}
