/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/


#include <iostream>
#include <algorithm>
#include <cassert>
#include <math.h>

#include "params.h"
#include "hdf5ParticleIO.h"
#include "user.h"
#include "lattice.h"

//
// help message
//
void print_usage() {
  using namespace std;
  clog_one(warn) << "Initial data generator for Kelvin-Helmholtz test in" <<
  gdimension << "D" << endl << "Usage: ./kh_generator <parameter-file.par>"
  << endl;
}

//
// derived parameters
//
static int64_t nparticlesproc;        // number of particles per proc
static double rho_1, rho_2;           // densities
static double vx_1, vx_2;             // velocities
static double dely;                   // perturbation in vy
static std::string initial_data_file; // = initial_data_prefix + ".h5part"

void set_derived_params() {
  using namespace std;
  using namespace param;

  // total number of particles
  if(gdimension==2){
    SET_PARAM(nparticles, sqrt_nparticles*sqrt_nparticles);
  } else if(gdimension==2){
    SET_PARAM(nparticles, cbrt_nparticles*cbrt_nparticles*cbrt_nparticles);
  }

  // particle spacing and smoothing length
  if(gdimension==2){
    SET_PARAM(sph_separation, 1.0/((double)sqrt_nparticles-1.));
    SET_PARAM(sph_smoothing_length, (sph_separation*4.)); // TODO: use sph_eta instead
  } else if(gdimension==3){
    SET_PARAM(sph_separation, 1.0/((double)cbrt_nparticles-1.));
    SET_PARAM(sph_smoothing_length, (sph_separation*3.)); // TODO: use sph_eta instead
  }

  // test selector
  switch (khtest_num) {
    case (1):
      // -- top and bottom | middle -- //
      rho_1      = 1.0;      rho_2      = 2.0;
      vx_1       = 0.5;      vx_2       =-0.5;
      dely       = 0.01;
      break;

    case (2):
      rho_1      = 32.0;     rho_2      = 64.0;
      vx_1       = 0.1;      vx_2       =-0.1;
      dely       = 0.01;
      break;

    case (3):
      rho_1      = 32.0;     rho_2      = 64.0;
      vx_1       = 0.1;      vx_2       =-0.1;
      dely       = 0.002;
      break;

    default:
      clog_one(error) << "ERROR: invalid test (" << khtest_num << ")." << endl;
      MPI_Finalize();
      exit(-1);
  }

  // file to be generated
  std::ostringstream oss;
  oss << initial_data_prefix << ".h5part";
  initial_data_file = oss.str();

}

//----------------------------------------------------------------------------//
int main(int argc, char * argv[]){
  using namespace std;
  using namespace param;

  // launch MPI
  int rank, size, provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);
  clog_set_output_rank(0);

  // check options list: exactly one option is allowed
  if (argc != 2) {
    print_usage();
    MPI_Finalize();
    exit(0);
  }

  // set simulation parameters
  param::mpi_read_params(argv[1]);
  set_derived_params();

  // screen output
  clog_one(info) << "Kelvin-Helmholtz test #" << khtest_num << " in " << gdimension
         << "D:" << endl << " - number of particles: " << nparticles
         << endl << " - particles per core:  " << nparticlesproc << endl
         << " - generated initial data file: " << initial_data_file << endl;

  // Central coordinates: in most cases this should be centered at (0,0,0)
  point_t cbox_max = 0.;
  point_t cbox_min = 0.;
  point_t tbox_max = 0.;
  point_t tbox_min = 0.;
  point_t bbox_max = 0.;
  point_t bbox_min = 0.;

  // set the dimensions of the rectangles: top, middle, bottom
  // currently 3 rectangles of within one cube length 1.0
  if(gdimension==2){
    cbox_max[0] = (sqrt_nparticles-1.)*sph_separation/2.;
    cbox_max[1] = cbox_max[0]/3.;

    cbox_min[0] = -cbox_max[0];
    cbox_min[1] = -cbox_max[1];
  } else{
    cbox_max[0] = (cbrt_nparticles-1.)*sph_separation/2.;
    cbox_max[1] = cbox_max[0]/3.;
    cbox_max[2] = cbox_max[0];

    cbox_min[0] = -cbox_max[0];
    cbox_min[1] = -cbox_max[1];
    cbox_min[2] = -cbox_max[0];
  }

  tbox_max[0] = cbox_max[0];
  tbox_min[0] = cbox_min[0];

  bbox_max[0] = cbox_max[0];
  bbox_min[0] = cbox_min[0];

  tbox_max[1] = 3.*cbox_max[1]+sph_separation/1000.;
  tbox_min[1] = cbox_max[1]+sph_separation/1000.;

  bbox_max[1] = cbox_min[1]-sph_separation/1000.;
  bbox_min[1] = 3.*cbox_min[1]-sph_separation/1000.;

  if(gdimension>2){
    tbox_max[2] = cbox_max[2];
    tbox_min[2] = cbox_min[2];

    bbox_max[2] = cbox_max[2];
    bbox_min[2] = cbox_min[2];
  }

  // allocate arrays
  int64_t tparticles = 0;
  int64_t parts_each_rect = 0;

  // set defined values
  const double lambda = 0.5;
  const double pressure = 2.5;

  bool count = true;
  tparticles = generate_lattice(lattice_type,2,cbox_min,cbox_max,sph_separation,0,count);
  parts_each_rect = tparticles;
  tparticles += generate_lattice(lattice_type,2,tbox_min,tbox_max,sph_separation,tparticles-1,count);
  tparticles += generate_lattice(lattice_type,2,bbox_min,bbox_max,sph_separation,tparticles-1,count);
  count = false;

  // Initialize the arrays to be filled later
  // Position
  double* x = new double[tparticles]();
  double* y = new double[tparticles]();
  double* z = new double[tparticles]();
  // Velocity
  double* vx = new double[tparticles]();
  double* vy = new double[tparticles]();
  double* vz = new double[tparticles]();
  // Acceleration
  double* ax = new double[tparticles]();
  double* ay = new double[tparticles]();
  double* az = new double[tparticles]();
  // Smoothing length
  double* h = new double[tparticles]();
  // Density
  double* rho = new double[tparticles]();
  // Internal Energy
  double* u = new double[tparticles]();
  // Pressure
  double* P = new double[tparticles]();
  // Mass
  double* m = new double[tparticles]();
  // Id
  int64_t* id = new int64_t[tparticles]();
  // Timestep
  double* dt = new double[tparticles]();

  tparticles = generate_lattice(lattice_type,2,cbox_min,cbox_max,sph_separation,0,count,x,y,z);
  tparticles += generate_lattice(lattice_type,2,tbox_min,tbox_max,sph_separation,tparticles-1,count,x,y,z);
  tparticles += generate_lattice(lattice_type,2,bbox_min,bbox_max,sph_separation,tparticles-1,count,x,y,z);

  // particle id number
  int64_t posid = 0;

  // max. value for the speed of sound
  double cs = sqrt(poly_gamma*max(pressure/rho_1,pressure/rho_2));

  // The value for constant timestep
  double timestep = 0.5*sph_separation/cs;

  for(int64_t part=0; part<tparticles; ++part){
    id[part] = posid++;
    P[part] = pressure;
    if(y[part] < cbox_min[1] || y[part] > cbox_max[1]){
      rho[part] = rho_2;
      vx[part] = vx_2;
      m[part] = rho[part]/(double)parts_each_rect/3.;
    } else{
      rho[part] = rho_1;
      vx[part] = vx_1;
      m[part] = rho[part]/(double)parts_each_rect/3.;
    }

    //perturbation from 0 vy
    vy[part] = dely*sin(2.*M_PI*x[part]/lambda);

    // compute internal energy using gamma-law eos
    u[part] = P[part]/(poly_gamma-1.)/rho[part];

    // particle smoothing length
    h[part] = sph_smoothing_length;


  } // for part=0..nparticles

  clog_one(info) << "Real number of particles: " << tparticles << std::endl;
  // delete the output file if exists
  remove(initial_data_file.c_str());

  // Header data
  // the number of particles = nparticles
  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(initial_data_file,MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",tparticles);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",gdimension);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  //testDataSet.writeDatasetAttributeArray("name","string",simName);
  testDataSet.closeFile();

  testDataSet.openFile(MPI_COMM_WORLD);
  testDataSet.setTimeStep(0);

  Flecsi_Sim_IO::Variable _d1,_d2,_d3;

  _d1.createVariable("x",Flecsi_Sim_IO::point,"double",tparticles,x);
  _d2.createVariable("y",Flecsi_Sim_IO::point,"double",tparticles,y);
  _d3.createVariable("z",Flecsi_Sim_IO::point,"double",tparticles,z);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("vx",Flecsi_Sim_IO::point,"double",tparticles,vx);
  _d2.createVariable("vy",Flecsi_Sim_IO::point,"double",tparticles,vy);
  //_d3.createVariable("vz",Flecsi_Sim_IO::point,"double",nparticlesproc,vz);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  //_d1.createVariable("ax",Flecsi_Sim_IO::point,"double",nparticlesproc,ax);
  //_d2.createVariable("ay",Flecsi_Sim_IO::point,"double",nparticlesproc,ay);
  //_d3.createVariable("az",Flecsi_Sim_IO::point,"double",nparticlesproc,az);

  //testDataSet.vars.push_back(_d1);
  //testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  //testDataSet.writeVariables();


  _d1.createVariable("h",Flecsi_Sim_IO::point,"double",tparticles,h);
  _d2.createVariable("rho",Flecsi_Sim_IO::point,"double",tparticles,rho);
  _d3.createVariable("u",Flecsi_Sim_IO::point,"double",tparticles,u);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("P",Flecsi_Sim_IO::point,"double",tparticles,P);
  _d2.createVariable("m",Flecsi_Sim_IO::point,"double",tparticles,m);
  _d3.createVariable("id",Flecsi_Sim_IO::point,"int64_t",tparticles,id);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  testDataSet.closeFile();

  delete[] x;
  delete[] y;
  delete[] z;
  delete[] vx;
  delete[] vy;
  delete[] vz;
  delete[] ax;
  delete[] ay;
  delete[] az;
  delete[] h;
  delete[] rho;
  delete[] u;
  delete[] P;
  delete[] m;
  delete[] id;
  delete[] dt;

  MPI_Finalize();
  return 0;
}
