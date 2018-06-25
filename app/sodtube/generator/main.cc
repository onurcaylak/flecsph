/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

#include <iostream>
#include <algorithm>
#include <cassert>

#include "hdf5ParticleIO.h"
#include "kernel.h"


namespace simulation_params {
  int64_t nparticles;        // global number of particles
  double ldistance;          // particles spacing
  double localgamma;         // polytropic index
  double smoothing_length;   // constant smoothing length
  double length;             // Length of the tube

  // test conditions for two sides of the domain
  int    sodtest_num;            // which Sod test to generate
  double rho_l, rho_r;           // densities
  double vx_l, vx_r;             // velocities
  double pressure_l, pressure_r; // pressures
  double distance_l, distance_r; // distances

  // output filename
  const char* fileprefix = "hdf5_sodtube";
  char output_filename[128];
}


double 
cubic_spline_kernel(
    double r, 
    double h)
{
  double rh = r/h;
  // Default 1D case
  double sigma = 2./(3.*h);
  if(gdimension == 2){
    sigma = 10./(7.*M_PI*h*h);
  }
  double result = 0.; 
  if (0.0 <= rh && rh <= 1.0) {
    result = 1.0 - 1.5*rh*rh + .75*rh*rh*rh;
    result *= sigma;  
  }else if (1.0 < rh && rh <= 2.0) {
    result = 0.25 * (2-rh)*(2-rh)*(2-rh);
    result *= sigma; 
  }
  return result;
} // kernel

double compute_density(int source,double*x,double*h,double*mass)
{
  using namespace std;
  using namespace simulation_params;
  
  double density = 0.;
  for(int part = 0; part < nparticles; ++part){
    double dist = fabs(x[part]-x[source]);
    density += mass[part]*cubic_spline_kernel(dist, 
      (h[part]+h[source])/2.);
  }
  return density;
}


double bisection_smoothing_length(double * x,double * h, double * mass,double * rho)
{
  using namespace std;
  using namespace simulation_params;

  double* min_h = new double[nparticles]();
  double* max_h = new double[nparticles]();

  int niter = 100;

  for(int part = 0 ; part < nparticles; part++){
    h[part] = 1.4*(mass[part]/rho[part]);
    min_h[part] = h[part]-1.;
    max_h[part] = h[part]+1.;
    h[part] = (min_h[part]+max_h[part])/2.;
    if(part == nparticles/2){
      std::cout<<"init h="<<h[part]<<" min="<<min_h[part]<<" max="<<max_h[part]<<std::endl;
    }
  }

  for(int i = 0 ; i < niter; ++i)
  {
    // For each particles 
    for(int part = 0; part < nparticles; ++part){
      // Compute the density 
      double density = compute_density(part,x,h,mass);
      if(part == nparticles/2){
        std::cout<<"density="<<density<<" rho="<<rho[part]<<" h="<<h[part]<<std::endl;
      }
      // If not ok, change the h 
      if(density > rho[part])
      {
        max_h[part] = h[part];
      }else{
        min_h[part] = h[part];
      }
    }
    for(int part = 0 ; part < nparticles; ++part){
      h[part] = (min_h[part]+max_h[part])/2.;
    }
  }
  delete[] min_h;
  delete[] max_h;
}


//
// setup parameter defaults
//
void set_default_param(int rank, int size) {
  using namespace simulation_params;

  // number of particles
  nparticles = 405;

  // equation of state parameters (one so far)
  localgamma = 1.4;

  // run Sod test 1 by default
  sodtest_num = 1;

  // Length of the tube 
  length = 1.;
}


//
// help message
//
void print_usage(int rank) {
  using namespace std;
  if (rank == 0)
    cout << "Initial data generator for Sod shocktube test in 1D" << endl
         << "Usage: ./sodtube_generator [OPTIONS]" << endl
         << " -h: this help" << endl
         << " -n <number of particles>" << endl
         << " -t <Sod test (integer from 1 to 5)>" << endl;
}


//
// option parser
//
void parse_command_line_options(int rank, int size, int argc, char* argv[]) {
  using namespace std;
  using namespace simulation_params;

  for (int i=1; i<argc; ++i)
    if (argv[i][0] == '-')
      switch(argv[i][1]) {
      case 'h':
        print_usage(rank);
        MPI_Finalize();
        exit(0);
        break;

      case 'n':
        nparticles = atoll(argv[++i]);
        break;

      case 't':
        sodtest_num = atoi(argv[++i]);
        assert(sodtest_num>=1 && sodtest_num<=6);
        break;

      default:
        if (rank == 0)
          cerr << "ERROR: unknown option '-" << argv[i][1] << "'" << endl;
        MPI_Finalize();
        exit(-1);

      } // switch argv[i][1]

}


//
// setup simulation parameters
//
void set_param(int rank, int size) {
  using namespace std;
  using namespace simulation_params;

  // test selector
  switch (sodtest_num) {
    case (1):
      // -- left side      | right side -- //
      rho_l      = 1.0;      rho_r      = 0.125;
      pressure_l = 1.0;      pressure_r = 0.1;
      vx_l       = 0.0;      vx_r       = 0.0;
      break;

    case (2):
      rho_l      = 1.0;      rho_r      = 1.0;
      pressure_l = 0.4;      pressure_r = 0.4;
      vx_l       =-2.0;      vx_r       = 2.0;
      break;

    case (3):
      rho_l      = 1.0;      rho_r      = 1.0;
      pressure_l = 1000.;    pressure_r = 0.01;
      vx_l       = 0.0;      vx_r       = 0.0;
      break;

    case (4):
      rho_l      = 1.0;      rho_r      = 1.0;
      pressure_l = 0.01;     pressure_r = 100.;
      vx_l       = 0.0;      vx_r       = 0.0;
      break;

    case (5):
      rho_l      = 5.99924;  rho_r      = 5.99242;
      pressure_l = 460.894;  pressure_r = 46.0950;
      vx_l       = 19.5975;  vx_r       =-6.19633;
      break;

    default:
      if (rank == 0)
        cerr << "ERROR: invalid test (" << sodtest_num << ")." << endl;
      MPI_Finalize();
      exit(-1);
  }

  distance_r = length/2./(nparticles/9.);
  distance_l = length/2./(8.*(nparticles/9.));

  std::cout<<"distance L = "<<distance_l<<" R = "<<distance_r<<std::endl;

  // output file
  sprintf(output_filename,"%s.h5part",fileprefix);

}

//----------------------------------------------------------------------------//
int main(int argc, char * argv[]){
  using namespace std;
  using namespace simulation_params;

  // launch MPI
  int rank, size, provided;
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  // set simulation parameters
  set_default_param(rank,size);
  parse_command_line_options(rank,size,argc,argv);
  set_param(rank,size);

  // screen output
  if(rank==0){
    cout << "Sod test #" << sodtest_num << " in 1D:" << endl
         << " - number of particles: " << nparticles << endl
         << " - output file: " << output_filename << endl;
  }

  // allocate arrays

  // Position
  double* x = new double[nparticles]();
  double* y = new double[nparticles]();
  double* z = new double[nparticles]();
  // Velocity
  double* vx = new double[nparticles]();
  double* vy = new double[nparticles]();
  double* vz = new double[nparticles]();
  // Acceleration
  double* ax = new double[nparticles]();
  double* ay = new double[nparticles]();
  double* az = new double[nparticles]();
  // Smoothing length
  double* h = new double[nparticles]();
  // Density
  double* rho = new double[nparticles]();
  // Internal Energy
  double* u = new double[nparticles]();
  // Pressure
  double* P = new double[nparticles]();
  // Mass
  double* m = new double[nparticles]();
  // Id
  int64_t* id = new int64_t[nparticles]();
  // Timestep
  double* dt = new double[nparticles]();

  // Generate data
  // Find middle to switch m, u and rho
  //double middle = nparticles*ldistance/2.;
  // Find my first particle position
  double lposition = -length/2.;//ldistance*nparticles*rank;
  // Id of my first particle
  int64_t posid = 0;//nparticles*rank;

  // max. value for the speed of sound
  double cs = sqrt(localgamma*max(pressure_l/rho_l,pressure_r/rho_r));

  // The value for constant timestep
  double timestep = 0.07*ldistance/cs;


  for(int64_t part=0 ; part<nparticles ; ++part){
    id[part] = posid++;
    x[part] = lposition;

    double dist;

    if(lposition <= 0.){
      P[part] = pressure_l;
      rho[part] = rho_l;
      vx[part] = vx_l;
      dist = distance_l;
    }else{
      P[part] = pressure_r;
      rho[part] = rho_r;
      vx[part] = vx_r;
      dist = distance_r;
    }

    lposition += dist;

    // compute internal energy using gamma-law eos
    u[part] = P[part]/((localgamma-1.)*rho[part]);

    // particle masses and smoothing length
    // h = gamma(m/rho)
    m[part] = 1.;
    // Guess to start the search 
    //h[part] = localgamma*(m[part]/rho[part]);
    //m[part] = rho[part]*h[part]/localgamma;

  } // for part=0..nparticles

  bisection_smoothing_length(x,h,m,rho);

  std::cout<<"Position = ["<<-length/2.<<";"<<lposition<<"]"<<std::endl;

  // delete the output file if exists
  remove(output_filename);

  // Header data
  // the number of particles = nparticles
  Flecsi_Sim_IO::HDF5ParticleIO testDataSet;
  testDataSet.createDataset(output_filename,MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",nparticles);
  testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",1);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  //testDataSet.writeDatasetAttributeArray("name","string",simName);
  testDataSet.closeFile();

  testDataSet.openFile(MPI_COMM_WORLD);
  testDataSet.setTimeStep(0);

  Flecsi_Sim_IO::Variable _d1,_d2,_d3;

  _d1.createVariable("x",Flecsi_Sim_IO::point,"double",nparticles,x);
  _d2.createVariable("y",Flecsi_Sim_IO::point,"double",nparticles,y);
  _d3.createVariable("z",Flecsi_Sim_IO::point,"double",nparticles,z);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("vx",Flecsi_Sim_IO::point,"double",nparticles,vx);
  //_d2.createVariable("vy",Flecsi_Sim_IO::point,"double",nparticles,vy);
  //_d3.createVariable("vz",Flecsi_Sim_IO::point,"double",nparticles,vz);

  testDataSet.vars.push_back(_d1);
  //testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("h",Flecsi_Sim_IO::point,"double",nparticles,h);
  _d2.createVariable("rho",Flecsi_Sim_IO::point,"double",nparticles,rho);
  _d3.createVariable("u",Flecsi_Sim_IO::point,"double",nparticles,u);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("P",Flecsi_Sim_IO::point,"double",nparticles,P);
  _d2.createVariable("m",Flecsi_Sim_IO::point,"double",nparticles,m);
  _d3.createVariable("id",Flecsi_Sim_IO::point,"int64_t",nparticles,id);

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
