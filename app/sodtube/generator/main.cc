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
  int64_t nparticles;        // Total number of particles
  double ldistance;          // Distance between the particles 
  double localgamma;         // polytropic index
  double smoothing_length;   // constant smoothing length

  // test conditions for two sides of the domain
  int    sodtest_num;            // which Sod test to generate
  double rho_1, rho_2;           // densities
  double pressure_1, pressure_2; // pressures
  double u_1, u_2;               // internal energies
  double m_1, m_2;               // particle masses

  // output filename
  const char* fileprefix = "hdf5_sodtube";
}

// 
// set default parameters 
//
void set_default_param(int rank) {
  using namespace simulation_params;

  // number of particles
  nparticles = 1000;

  // setup Sod test #1
  sodtest_num = 1;
  ldistance = 0.001;  // Distance between the particles 
  localgamma = 5./3.;
  rho_1 = 1;          rho_2 = 0.125;
  pressure_1 = 1;     pressure_2 = 0.1;
  u_1 = 2.5;          u_2 = 2;
  m_1 = 1.0e-4;       m_2 = 1.0e-5;
  smoothing_length = 1.0e-2;
}

// 
// parse command-line options 
//
int option_parser(int rank, int argc, char * argv[]) {
  using namespace std;
  for (int i=1; i<argc; ++i) 
    if (argv[i][0] == '-') 
      switch(argv[i][1]) {
      case 'h':
        if (rank == 0)
          cout << "Usage: ./" << argv[0] 
               << " -n <number of particles>" << endl;
        return (-1);
      } // switch argv[i][1]

  return 0;

}

int main(int argc, char * argv[]){
  using namespace simulation_params;

  int rank, size;
  int provided;  
  MPI_Init_thread(&argc,&argv,MPI_THREAD_MULTIPLE,&provided);
  assert(provided>=MPI_THREAD_MULTIPLE);
  MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  MPI_Comm_size(MPI_COMM_WORLD,&size);

  // TODO
  set_default_param(rank);

  // 2. parse options
  int parsed= option_parser(rank,argc,argv);
  if (parsed != 0) {
    MPI_Finalize();
    return (0);
  }

  if(argc!=2){
    if(rank==0){
      printf("./sodtube_generator [nParticles]\n");
      fprintf(stderr,"Running with default number of particles=1000\n");
    }
  }else{
    nparticles = atoll(argv[1]);
  }

  int64_t nparticlesproc = nparticles/size;
  if(rank==size-1){
    nparticlesproc = nparticles - nparticlesproc*(size-1);
  }

  if(rank==0){
    printf("Generating %ld particles\n",nparticles);
    printf("%ld particles per proc (last %ld)\n",nparticlesproc,
        nparticles-nparticlesproc*(size-1));
  }

  // Position
  double* x = new double[nparticlesproc]();
  double* y = new double[nparticlesproc]();
  double* z = new double[nparticlesproc]();
  // Velocity
  double* vx = new double[nparticlesproc]();
  double* vy = new double[nparticlesproc]();
  double* vz = new double[nparticlesproc]();
  // Acceleration
  double* ax = new double[nparticlesproc]();
  double* ay = new double[nparticlesproc]();
  double* az = new double[nparticlesproc]();
  // Smoothing length 
  double* h = new double[nparticlesproc]();
  // Density 
  double* rho = new double[nparticlesproc]();
  // Internal Energy 
  double* u = new double[nparticlesproc]();
  // Pressure
  double* P = new double[nparticlesproc]();
  // Mass
  double* m = new double[nparticlesproc]();
  // Id
  int64_t* id = new int64_t[nparticlesproc]();
  // Timestep 
  double* dt = new double[nparticlesproc]();
  
  // Generate data
  // Find middle to switch m, u and rho  
  double middle = nparticles*ldistance/2.;
  // Find my first particle position 
  double lposition = ldistance*nparticlesproc*rank;
  // Id of my first particle 
  int64_t posid = nparticlesproc*rank;

  // Header data 
  // the number of particles = nparticles 
  // The value for constant timestep 
  double timestep = 0.001;
  int dimension = 1;
  
  
  for(int64_t part=0; part<nparticlesproc; ++part){
    x[part] = lposition;

    if(x[part] > middle){
      P[part] = pressure_2;
      rho[part] = rho_2; 
      u[part] = u_2;
      //m[part] = m_2;
    }else{
      P[part] = pressure_1;
      rho[part] = rho_1;
      u[part] = u_1;
      //m[part] = m_1;
    }

    m[part] = rho[part]*middle/(nparticles/2.);

    //m[part] = 0.;
    // Y and Z not used 
    // VX, VY, VZ and AX, AY, AZ stay to 0
    h[part] = smoothing_length;
    
    // P stay to 0
    id[part] = posid++; 
    // Move to the next particle 
    lposition += ldistance;
    //std::cout<<x[part]<<": "<<h[part]<<std::endl;
  }

  char filename[128];
  //sprintf(filename,"%s_%d.h5part",fileprefix,nparticles);
  sprintf(filename,"%s.h5part",fileprefix);

  // Destroy the file if exists 
  remove(filename);

  Flecsi_Sim_IO::HDF5ParticleIO testDataSet; 
  testDataSet.createDataset(filename,MPI_COMM_WORLD);

  // add the global attributes
  testDataSet.writeDatasetAttribute("nparticles","int64_t",nparticles);
  //testDataSet.writeDatasetAttribute("timestep","double",timestep);
  testDataSet.writeDatasetAttribute("dimension","int32_t",1);
  testDataSet.writeDatasetAttribute("use_fixed_timestep","int32_t",1);

  //testDataSet.writeDatasetAttributeArray("name","string",simName);
  testDataSet.closeFile();

  testDataSet.openFile(MPI_COMM_WORLD);
  testDataSet.setTimeStep(0);

  Flecsi_Sim_IO::Variable _d1,_d2,_d3;

  _d1.createVariable("x",Flecsi_Sim_IO::point,"double",nparticlesproc,x);
  _d2.createVariable("y",Flecsi_Sim_IO::point,"double",nparticlesproc,y);
  _d3.createVariable("z",Flecsi_Sim_IO::point,"double",nparticlesproc,z);

  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  //_d1.createVariable("vx",Flecsi_Sim_IO::point,"double",nparticlesproc,vx);
  //_d2.createVariable("vy",Flecsi_Sim_IO::point,"double",nparticlesproc,vy);
  //_d3.createVariable("vz",Flecsi_Sim_IO::point,"double",nparticlesproc,vz);

  //testDataSet.vars.push_back(_d1);
  //testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  //testDataSet.writeVariables();

  //_d1.createVariable("ax",Flecsi_Sim_IO::point,"double",nparticlesproc,ax);
  //_d2.createVariable("ay",Flecsi_Sim_IO::point,"double",nparticlesproc,ay);
  //_d3.createVariable("az",Flecsi_Sim_IO::point,"double",nparticlesproc,az);

  //testDataSet.vars.push_back(_d1);
  //testDataSet.vars.push_back(_d2);
  //testDataSet.vars.push_back(_d3);

  //testDataSet.writeVariables();


  _d1.createVariable("h",Flecsi_Sim_IO::point,"double",nparticlesproc,h);
  _d2.createVariable("rho",Flecsi_Sim_IO::point,"double",nparticlesproc,rho);
  _d3.createVariable("u",Flecsi_Sim_IO::point,"double",nparticlesproc,u);
  
  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  _d1.createVariable("P",Flecsi_Sim_IO::point,"double",nparticlesproc,P);
  _d2.createVariable("m",Flecsi_Sim_IO::point,"double",nparticlesproc,m);
  _d3.createVariable("id",Flecsi_Sim_IO::point,"int64_t",nparticlesproc,id);
  
  testDataSet.vars.push_back(_d1);
  testDataSet.vars.push_back(_d2);
  testDataSet.vars.push_back(_d3);

  testDataSet.writeVariables();

  testDataSet.closeFile(); 
  
  delete[] x,y,z,vx,vy,vz,ax,ay,az,h,rho,u,P,m,id,dt;
 
  MPI_Finalize();
  return 0;
}
