/*~--------------------------------------------------------------------------~*
 * Copyright (c) 2017 Los Alamos National Security, LLC
 * All rights reserved.
 *~--------------------------------------------------------------------------~*/

 /*~--------------------------------------------------------------------------~*
 * 
 * /@@@@@@@@  @@           @@@@@@   @@@@@@@@ @@@@@@@  @@      @@
 * /@@/////  /@@          @@////@@ @@////// /@@////@@/@@     /@@
 * /@@       /@@  @@@@@  @@    // /@@       /@@   /@@/@@     /@@
 * /@@@@@@@  /@@ @@///@@/@@       /@@@@@@@@@/@@@@@@@ /@@@@@@@@@@
 * /@@////   /@@/@@@@@@@/@@       ////////@@/@@////  /@@//////@@
 * /@@       /@@/@@//// //@@    @@       /@@/@@      /@@     /@@
 * /@@       @@@//@@@@@@ //@@@@@@  @@@@@@@@ /@@      /@@     /@@
 * //       ///  //////   //////  ////////  //       //      //  
 *
 *~--------------------------------------------------------------------------~*/

/**
 * @file physics.h
 * @author Julien Loiseau
 * @date April 2017
 * @brief Basic physics implementation
 */

#ifndef _physics_physics_h_
#define _physics_physics_h_

#include <vector>

#include "kernel.h"
#include "tree.h"

namespace physics{

  bool fixed_timestep = true;
  bool do_boundaries = false;
  bool stop_boundaries = false; 
  bool reflect_boundaries = false;
  point_t max_boundary = {};
  point_t min_boundary = {};
  double dt = 0.0;
  double alpha = 2; 
  double beta = 2; // 1; 
  double gamma = 2.0; // 1.4; 
  double K = 1;
  double epsilon = 1;
  double g_strength = 1; 
  double damp = 1;
  double totaltime = 0.0;
  double MAC = 0.;
  double eta = 0.01;
  int64_t nparticles = 0;

  int verlet_cstep = 10.;
  int verlet_current = 0.;
  bool do_verlet_cor = false; 

  // Default configuration for kernel
  int kernel_choice = 0;
  auto kernel = kernel::cubic_spline_kernel;
  auto gradKernel = kernel::cubic_spline_gradKernel;

  // Compute density based on neighbors
  // The formula used here is:
  // rho_i = \sum_j m_j*cubic\_spline\_kernel\_3D_ij
  //static
  void 
  compute_density(
      body_holder* srch, 
      std::vector<body_holder*>& nbsh)
  {
    body* source = srch->getBody();
    double density = 0.;
    mpi_assert(nbsh.size()>0);

    // Reset the smoothing length if out the boundaries
    for(auto nbh : nbsh){
      body* nb = nbh->getBody();
      double dist = flecsi::distance(source->getPosition(),nb->getPosition());
      mpi_assert(dist>=0.0);
      density += nb->getMass()*
        kernel(
        dist,
        1./2.*(source->getSmoothinglength()+nb->getSmoothinglength()));
    } // for
    mpi_assert(density>0);
    source->setDensity(density);
  } // compute_density

  // Computre pressure on a body 
  // This function does not need the neighbors
  // Formula is:
  // P_i = u_i*rho_i^{\Gamma}
  //static
  void 
  compute_pressure_sodtube(
      body_holder* srch)
  { 
    body* source = srch->getBody();
    double pressure = (gamma-1.0)*
      source->getDensity()*source->getInternalenergy();
    source->setPressure(pressure);
  } // compute_pressure


  void 
  init_energy_velocity(
    body_holder* srch)
  {
    body* source = srch->getBody();
    source->setVelocityTmp(source->getVelocity());
    source->setInternalenergyTmp(source->getInternalenergy());
  }

    // Computre pressure on a body 
  // This function does not need the neighbors
  // Formula is:
  // P_i = u_i*rho_i^{\Gamma}
  //static
  void 
  compute_pressure_sedov(
      body_holder* srch)
  { 
    body* source = srch->getBody();
    //double pressure = (gamma-1.0)*
    //  source->getDensity()*source->getInternalenergy();
    double A = (gamma-1)*source->getInternalenergy()/
      (pow(source->getDensity(),gamma-1));
    //double A = 0.6366;
    double pressure = A*
      pow(source->getDensity(),gamma);
    //assert(pressure>=0);
    source->setPressure(pressure);
  } // compute_pressure


  //For zero temperature white dwarf EOS
  void 
  compute_pressure_wd(
      body_holder* srch)
  { 
    body* source = srch->getBody();
    double A_dwd = 6.00288e22;
    double B_dwd = 9.81011e5;

    double x_dwd = pow((source->getDensity())/B_dwd,1.0/3.0);
    double pressure = A_dwd*(x_dwd*(2.0*x_dwd*x_dwd-3.0)*
 		      pow(x_dwd*x_dwd+1.0,1.0/2.0)+3.0*asinh(x_dwd));
    source->setPressure(pressure);
  } // compute_pressure_wd

  // Compute the sound speed on a body 
  // This function does not need the neighbors
  // Formula is:
  // CS_i = (P_i/rho_i)^{1/2}
  //static
  void 
  compute_soundspeed(
      body_holder* srch)
  {
    body* source = srch->getBody();
    double soundspeed = sqrt(gamma*
        source->getPressure()/source->getDensity());
    source->setSoundspeed(soundspeed);
  } // computeSoundspeed

  void 
  compute_smoothing_length(
    body_holder* srch)
  {
    body* source = srch->getBody(); 
    double smoothing_length = gamma*source->getMass()/source->getDensity();
    source->setSmoothinglength(smoothing_length);
  }

  void 
  compute_density_pressure_soundspeed_sedov(
    body_holder* srch, 
    std::vector<body_holder*>& nbsh)
  {
    compute_density(srch,nbsh);
    compute_pressure_sedov(srch);
    compute_soundspeed(srch); 
  }

  void 
  compute_density_pressure_soundspeed_sodtube(
    body_holder* srch, 
    std::vector<body_holder*>& nbsh)
  {
    if(do_boundaries && stop_boundaries){
      if(srch->getPosition() < min_boundary || 
          srch->getPosition() > max_boundary){
        return;
      }
    }
    compute_density(srch,nbsh);
    compute_pressure_sodtube(srch);
    compute_soundspeed(srch);
  }

  // mu used in artificial viscosity 
  // Formula is:
  // h_ij*(v_i-v_j).(p_i-p_j)/(dist*dist + h_ij+epsi)
  // with h_ij = (h_i+h_j)/2
  //static
  double 
  mu(
      body* source, 
      body* nb,
      double gamma)
  {  
    double result = 0.0;
    double h_ij = (1./2.)*
      (source->getSmoothinglength()+nb->getSmoothinglength()); 
    space_vector_t vecVelocity = flecsi::point_to_vector(
        source->getVelocity() - nb->getVelocity());
    space_vector_t vecPosition = flecsi::point_to_vector(
        source->getPosition() - nb->getPosition());
    double dotproduct = flecsi::dot(vecVelocity,vecPosition);
    if(dotproduct >= 0.0)
      return result;
    // Should add norm to space_vector
    double dist = flecsi::distance(source->getPosition(),nb->getPosition());
    result = h_ij * dotproduct / (dist*dist + 0.01*h_ij*h_ij);
    //result = h_ij * dotproduct / (dist*dist+eta*eta);
    //  res = dotRes / (h_ij * (((dist*dist)/(h_ij*h_ij))+visc_eta*visc_eta));  

    mpi_assert(result < 0.0);
    return result; 
  } // mu

  // Compute the hydro force with the artificial viscosity 
  // Formula is:
  // f_{hydro}_i = m_j * (\sum_j P_i/rho_i^2 + P_j/rho_j^2 + Pi_{ij})*
  // gradKernel_{ij}
  // with Pi_{ij} the artificial viscosity 
  //static
  void 
  compute_hydroacc_internalenergy(
    body_holder* srch, 
    std::vector<body_holder*>& ngbsh)
  { 
    body* source = srch->getBody();

    // For Verlet 
    source->setVelocityNM1(source->getVelocityTmp());
    source->setVelocityTmp(source->getVelocity());
    source->setInternalenergyNM1(source->getInternalenergyTmp());
    source->setInternalenergyTmp(source->getInternalenergy());

    point_t hydro = {};
    point_t velCor; 
    double dudt_pressure = 0.;
    double dudt_viscosity = 0.;
    double dudt = 0.;

    for(auto nbh : ngbsh){ 
      body* nb = nbh->getBody();

      // Ignore the particle itself
      if(nb->getId() == source->getId()){
        continue;
      }

      // Artificial viscosity
      double density_ij = (1./2.)*(source->getDensity()+nb->getDensity());
      double soundspeed_ij = (1./2.)*
        (source->getSoundspeed()+nb->getSoundspeed());
      double mu_ij = mu(source,nb,epsilon);
      double viscosity = (-alpha*mu_ij*soundspeed_ij+beta*mu_ij*mu_ij)
        /density_ij;
      mpi_assert(viscosity>=0.0);

      double dist = flecsi::distance(source->getPosition(),nb->getPosition());
      double wab =kernel(dist,
          1./2.*(source->getSmoothinglength()+nb->getSmoothinglength()));
        
      // Gradient
      space_vector_t vecVelocity = flecsi::point_to_vector(
          source->getVelocity()-nb->getVelocity());
      point_t vecPosition = source->getPosition()-nb->getPosition();
      double pressureDensity = 
          source->getPressure()/(source->getDensity()*source->getDensity())
          + nb->getPressure()/(nb->getDensity()*nb->getDensity());
      
      //point_t sourcekernelgradient = gradKernel(
      //    vecPosition,source->getSmoothinglength());
      //point_t nbkernelgradient = gradKernel(
      //    vecPosition,nb->getSmoothinglength());
      //point_t resultkernelgradient = (1./2.)*
      //  (sourcekernelgradient+nbkernelgradient);
      point_t resultkernelgradient = gradKernel(
          vecPosition,
          1./2.*(source->getSmoothinglength()+nb->getSmoothinglength()));
      
      hydro = nb->getMass()*(pressureDensity+viscosity)
        *resultkernelgradient;

      dudt_pressure += nb->getMass()*
        flecsi::dot(vecVelocity,flecsi::point_to_vector(resultkernelgradient));
      dudt_viscosity += nb->getMass()*
        viscosity*
        flecsi::dot(vecVelocity,flecsi::point_to_vector(resultkernelgradient));

        velCor = velCor - wab * (source->getVelocity()-nb->getVelocity()) * nb->getMass() / density_ij; 
    }

    dudt = source->getPressure()/
      (source->getDensity()*source->getDensity())*dudt_pressure+
      1./2.*dudt_viscosity;
    
    source->setDudt(dudt);
    source->setAcceleration(-1.*hydro);
    source->setVelocityCor(source->getVelocity()+epsilon*velCor);

  } // compute_hydro_acceleration

  // Apply boundaries 
  bool
  compute_boundaries(
      body_holder* srch)
  {
    body* source = srch->getBody();
    point_t velocity = source->getVelocity();
    point_t position = source->getPosition();
    point_t velocityHalf = source->getVelocityhalf();

    bool considered = false;

    if(stop_boundaries){
      if(position < min_boundary || 
          position > max_boundary){
        velocity = point_t{};
        velocityHalf = point_t{};
        considered = true;
      }
    }else if(reflect_boundaries){
      for(size_t dim=0;dim < gdimension ; ++dim){
        if(position[dim] < min_boundary[dim] || 
            position[dim] > max_boundary[dim]){
          double barrier = max_boundary[dim];
          if(position[dim] < min_boundary[dim]){
            barrier = min_boundary[dim];
          }
          // Here just invert the velocity vector and velocityHalf 
          double tbounce = (position[dim]-barrier)/velocity[dim];
          position -= velocity*(1-damp)*tbounce;

          position[dim] = 2*barrier-position[dim];
          velocity[dim] = -velocity[dim];
          velocityHalf[dim] = -velocityHalf[dim];

          velocity *= damp;
          velocityHalf *= damp;
          considered = true;
        }
      }
    }
    if(considered){
      source->setPosition(position);
      source->setVelocity(velocity);
      source->setVelocityhalf(velocityHalf);
    }
    return considered;
  }

    void 
  leapfrog_integration_first_step(
      body_holder* srch)
  {
    body* source = srch->getBody();

    point_t velocityHalf = source->getVelocity() + 
        dt/2.*source->getAcceleration();
    point_t position = source->getPosition()+velocityHalf*dt;
    point_t velocity = 1./2.*(source->getVelocityhalf()+velocityHalf);

    compute_smoothing_length(srch); 

    // Internal energy 
    source->setInternalenergy(source->getInternalenergy()+
        dt*source->getDudt());


    if(do_boundaries){
      if(physics::compute_boundaries(srch)){
        return;
      }
    }
 
    source->setVelocityhalf(velocityHalf);
    source->setVelocity(velocity);
    source->setPosition(position);
    
    mpi_assert(!std::isnan(position[0])); 
  }

  void 
  leapfrog_integration(
      body_holder* srch)
  {
    body* source = srch->getBody();

    point_t velocityHalf = source->getVelocityhalf() + 
        dt/2.*source->getAcceleration();
    point_t position = source->getPosition()+velocityHalf*dt;
    point_t velocity = 1./2.*(source->getVelocityhalf()+velocityHalf);

    compute_smoothing_length(srch); 

    // Internal energy 
    source->setInternalenergy(source->getInternalenergy()+
        dt*source->getDudt());


    if(do_boundaries){
      if(physics::compute_boundaries(srch)){
        return;
      }
    }
    
    source->setVelocityhalf(velocityHalf);
    source->setVelocity(velocity);
    source->setPosition(position);
    
    mpi_assert(!std::isnan(position[0])); 
  }

  void 
  compute_dt(
      body_holder* srch,
      std::vector<body_holder*>& ngbhs)
  {
    body* source = srch->getBody();
    
    // First compute dt based on acceleration norm 
    double accelNorm = 0.0;
    for(size_t i=0;i<gdimension;++i){
      accelNorm += source->getAcceleration()[i]*source->getAcceleration()[i];
    }
    accelNorm = sqrt(accelNorm);
    double dt1 = pow(source->getSmoothinglength()/accelNorm,1.0/2.0);
    //std::cout<<"dt1 = "<<dt1<<std::endl;
  
    // Second based on max mu 
    double max_mu_ij = -9999999;
    for(auto nbh: ngbhs){
      body* nb = nbh->getBody(); 
      double local_mu = mu(source,nb,gamma);
      max_mu_ij = std::max(max_mu_ij,local_mu);
    }
    double dt2 = source->getSmoothinglength()/
        (source->getSoundspeed()+
         .6*alpha*source->getSoundspeed()+
         .6*beta*max_mu_ij);
    
    double min = std::min(dt1,dt2)*0.25;

    #pragma omp critical 
    physics::dt = std::min(dt,min);
    //return std::min(physics::dt,dt1); 
  }

  void 
  verlet_integration(
      body_holder* srch)
  {

    verlet_current++;
    do_verlet_cor = false;
    if(verlet_current % verlet_cstep == 0){
      do_verlet_cor = true;
      verlet_current = 0;
    }

    body* source = srch->getBody();

    point_t position = source->getPosition();
    point_t velocity = {};
    double internalenergy = 0.;

    if(!do_verlet_cor){
      position += source->getVelocityCor()*dt+
        source->getAcceleration()*dt*dt*0.5;
      velocity = source->getVelocityNM1()+source->getAcceleration()*2*dt;
      internalenergy = source->getInternalenergyNM1() + 
        source->getDudt()*2.*dt;
    }else{
      position += source->getVelocityCor()*dt+
        source->getAcceleration()*dt*dt*0.5;
      velocity = source->getVelocity()+source->getAcceleration()*dt;
      internalenergy = source->getInternalenergy() + source->getDudt() * dt;
    }
    
    if(do_boundaries){
      if(physics::compute_boundaries(srch)){
        return;
      }
    }

    source->setVelocity(velocity);
    source->setPosition(position);
    source->setInternalenergy(internalenergy);

    compute_smoothing_length(srch); 


  } // verlet_integration_EOS


}; // physics

#endif // _physics_physics_h_
