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
 * @file lattice_type.h
 * @author Alexander Kaltenborn
 * @date July 2018
 * @brief function to choose lattice type and populate initial conditions
 *
 * Module for determining lattice types for different tests. Making it modular
 * for better readability of code and more generic applicability to different
 * tests.
 *
 * Having dimension set in the user.h file, this function will populate either
 * a rectangular lattice, a triangular lattice, an FCC lattice, or a HCP lattice.
 *
 * Inputs:
 *  lattice_type - int value, corresponding to 0:rect, 1:HCP, 2:FCC,
 *  domain_type  - int value, 0:cube, 1:sphere, 2:full box (? what did oleg want different between 0 and 3)
 *  bbox_min     - bottom corner location for smallest cube covering the domain
 *  bbox_max     - top corner location for smallest cube covering the domain
 *  count_only   - boolean for determing particle number or writing positions
 *
 *
 * TODO:
 */

 #include <stdlib.h>
 #include "mpi.h"
 #include "user.h"
 #include "params.h"

bool in_domain(
    const double x,
    const double y,
    const double z,
    const double x0,
    const double y0,
    const double z0,
    const double r,
    const int domain_type)
{
  if(domain_type==0){
    if(std::abs(x-x0)<=r && std::abs(y-y0)<=r && std::abs(z-z0)<=r){
      return true;
    } else{
      return false;
    }
  } else if(domain_type==1){
    return (x-x0)*(x-x0)+(y-y0)*(y-y0)+(z-z0)*(z-z0)<r*r;
  } else if(domain_type==2){
    return true;
  }
}
int64_t
generate_lattice(
    const int lattice_type,
    const int domain_type,
    const point_t& bbox_min,
    const point_t& bbox_max,
    bool count_only,
    double x[] = NULL,
    double y[] = NULL,
    double z[] = NULL)
 {
   using namespace param;

   // Determine radius of the domain as a radius of inscribed sphere
   double radius = 0.0;
   if(domain_type==1 || domain_type==0) {
      radius = bbox_max[0] - bbox_min[0];
      for(int j=1;j<gdimension;j++) 
        if (radius > bbox_max[j] - bbox_min[j])
          radius = bbox_max[j] - bbox_min[j];
   } 
   radius = 0.5*radius;

   ////double radius = 0;
   ////if(domain_type==1 || domain_type==0) radius = std::min_element(std::begin(bbox_max),std::end(bbox_max));

   // Central coordinates: in most cases this should be centered at (0,0,0)
   double x_c = (bbox_max[0]+bbox_min[0])/2.;
   double y_c = 0;
   double z_c = 0;
   if(gdimension > 1) y_c = (bbox_max[1]+bbox_min[1])/2.;
   if(gdimension > 2) z_c = (bbox_max[2]+bbox_min[2])/2.;

   // Setting the coordinates to start building the lattice
   double x_topproc = bbox_min[0];
   double y_topproc = 0;
   double z_topproc = 0;
   if(gdimension > 1) y_topproc = bbox_min[1];
   if(gdimension > 2) z_topproc = bbox_min[2];

   int64_t Nx = int((bbox_max[0]-bbox_max[0])/sph_separation);
   int64_t Ny = 1;
   int64_t Nz = 1;
   if(gdimension>1) Ny = int((bbox_max[1]-bbox_max[1])/sph_separation);
   if(gdimension>2) Nz = int((bbox_max[2]-bbox_max[2])/sph_separation);

   // True number of particles to be determined and returned
   int64_t tparticles = 0;
   // Id of my first particle
   int64_t posid = 0;
   double x_position;
   double y_position;
   double z_position;

   // The loop for lattice_type==0
   if(lattice_type==0){
     for(int i=0;i<Nz;i++){
       for(int j=0;j<Ny;j++){
         for(int k=0;k<Nx;k++){
            //rectangular lattice in 2 or 3D
           z_position = z_topproc + i*sph_separation;
           y_position = y_topproc + j*sph_separation;
           x_position = x_topproc + k*sph_separation;
           if(in_domain(x_position,y_position,z_position,x_c,y_c,z_c,radius,domain_type)){
             tparticles++;
             if(!count_only){
               x[posid] = x_position;
               y[posid] = y_position;
               z[posid] = z_position;
             }
             posid++;
           }
         }
       }
     }
   } else if(lattice_type==1){//hcp lattice in 2 or 3D (in 2D it is just triangular)
     for(int i=0;i<Nz;i++){
       for(int j=0;j<Ny;j++){
         for(int k=0;k<Nx;k++){
            //hcp lattice in 2 or 3D
            z_position = z_topproc + i*sqrt(2./3.)*sph_separation;
            if(i%2==1){
              y_position = y_topproc - sqrt(3)*sph_separation/6. + j*sqrt(3.)/2.*sph_separation;
              if(j%2==1){
                x_position = x_topproc + k*sph_separation;
              } else {
                x_position = x_topproc - sph_separation/2. + k*sph_separation;
              }
            } else {
              y_position = y_topproc + j*sqrt(3.)/2.*sph_separation;
              if(j%2==1){
                x_position = x_topproc - sph_separation/2. + k*sph_separation;
              } else {
                x_position = x_topproc + k*sph_separation;
              }
            }
           if(in_domain(x_position,y_position,z_position,x_c,y_c,z_c,radius,domain_type)){
             tparticles++;
             if(!count_only){
               x[posid] = x_position;
               y[posid] = y_position;
               z[posid] = z_position;
             }
             posid++;
           }
         }
       }
     }
   } else if(lattice_type==2){//fcc lattice in 2 or 3D (in 2D it is just triangular)
     for(int i=0;i<Nz;i++){
       for(int j=0;j<Ny;j++){
         for(int k=0;k<Nx;k++){
            //fcc lattice in 2 or 3D
            z_position = z_topproc + i*sqrt(2./3.)*sph_separation;
            if(i%3==0){
              y_position = y_topproc + j*sqrt(3.)/2.*sph_separation;
              if(j%2==1){
                x_position = x_topproc - sph_separation/2. + k*sph_separation;
              } else {
                x_position = x_topproc + k*sph_separation;
              }
            } else if(i%3==1) {
              y_position = y_topproc - sqrt(3)*sph_separation/6. + j*sqrt(3.)/2.*sph_separation;
              if(j%2==1){
                x_position = x_topproc + k*sph_separation;
              } else {
                x_position = x_topproc - sph_separation/2. + k*sph_separation;
              }
            } else if(i%3==2){
              y_position = y_topproc + sqrt(3)*sph_separation/6. + j*sqrt(3.)/2.*sph_separation;
              if(j%2==1){
                x_position = x_topproc + k*sph_separation;
              } else {
                x_position = x_topproc - sph_separation/2. + k*sph_separation;
              }
            }
           if(in_domain(x_position,y_position,z_position,x_c,y_c,z_c,radius,domain_type)){
             tparticles++;
             if(!count_only){
               x[posid] = x_position;
               y[posid] = y_position;
               z[posid] = z_position;
             }
             posid++;
           }
         }
       }
     }
   }

   return tparticles;
 }

