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

#ifndef _physics_kernel_h_
#define _physics_kernel_h_

#include <vector>

#include "tree.h"

namespace kernel{


  // Standard spline kernel
  // for 1-2-3D
  //static 
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

  // Standard gradient of spline kernel
  // for 1-2-3D  
  //static
  point_t 
  cubic_spline_gradKernel(
      point_t vecP, 
      double h)
  {
    // Default 1D case
    double sigma = 2./(3.*h*h);
    if(gdimension == 2){
      sigma = 10./(7.*M_PI*h*h*h);
    }
    // Compute distance particles 
    double r = 0;
    for(size_t i=0;i<gdimension;++i){
      r += vecP[i]*vecP[i];
    }
    r = sqrt(r);

    // Normalize vector 
    point_t eab = vecP / r;

    double rh = r/h;

    point_t result{};
    if (0.0 <= rh && rh <= 1.0){
      result = sigma*eab;
      result *= -3.0*rh + 9./4.*rh*rh;
    }else if(1.0 < rh && rh <= 2.0){
      result = sigma*eab;
      result *= -.75*(2-rh)*(2-rh);
    }
    return result;
  } // gradKernel 

}; // kernel

#endif // _physics_kernel_h_
