# XEOS library

The purpose of this project is to create a uniform user-friendly
C++ API to a variety of equations of state, from simple analytic 
ones such as a polytrope or an ideal gas, to more complex multi-
parametric EoS tables such as 
[CompOSE](https://compose.obspm.fr/).

Use case examples:

1. Set global physical unitrs to e.g. CGS
All physical quantities created hereafter will be in CGS
by default, unless specified otherwise.
```cpp
   PhysicalQuantity::SetGlobalUnits (PhUnits::CGS);
```


2. Create the simplest analytic EOS - polytropic.
Notice how you can specify polytropic constant in units which
are different from global. It is automatically converted to 
global units.
```cpp
  XeosPolytropic eos1 = new XeosPolytropic(
   4./3.,             // Gamma
   0.01,              // K in P = K*rho^Gamma 
   PhUnits::NUCLEAR); // units of K (different from global)
```

3. Create CompOSE tabulated EOS with three parameters:
(TODO: implement!!)
```cpp
  XeosCompose3D eos2 = new XeosCompose3D(
   "data/compose/eos.thermo",  // path to the data file
   {Pq::Pressure,              // primary variables
    Pq::ElectronFraction,
    Pq::Temperature,
    Pq::BaryonMassDensity},  
   PhUnits::NUCLEAR);  // CompOSE tables are in nuclear units
```

4. We can use xeos_base class to access both equations of 
state in a unified manner through e.g. a pointer:
```cpp  
  XeosBase *eosptr = { &eos1, &eos2 };
```

5. In this example, eos1 has overloaded `operator()`, which has two 
arguments: input and output. The first one is a const reference, 
while the second a pointer. The following code computes pressure
from density using eos1 object:
```cpp
  PqDensity  rho(12.0); // specify density in [g/cm3]
  PqPressure P;         // declare pressure variable [dynes]
  eos1(rho, &P);        // arg1: input, arg2: output.
```

6. The following code compute pressure at {rho, eps, Ye} for the 
CompOSE eos (TODO: not implemented yet):
```cpp
  PqSpecificInternalEnergy eps(0.01); // specify eps [erg/g]
  PqElectronFraction         Ye(0.5); // specify Ye [mol/g]
  eos2(rho,eps,Ye, P);  // args 1-3: input, arg 4: output
```

7. Class `PhysicalNVector<Kind,ScalarType>` is a template class for 
an N-dimensional array of physical quantities, handy for vector 
operations. E.g.: computing temperature array:
```cpp
  typedef 
    PhysicalNVector<Pq::Temperature,double> NvTemperature;

  const int Np = 1000;
  NvTemperature Temps(Np);
  NvDensity rhos(Np);
  NvSpecificInternalEnergy epss(Np);
  NvElectronFraction Yess(Np);
  NvPressure Ps(Np);

  for (i=0; i<Np; ++i) {
    rhos[i] = <...> // fill in arrays
    epss[i] = <...>
    Yess[i] = <...>
  }

  // now everything is ready for an EoS function call, 
  // which is just like before for the scalars:
  eos1(rhos, Ps);
  eos2(rhos, epss, Yess, Temps);
```

8. TODO
Compute multiple quantities: some complex equations of
state have interface which allows computation of several
physical quantities at once. 
```cpp
  NvEntropy Ents(Np);
  eos2 (3, (eos_in) {&rhos, &epss, &Yess},  // input arrays
           (eos_out){&Temps, &Ps, &Ents);  // output arrays
```


