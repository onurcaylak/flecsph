/*~-------------------------------------------------------------------~*
 * Copyright (c) 2018 Los Alamos National Security, LLC
 * All rights reserved.
 *~-------------------------------------------------------------------~*/

/**
 * @file xeoscoldn.h
 * @author Oleg Korobkin
 * @date 2018-01-05
 * @brief XeosColdNuclear: tabulated 1D equation of state 
 *        for cold nuclear matter
 */
#ifndef XEOS_XEOSCOLDN_H_
#define XEOS_XEOSCOLDN_H_

#include <math.h>
#include "physquan.h"
#include "xeosbase.h"
#include "filereader.h"

namespace xeos {

/**
 * Single-parameter tabulated EoS for cold nuclear matter
 * TODO: implement details, give usage example here
 */
template <class FReader>
class XeosColdNuclear : public XeosTabulated {
 private:

  std::vector<Pq> primary_vars;
  std::vector<Pq> additional_vars;
  FReader freader;

 protected:
  // units in which the tables are given
  PhUnits eos_units = PhUnits::CGS;

  using XeosTabulated::data_path;
  using XeosTabulated::format;

 public:

  // two different formats for now..
  static const format_type
    rns_4col  = 0x4C01,
    rns_16col = 0x16C01;

  // default constructor
  XeosColdNuclear()
      : XeosColdNuclear( "",XeosColdNuclear::rns_4col,PhUnits::CGS) {}

  // workhorse constructor
  XeosColdNuclear (const std::string _path,
      const format_type _fmt,
      const PhUnits _u);

  // report EoS units
  PhUnits GetUnits() const { return eos_units; }

  virtual void ReadEosTables() override;

  virtual std::vector<Pq>
  GetPrimaryQuantities() const override { return primary_vars; }

  virtual std::vector<Pq>
  GetAdditionalQuantities() const override { return additional_vars; }

  virtual void
  operator() (const int num_out, eos_in in_array[],
      eos_out out_array[]) override;

  virtual void
  DfDx (const int num_out, const Pq dX, eos_in  in_array[],
      eos_out dF_array[]) override;

  virtual bool
  ConsistencyCheck() const override;

}; // class XeosColdNuclear


//
// XeosColdNuclear: workhorse constructor definition
//
template <class FReader>
XeosColdNuclear<FReader>::XeosColdNuclear (
    const std::string _path,
    const format_type _fmt,
    const PhUnits _u) {

  XeosTabulated::data_path = _path;
  XeosTabulated::format = _fmt;
  eos_units = PhUnits::CGS;

  // setup primary and secondary variables
  switch(_fmt) {
    case rns_4col:
      primary_vars.push_back(Pq::RestMassDensity);
      primary_vars.push_back(Pq::Pressure);
      primary_vars.push_back(Pq::Enthalpy);
      primary_vars.push_back(Pq::BaryonNumberDensity);
      break;

    case rns_16col:
    default:
      assert(false); // format not implemented
  }
}

//
// File reader (TODO: implement with templates)
//
template <class FReader>
void XeosColdNuclear<FReader>::ReadEosTables() {
  /*      +
     TODO
   +      */
}

//
// General EoS function
//
template <class FReader>
void XeosColdNuclear<FReader>::operator() (const int num_out,
    eos_in in_array[],
    eos_out out_array[]) {
  /*      +
     TODO
   +      */
}

//
// General derivative(s), (dF/dX)_const
//
template <class FReader>
void XeosColdNuclear<FReader>::DfDx (const int num_out, const Pq dX,
      eos_in  in_array[],
      eos_out dF_array[]) {
  /*      +
     TODO
   +      */
}

template <class FReader>
bool XeosColdNuclear<FReader>::ConsistencyCheck() const {
  /*      +
     TODO
   +      */
  return false;
}

} // namespace xeos
#endif // XEOS_XEOSCOLDN_H_

#if 0
#include <iostream>

int main() {
  using namespace std;
  using namespace xeos;
  double fields[20];
  AsciiFileReader A("../data/rnsid/sfho_0.1MeV_beta.txt");
  cout << "A.GetFilename() -> " << A.GetFilename() << endl;
  int file_len = A.NumLines();
  cout << "A.NumLines() ->    " << file_len << endl;
  A.Open();
  int header_len = A.SkipHashHeader();
  cout << "A.SkipHashHeader() -> " << header_len << endl;
  int n_fields;
  cout << scientific << setprecision(6);
  A.Rewind(); // rewind;
  //A.SkipHashHeader();
  A.SkipHeader(header_len);
  for (int l=0; l < file_len - header_len; ++l) {
    n_fields = A.ReadFields(20,fields);
    if (n_fields<16)
      break;
    //cout << setw(4) << l << ":";
    cout << fields[0];
    for (int i=1;i<n_fields;++i)
      cout << " " << fields[i];
    cout << endl;
  }
  A.Close();
}
#endif
