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

/*
 * enum for format types
 */
enum class XeosColdNuclearFormat {
  kFourColumnsCGS,
  kSixteenColumnsNuclear
};


/**
 * Single-parameter tabulated EoS for cold nuclear matter
 * TODO: implement details, give usage example here
 */
template <class FR>
class XeosColdNuclear : public XeosTabulated {
 public:
  typedef enum XeosColdNuclearFormat format_type;

  // default constructor
  XeosColdNuclear<FR>()
      : XeosColdNuclear<FR>( "", 
        format_type::kFourColumnsCGS,
        PhysicalQuantity::GetGlobalUnits()) {}

  // constructor which omits units 
  XeosColdNuclear<FR> (const std::string _path, const format_type _fmt)
      : XeosColdNuclear<FR>( _path, _fmt,
        PhysicalQuantity::GetGlobalUnits()) {}

  // workhorse constructor declaration
  XeosColdNuclear<FR> (const std::string _path, const format_type _fmt,
      const PhUnits _u);

  // report EoS units
  PhUnits GetUnits() const { return eos_units; }
  format_type GetFormat() const { return format; }

  virtual void ReadEosTables() override;

  virtual double
  EosTableGridPoint(const Pq var, const int i) const override;

  virtual const int
  EosTableGridSize(const Pq) const override { return num_records; }

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

 protected:
  // units in which to keep the tables
  PhUnits eos_units = PhysicalQuantity::GetGlobalUnits();

  using XeosTabulated::data_path;

 private:

  std::vector<Pq> primary_vars;
  std::vector<Pq> additional_vars;
  FR freader;
  int num_records;
  int num_fields;
  format_type format;

  NvPressure pres_table {PhUnits::CGS};
  NvDensity  rho_table  {PhUnits::CGS};
  NvEnthalpy enth_table {PhUnits::CGS};
  NvSpecificInternalEnergy eps_table {PhUnits::CGS};
  NvBaryonNumberDensity    nbar_table {PhUnits::CGS};
}; // class XeosColdNuclear


//
// XeosColdNuclear: workhorse constructor definition
//
template <class FR>
XeosColdNuclear<FR>::XeosColdNuclear (
    const std::string _path,
    const format_type _fmt,
    const PhUnits _u) {

  XeosTabulated::data_path = _path;
  format = _fmt;
  eos_units = _u;
  freader.SetFilename(_path);

  // setup primary and secondary variables
  switch(_fmt) {
    case format_type::kFourColumnsCGS:
      primary_vars.push_back(Pq::RestMassDensity);
      primary_vars.push_back(Pq::Pressure);
      primary_vars.push_back(Pq::Enthalpy);
      primary_vars.push_back(Pq::BaryonNumberDensity);
      break;

    case format_type::kSixteenColumnsNuclear:
      primary_vars.push_back(Pq::Pressure);
      primary_vars.push_back(Pq::RestMassDensity);
      primary_vars.push_back(Pq::Enthalpy);
      primary_vars.push_back(Pq::BaryonNumberDensity);
      break;

    default:
      assert(false); // format not implemented
  }
}

//
// File reader (TODO: implement with templates)
//
template <class FR>
void XeosColdNuclear<FR>::ReadEosTables() {
  freader.Open();
  int file_len = freader.NumLines();
  int header_len;
  double fields[16];
  
  switch (format) {
    case format_type::kFourColumnsCGS:
      header_len = 1;
      num_fields = 4; // number of fields depends on format
      rho_table.ConvertTo(PhUnits::CGS);
      pres_table.ConvertTo(PhUnits::CGS);
      enth_table.ConvertTo(PhUnits::CGS);
      nbar_table.ConvertTo(PhUnits::CGS);

      // establish number of records in 1D table and resize the table    
      num_records = file_len - header_len;
      rho_table.Resize(num_records);
      pres_table.Resize(num_records);
      enth_table.Resize(num_records);
      nbar_table.Resize(num_records);
      break;

    case format_type::kSixteenColumnsNuclear:
      freader.Rewind();
      header_len = freader.SkipHashHeader();
      num_fields = 16;
      pres_table.ConvertTo(PhUnits::NUCLEAR);
      rho_table.ConvertTo(PhUnits::CGS);
      enth_table.ConvertTo(PhUnits::NUCLEAR);
      nbar_table.ConvertTo(PhUnits::NUCLEAR);

      // establish number of records in 1D table and resize the table    
      num_records = file_len - header_len;
      rho_table.Resize(num_records);
      pres_table.Resize(num_records);
      enth_table.Resize(num_records);
      nbar_table.Resize(num_records);
      break;

    default:
      assert (false);
  }
      
  freader.Rewind();
  freader.SkipHeader(header_len);
  for (int i=0; i<num_records; ++i) {
    if (freader.ReadFields(num_fields,fields) < num_fields) {
      num_records = i;
      break; // skip empty / incomplete lines in the end of file
    }

    switch (format) {
      case format_type::kFourColumnsCGS:
        rho_table(i)= fields[0];
        pres_table(i)= fields[1];
        enth_table(i)= fields[2];
        nbar_table(i)= fields[3];
        break;

      case format_type::kSixteenColumnsNuclear:
        pres_table(i)= fields[0];
        rho_table(i) = pow(10.,fields[1]);
        enth_table(i)= fields[2];
        nbar_table(i)= fields[3];
        break;

      default:
        assert (false);
    }
  }

  // convert to eos_units
  rho_table.ConvertTo(eos_units);
  pres_table.ConvertTo(eos_units);
  enth_table.ConvertTo(eos_units);
  eps_table.ConvertTo(eos_units);
  nbar_table.ConvertTo(eos_units);

} // ReadEosTable()


//
// Return a point of the one-dimensional grid
//
template <class FR>
double XeosColdNuclear<FR>::EosTableGridPoint (
    const Pq var, const int i) const {
  /* TODO */
  double retval;
  const NvDensity::size_type j = i;
  assert (i>=0 && i<num_records);
  switch(var) {
    case Pq::Density:
      retval = rho_table(j);
      break;

    case Pq::Pressure:
      retval = pres_table(j);
      break;

    case Pq::Enthalpy:
      retval = enth_table(j);
      break;

    case Pq::SpecificInternalEnergy:
      retval = eps_table(j);
      break;

    case Pq::BaryonNumberDensity:
      retval = nbar_table(j);
      break;

    default:
      assert (false);

  }
  return retval;
}

//
// General EoS function
//
template <class FR>
void XeosColdNuclear<FR>::operator() (const int num_out,
    eos_in in_array[],
    eos_out out_array[]) {
  /*      +
     TODO
   +      */
}

//
// General derivative(s), (dF/dX)_const
//
template <class FR>
void XeosColdNuclear<FR>::DfDx (const int num_out, const Pq dX,
    eos_in  in_array[],
    eos_out dF_array[]) {
  /*      +
     TODO
   +      */
}

template <class FR>
bool XeosColdNuclear<FR>::ConsistencyCheck() const {
  /*      +
     TODO
   +      */
  return false;
}

} // namespace xeos
#endif // XEOS_XEOSCOLDN_H_

#if 1
#include <iostream>

int main() {
  using namespace std;
  using namespace xeos;
  PhysicalQuantity::SetGlobalUnits(PhUnits::CGS);
  XeosColdNuclear<AsciiFileReader> 
      eos_SFHO("../data/rnsid/sfho_0.1MeV_beta.txt", 
      XeosColdNuclearFormat::kSixteenColumnsNuclear);
  XeosColdNuclear<AsciiFileReader> 
      eosA("../data/rnsid/eosA", 
      XeosColdNuclearFormat::kFourColumnsCGS);
  eos_SFHO.ReadEosTables();
  for (int i=0;i<eos_SFHO.EosTableGridSize(Pq::Density);++i) 
    cout << eos_SFHO.EosTableGridPoint(Pq::Density, i) << " " 
         << eos_SFHO.EosTableGridPoint(Pq::Pressure,i) << endl;
  cout << endl << endl;

  eosA.ReadEosTables();
  for (int i=0;i<eosA.EosTableGridSize(Pq::Density);++i) 
    cout << eosA.EosTableGridPoint(Pq::Density, i) << " " 
         << eosA.EosTableGridPoint(Pq::Pressure,i) << endl;

/* double fields[20];
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
 */
}
#endif
