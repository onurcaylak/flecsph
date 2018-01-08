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
  NvBaryonNumberDensity nbar_table {PhUnits::CGS};
  NvElectronFraction  ye_table {PhUnits::NUCLEAR};
  NvSpecificInternalEnergy eps_table {PhUnits::NUCLEAR};
  NvSoundSpeedSq sc2_table {PhUnits::CGS};
  NvChemPotentialDiff muhat_table {PhUnits::NUCLEAR};
  NvSpecificEntropy ent_table {PhUnits::NUCLEAR};

  NvNeutronMassFraction   Xn_table {PhUnits::NUCLEAR};
  NvProtonMassFraction    Xp_table {PhUnits::NUCLEAR};
  NvDeuteronMassFraction  Xd_table {PhUnits::NUCLEAR};
  NvTritonMassFraction    Xt_table {PhUnits::NUCLEAR};
  NvHelionMassFraction    Xh_table {PhUnits::NUCLEAR};
  NvAlphaMassFraction     Xa_table {PhUnits::NUCLEAR};
  NvHeavyNucMassFraction  Xheavy_table {PhUnits::NUCLEAR};
  NvAverageAtomicMass     Abar_table {PhUnits::NUCLEAR};
  NvAverageNuclearCharge  Zbar_table {PhUnits::NUCLEAR};
  
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
      primary_vars.push_back(Pq::Density);
      primary_vars.push_back(Pq::Pressure);
      primary_vars.push_back(Pq::Enthalpy);
      primary_vars.push_back(Pq::BaryonNumberDensity);

      num_fields = 4; // number of fields depends on format
      rho_table.ConvertTo(PhUnits::CGS);
      pres_table.ConvertTo(PhUnits::CGS);
      enth_table.ConvertTo(PhUnits::CGS);
      nbar_table.ConvertTo(PhUnits::CGS);
      break;

    case format_type::kSixteenColumnsNuclear:
      primary_vars.push_back(Pq::Pressure);
      primary_vars.push_back(Pq::Density);
      primary_vars.push_back(Pq::ElectronFraction);
      primary_vars.push_back(Pq::SpecificInternalEnergy);
      
      additional_vars.push_back(Pq::ChemPotentialDiff);
      additional_vars.push_back(Pq::SoundSpeedSq);
      additional_vars.push_back(Pq::SpecificEntropy);
      additional_vars.push_back(Pq::NeutronMassFraction);

      additional_vars.push_back(Pq::ProtonMassFraction);
      additional_vars.push_back(Pq::DeuteronMassFraction);
      additional_vars.push_back(Pq::TritonMassFraction);
      additional_vars.push_back(Pq::HelionMassFraction);

      additional_vars.push_back(Pq::AlphaMassFraction);
      additional_vars.push_back(Pq::HeavyNucMassFraction);
      additional_vars.push_back(Pq::AverageAtomicMass);
      additional_vars.push_back(Pq::AverageNuclearCharge);

      num_fields = 16;
      pres_table.ConvertTo(PhUnits::NUCLEAR);
      rho_table.ConvertTo(PhUnits::CGS);
      // ye_table.ConvertTo(PhUnits::NULEAR); // no need to convert Ye
      eps_table.ConvertTo(PhUnits::NUCLEAR);
      muhat_table.ConvertTo(PhUnits::NUCLEAR);
      sc2_table.ConvertTo(PhUnits::CGS);
      ent_table.ConvertTo(PhUnits::NUCLEAR);
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

      // establish number of records in 1D table and resize the table
      num_records = file_len - header_len;

      pres_table.Resize(num_records);
      rho_table.Resize(num_records);
      ye_table.Resize(num_records);
      eps_table.Resize(num_records);
      
      muhat_table.Resize(num_records);
      sc2_table.Resize(num_records);
      ent_table.Resize(num_records);
      Xn_table.Resize(num_records);
      
      Xp_table.Resize(num_records);
      Xd_table.Resize(num_records);
      Xt_table.Resize(num_records);
      Xh_table.Resize(num_records);
      
      Xa_table.Resize(num_records);
      Xheavy_table.Resize(num_records);
      Abar_table.Resize(num_records);
      Zbar_table.Resize(num_records);
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
        ye_table(i) =  fields[2];
        eps_table(i)=  fields[3];

        muhat_table(i)=fields[4];
        sc2_table(i)= fields[5];
        ent_table(i)= fields[6];

        Xn_table(i) = fields[7];
        Xp_table(i) = fields[8];
        Xd_table(i) = fields[9];
        Xt_table(i) = fields[10];
        Xh_table(i) = fields[11];
        Xa_table(i) = fields[12];
        Xheavy_table(i) = fields[13];
        Abar_table(i) = fields[14];
        Zbar_table(i) = fields[15];
        break;

      default:
        assert (false);
    }
  }

  // convert to eos_units
  rho_table.ConvertTo(eos_units);
  pres_table.ConvertTo(eos_units);
  enth_table.ConvertTo(eos_units);
  nbar_table.ConvertTo(eos_units);
  ye_table.ConvertTo(eos_units);
  eps_table.ConvertTo(eos_units);
  muhat_table.ConvertTo(eos_units);
  sc2_table.ConvertTo(eos_units);
  ent_table.ConvertTo(eos_units);

  Xn_table.ConvertTo(eos_units);
  Xp_table.ConvertTo(eos_units);
  Xd_table.ConvertTo(eos_units);
  Xt_table.ConvertTo(eos_units);
  Xh_table.ConvertTo(eos_units);
  Xa_table.ConvertTo(eos_units);
  Xheavy_table.ConvertTo(eos_units);
  Abar_table.ConvertTo(eos_units);
  Zbar_table.ConvertTo(eos_units);

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
    case Pq::BaryonNumberDensity:
      retval = nbar_table(j);
      break;

    case Pq::ChemPotentialDiff:
      retval = muhat_table(j);
      break;

    case Pq::Density:
      retval = rho_table(j);
      break;

    case Pq::ElectronFraction:
      retval = ye_table(j);
      break;

    case Pq::Enthalpy:
      retval = enth_table(j);
      break;

    case Pq::Pressure:
      retval = pres_table(j);
      break;

    case Pq::SoundSpeedSq:
      retval = sc2_table(j);
      break;

    case Pq::SpecificEntropy:
      retval = ent_table(j);
      break;

    case Pq::SpecificInternalEnergy:
      retval = eps_table(j);
      break;

    case Pq::NeutronMassFraction:
      retval = Xn_table(j);
      break;

    case Pq::ProtonMassFraction:
      retval = Xp_table(j);
      break;

    case Pq::DeuteronMassFraction:
      retval = Xd_table(j);
      break;

    case Pq::TritonMassFraction:
      retval = Xt_table(j);
      break;

    case Pq::HelionMassFraction:
      retval = Xh_table(j);
      break;

    case Pq::AlphaMassFraction:
      retval = Xa_table(j);
      break;

    case Pq::HeavyNucMassFraction:
      retval = Xheavy_table(j);
      break;

    case Pq::AverageAtomicMass:
      retval = Abar_table(j);
      break;

    case Pq::AverageNuclearCharge:
      retval = Zbar_table(j);
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

  PhysicalQuantity::SetGlobalUnits(PhUnits::NUCLEAR);

  XeosColdNuclear<AsciiFileReader> 
      eSFHo("../data/rnsid/sfho_0.1MeV_beta.txt", 
      XeosColdNuclearFormat::kSixteenColumnsNuclear);
  XeosColdNuclear<AsciiFileReader> 
      eosA("../data/rnsid/eosA", 
      XeosColdNuclearFormat::kFourColumnsCGS);
  
  eSFHo.ReadEosTables();
  cout << scientific << setprecision(12) << setw(18);
  for (int i=0;i<eSFHo.EosTableGridSize(Pq::Density);++i) 
    cout << eSFHo.EosTableGridPoint(Pq::Density, i) << " " 
         << eSFHo.EosTableGridPoint(Pq::Pressure,i) << " "
         << eSFHo.EosTableGridPoint(Pq::ElectronFraction,i) << " "
         << eSFHo.EosTableGridPoint(Pq::ChemPotentialDiff,i) << " " 
         << eSFHo.EosTableGridPoint(Pq::SoundSpeedSq,i) << " "
         << eSFHo.EosTableGridPoint(Pq::SpecificEntropy,i) << " "
         << eSFHo.EosTableGridPoint(Pq::NeutronMassFraction,i) << " "
         << eSFHo.EosTableGridPoint(Pq::ProtonMassFraction,i) << " "
         << eSFHo.EosTableGridPoint(Pq::DeuteronMassFraction,i) << " "
         << eSFHo.EosTableGridPoint(Pq::TritonMassFraction,i) << " "
         << eSFHo.EosTableGridPoint(Pq::HelionMassFraction,i) << " "
         << eSFHo.EosTableGridPoint(Pq::AlphaMassFraction,i) << " "
         << eSFHo.EosTableGridPoint(Pq::HeavyNucMassFraction,i) << " "
         << eSFHo.EosTableGridPoint(Pq::AverageAtomicMass,i) << " "
         << eSFHo.EosTableGridPoint(Pq::AverageNuclearCharge,i) << " "
         << endl;
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
