// Author: A. Rustamov
// Written: 
// Revised: 

#ifndef _PRESONANCEDALITZS_H_
#define _PRESONANCEDALITZS_H_

#include "PDalitzDecay.h"

class PResonanceDalitz : public PDalitzDecay  
{
  
public:

  using PDalitzDecay::GetWeight;
  PResonanceDalitz(const Char_t *id, const Char_t *de, Int_t key);
  PDistribution *Clone(const char *delme=NULL) const;
  virtual double dGdM(const int& id, const double& m, const double& ecm);
  
  ////////////////////
  double GetMatrixT(int&, int&, int&, const double&, const double&);
  double GetMatrixL(int&, int&, int&, const double&, const double&);
  double getLambda(double, double, double);

  
  double ecm2;
  double ecm4;
  double ecm3;
  double ecmPmn2;
  double ecmMmn2;
  double M2;
  double M4;  
  
  double mn;
  double mn2;
  double mn3;
  double mn4;

  Double_t  g_Em;
  Int_t spin,par,charge;
  
  void setGm(Double_t my_g_Em) {
      g_Em = my_g_Em;
  }


  ////////////////////    
  ClassDef(PResonanceDalitz, 0)  //Dalitz decay of N* resonances
};

#endif
