#ifndef _btagsf_hh_
#define _btagsf_hh_

#if !defined(__CINT__) || defined(__MAKECINT__)

#include "TRandom3.h"
#include <iostream>
#include <stdio.h>
#include <math.h>
#endif
using namespace std;
class BtagSF{
  
public:
  BtagSF( int seed=0 )
  { randm = new TRandom3(seed);};
  
  ~BtagSF() {delete randm;};
  
  Bool_t isbtagged(Float_t pt, Float_t eta, Float_t csv, Int_t jetflavor, Bool_t isdata, UInt_t btagsys, UInt_t mistagsys, Bool_t is2012);
  Double_t getSFb(Float_t pt, UInt_t btagsys);
  Double_t getSFl(Float_t pt, Float_t eta, UInt_t mistagsys, Bool_t is2012);
  Double_t getMistag(Float_t pt, Float_t eta);


  enum { kNo, kDown, kUp };                     // systematic variations 


private:
  
  TRandom3* randm;
  ClassDef(BtagSF, 1)
};
#endif
