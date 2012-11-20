#ifndef _uncertaintycomputer_h_
#define _uncertaintycomputer_h_

#if !defined(__CINT__) || defined(__MAKECINT__)

#include <iomanip>
#include <iostream>
#include <fstream>

#include "TRandom2.h"
#include "TMatrixD.h"
#include "TF1.h"
#include "TProfile.h"
#include "TObjArray.h"
#include "TMatrixD.h"
#include "TH1.h"
#include "TTimeStamp.h"
#include <exception>

#ifdef _STANDALONE
#include "Reader.h"
#else
#include "MiniTree/Selection/interface/Reader.h"
#endif
#include "MiniTree/Selection/interface/BtagSF.hh"

#endif

const double JEREtaMap[6] = {0., 0.5, 1.1, 1.7, 2.3, 5.0}; 
const double JERSF[5] = {1.052, 1.057, 1.096, 1.134, 1.288}; 
const double JERSFStatUn[5] = {0.012, 0.012, 0.017, 0.035, 0.127}; 
const double JERSFUp[5] = {0.062, 0.056, 0.063, 0.087, 0.155}; 
const double JERSFDown[5] = {0.061, 0.055, 0.062, 0.085, 0.153}; 

class UncertaintyComputer{

public :
  UncertaintyComputer()
  {
    btsf = new BtagSF(12345);
  }

  ~UncertaintyComputer(){
    delete btsf;
  }
  
  double metWithJES(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes=0);
  double metWithJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jer=0);
  double metWithJESJER(const vector<MyJet> & vJ, vector<int> *j, MyMET MET, int jes=0, int jer=0);
  double metWithUncl(const vector<MyJet> & vJ, vector<int> *j, const vector<MyMuon> &vMu, vector<int> *m, const vector<MyElectron> &vEle, vector<int> *el, MyMET MET, int unc=0);
  double getJERSF(double eta, int jer=0);
  double jetPtWithJESJER(MyJet jet, int jes=0, int jer=0); 
  bool getBtagWithSF(MyJet jet, bool isData, int scale, bool is2012);
  
private :
  BtagSF* btsf;
  enum BVariation{kNo = 0, kDown = 1, kUp = 2};

  ClassDef(UncertaintyComputer, 1)
};
#endif
