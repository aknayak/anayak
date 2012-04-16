#ifndef _objectselector_h_
#define _objectselector_h_

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

#endif

class ObjectSelector{

public :
  ObjectSelector()
  {
    setDefaultSelection();
  }
  ~ObjectSelector(){}

  void setDefaultSelection(){
    defaultSelection_=true;

    //Vertex //////
    ZMAX_    = 1;
    ZTAUMAX_ = 0.2;
    ///////////////

    // Jets ///////////////////////////////////////////
    JET_PT_MIN_        = 20;
    JET_ETA_MAX_       = 2.3; 
    JET_EMF_MIN_       = 0.01;
    JET_LEPTON_DRMIN_  = 0.5;
    //////////////////////////////////////////////////

    // electron //////////////////////////////////////
    E_RELISO_MAX_       = 0.1;
    E_ETA_MAX_          = 2.5;
    E_ET_MIN_           = 20;
    E_D0_MAX_           = 0.04;

    LOOSE_E_RELISO_MAX_ = 0.2;
    LOOSE_E_ETA_MAX_    = 2.5;
    LOOSE_E_ET_MIN_     = 10;
    RHO_AEFF_E_         = 0.24;
    //////////////////////////////////////////////////


    // muon //////////////////////////////////////////
    M_RELISO_MAX_  = 0.1; 
    M_PT_MIN_      = 30;  
    M_ETA_MAX_     = 2.1; 
    M_D0_MAX_      = 0.02;
    
    LOOSE_M_RELISO_MAX_ = 0.2;
    LOOSE_M_ETA_MAX_    = 2.5;
    LOOSE_M_PT_MIN_     = 10;
    RHO_AEFF_M_         = 0.112;
    ///////////////////////////////////////////////////


    // taus ///////////////////////////////////////////////////////
    TAU_ETA_MAX_          = 2.3;
    TAU_PT_MIN_           = 15.0;
    TAU_D0_MAX_           = 400;// was 200
    
  }

  // preselection of objects
  void preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex, bool isPFlow=false);
  void preSelectMuons(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isPFlow=false);
  void preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ);
  void preSelectTaus( vector<int> * t_i, const vector<MyTau> & vT, int requiredProngs ,string type, MyVertex & vertex);

  //Loose Lepton veto
  bool looseElectronVeto(int selectedElectron, const vector<MyElectron> & vE, bool isPFlow=false);
  bool looseMuonVeto( int selectedMuon, const vector<MyMuon> & vM, bool isPFlow=false);


  // object cleaning
  void ElectronCleaning( const vector<MyElectron> & vE, const vector<MyMuon> & vM, vector<int> * e_old, vector<int> * e_new, vector<int> * mu, double DR );

  void TauCleaning(const vector<MyTau> & vT, const vector<MyMuon> & vM, const vector<MyElectron> & vE, vector<int> * t_old, vector<int> * t_new, vector<int> * mu, vector<int> * el, double DR );
  void JetCleaning(const vector<MyJet> & vJ, const vector<MyMuon> & vM, const vector<MyElectron> & vE, const vector<MyTau> & vT, vector<int> * j_old, vector<int> * j_new, vector<int> * mu, vector<int> * el, vector<int> * tau, double DR, bool cleanTaus = false);
  double DeltaR(MyLorentzVector aV, MyLorentzVector bV); 
    
  
private :
  // Vertex
  double ZMAX_, ZTAUMAX_ ;

  // jet
  double JET_PT_MIN_,JET_ETA_MAX_,JET_LTK_PT_MIN_, JET_BTAGGING_, JET_EMF_MIN_, JET_LTK_,JET_LEPTON_DRMIN_  ;
  
  // electron
  double E_RELISO_MAX_, E_ETA_MAX_, E_D0_MAX_, E_PT_MIN_, E_ET_MIN_,RHO_AEFF_E_;
  double LOOSE_E_RELISO_MAX_, LOOSE_E_ETA_MAX_, LOOSE_E_D0_MAX_, LOOSE_E_ET_MIN_;
  
  
  // muon
  double M_RELISO_MAX_, M_ETA_MAX_, M_D0_MAX_, M_PT_MIN_,RHO_AEFF_M_;
  double LOOSE_M_RELISO_MAX_, LOOSE_M_ETA_MAX_, LOOSE_M_D0_MAX_, LOOSE_M_PT_MIN_;
  
  // taus
  double TAU_ETA_MAX_,TAU_D0_MAX_, TAU_PT_MIN_;
  int HPS_ISO_;
  
  bool defaultSelection_;

  ClassDef(ObjectSelector, 1)
};
#endif
