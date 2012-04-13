#include "MiniTree/Selection/interface/ObjectSelector.hh"
#include <iostream>
#include <iomanip>

ClassImp(ObjectSelector)

using namespace std;
void ObjectSelector::preSelectElectrons(vector<int> * e_i, const vector<MyElectron> & vE , MyVertex & vertex, bool isPFlow){
  
  for(unsigned int i=0;i<vE.size();i++){

    const MyElectron * e = &vE[i];
    double scEta    = e->electronSCEta;
    double eEta     = TMath::Abs(e->p4.eta());
    double eEt      = TMath::Abs(e->p4.Et());
    double ePt      = TMath::Abs(e->p4.pt());
    double eRelIso  = e->RelIso;
    double d0       = fabs(e->D0);
    std::map<std::string, float>idWPs = e->eidWPs;
    float eid = (isPFlow) ? idWPs["eidTightMC"] : idWPs["simpleEleId70cIso"];
    bool passId = (int(eid) & 0x1);
    bool isNotFromConversion = ((int(eid) >> 2) & 0x1);
    bool isEcalDriven = (e->isEcalDriven > 0 || e->isPFlow > 0);
    
    double zvertex   = vertex.XYZ.z();
    double zelectron = e->vertex.z();
    double dz = fabs(zvertex - zelectron);
    
    if( defaultSelection_  && ( fabs(scEta)> 1.4442 && fabs(scEta)<1.5660)  ){ continue; }
    if( defaultSelection_  && eEt < E_ET_MIN_ ){ continue; }
    if( !defaultSelection_ && ePt < E_PT_MIN_ ){ continue; }

    if(dz > ZMAX_)continue;

    if(passId && isNotFromConversion && isEcalDriven  && d0 < E_D0_MAX_ && eEta < E_ETA_MAX_ && eRelIso < E_RELISO_MAX_){ e_i->push_back(i);}

  }
}

void ObjectSelector::preSelectMuons(vector<int> * m_i, const vector<MyMuon> & vM , MyVertex & vertex, bool isPFlow){
  
  for( int i=0;i< (int) vM.size();i++){
    
    const MyMuon * m = &vM[i];

    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    double mD0      = fabs(m->D0);
    double mRelIso  = (isPFlow) ? m->RelIso : m->UserPFRelIso;

    bool passId = (m->GlobalMuonPromptTight > 0);

    double zvertex   = vertex.XYZ.z();
    double zmuon     = m->vertex.z();
    double dz =  fabs(zvertex-zmuon);
    
    if(dz > ZMAX_)continue;
    
    if(passId && mD0 < M_D0_MAX_&&  mPt > M_PT_MIN_ && mEta < M_ETA_MAX_ && mRelIso < M_RELISO_MAX_){ m_i->push_back(i);}
    
  }
  
}

void ObjectSelector::preSelectTaus( vector<int> * t_i, const vector<MyTau> & vT, int requiredProngs ,string myKey, MyVertex & vertex){
  
  
  for(unsigned int i=0;i<vT.size();i++){
    const MyTau *tau = &vT[i];
    double tauEta  = fabs(tau->p4.eta());
    double tauPt   = tau->p4.Pt();
    //double d0      = fabs( tau->d0_max );
    int charge     = fabs( tau->charge  );

    double zvertex  = vertex.XYZ.z();
    double ztau     = tau->vertex.z();
    double dz       = fabs(zvertex - ztau);
    int prongs      = int(tau->nSignalChargeHadron);

    if( tauPt < TAU_PT_MIN_ || tauEta > TAU_ETA_MAX_ ) continue;
    //if ( dz > ZTAUMAX_ )continue;
    //if(requiredProngs != 0 && prongs != requiredProngs)continue;
    if(tau->decayModeFinding < 0.5)continue;
    
    int tauIsolation = 0;
    if(      myKey =="Loose"    ){tauIsolation = int(tau->LooseIsolation);   }
    else if( myKey =="Medium"   ){tauIsolation = int(tau->MediumIsolation);  }
    else if( myKey =="Tight"    ){tauIsolation = int(tau->TightIsolation);   }
    else if( myKey =="LooseDB"  ){tauIsolation = int(tau->LooseIsolationDeltaBetaCorr);   }
    else if( myKey =="MediumDB"   ){tauIsolation = int(tau->MediumIsolationDeltaBetaCorr);  }
    else if( myKey =="TightDB"    ){tauIsolation = int(tau->TightIsolationDeltaBetaCorr);   }
    else if( myKey =="LooseCombDB"){tauIsolation = int(tau->LooseCombinedIsolationDeltaBetaCorr);   }
    else if( myKey =="MediumCombDB"){tauIsolation = int(tau->MediumCombinedIsolationDeltaBetaCorr);  }
    else if( myKey =="TightCombDB"){tauIsolation = int(tau->TightCombinedIsolationDeltaBetaCorr);   }

    if(tauIsolation < 0.5) continue;
    
    //if(tau->againstElectronMedium < 0.5 || tau->againstMuonLoose < 0.5)
    //  continue;
    
    //if(charge != 1)continue;

    t_i->push_back(i);
  }
}
    
void ObjectSelector::preSelectJets( string jetAlgo, vector<int> * j_i, const vector<MyJet> & vJ){
  
  for(unsigned int i=0;i<vJ.size();i++){
    
    const MyJet *jet = &vJ[i];
    double jetEta     = TMath::Abs(jet->p4.eta());
    double jetPt      = TMath::Abs(jet->p4.pt());

    if(jetPt > JET_PT_MIN_ && jetEta < JET_ETA_MAX_ )
      { j_i->push_back(i);}
  }
}
bool ObjectSelector::looseMuonVeto( int selectedMuon, const vector<MyMuon> & vM, bool isPFlow){

  bool looseVeto(false);
  
  for(int i=0;i< (int)vM.size();i++){
    
    if( i==selectedMuon ){continue;}
    
    const MyMuon * m = &vM[i];
    
    double mEta     = TMath::Abs(m->p4.eta());
    double mPt      = TMath::Abs(m->p4.pt());
    double mRelIso  = (isPFlow) ? m->RelIso : m->UserPFRelIso;
    
    //see if this muon is glogal
    static const unsigned int GlobalMuon     =  1<<1;
    bool isGlobal = (m->type & GlobalMuon);
    
    if(! isGlobal ) continue;
    
    if( mEta<LOOSE_M_ETA_MAX_  && mPt> LOOSE_M_PT_MIN_ && mRelIso < LOOSE_M_RELISO_MAX_ ){ looseVeto = true; }
    
  }
  
  return looseVeto;
    
}

bool ObjectSelector::looseElectronVeto(int selectedElectron, const vector<MyElectron> & vE, bool isPFlow){
  bool looseVeto(false);
  for( int i=0;i< (int) vE.size();i++ ){
    if(i==selectedElectron){ continue; }
    
    const MyElectron * e = &vE[i];
    
    double eEta      = TMath::Abs(e->p4.eta());
    double eEt       = TMath::Abs(e->p4.Et());
    double eRelIso   = e->RelIso;
    std::map<std::string, float>idWPs = e->eidWPs;
    float eid = (isPFlow) ? idWPs["eidLooseMC"] : idWPs["simpleEleId90cIso"];
    bool id = (int(eid) & 0x1);
    double minDR2mu  = 0.4;
        
    if(id && eEt > LOOSE_E_ET_MIN_ && eEta < LOOSE_E_ETA_MAX_ && eRelIso < LOOSE_E_RELISO_MAX_ && minDR2mu >0.1){ looseVeto = true;}
    
     
    
     
  }
  
  return looseVeto;
  
}

void ObjectSelector::ElectronCleaning( const vector<MyElectron> & vE, const vector<MyMuon> & vM, vector<int> * e_old, vector<int> * e_new, vector<int> * mu, double DR ){
  
  for(size_t i = 0; i < e_old->size(); i++){
    int iele = (*e_old)[i];
    double delR2Mu = 5.0;
    for(size_t j = 0; j < mu->size(); j++){
      int imu = (*mu)[j];
      double delR = DeltaR(vE[iele].p4, vM[imu].p4);
      if(delR < delR2Mu)delR2Mu = delR;
    }

    if(delR2Mu > DR) e_new->push_back(iele);
  }
}

void ObjectSelector::TauCleaning(const vector<MyTau> & vT, const vector<MyMuon> & vM, const vector<MyElectron> & vE, vector<int> * t_old, vector<int> * t_new, vector<int> * mu, vector<int> * el, double DR ){
  
  for(size_t i = 0; i < t_old->size(); i++){
    int itau = (*t_old)[i];

    double delR2Mu = 5.0, delR2Ele = 5.0;
    
    for(size_t j = 0; j < mu->size(); j++){
      int imu = (*mu)[j];
      double delR = DeltaR(vT[itau].p4, vM[imu].p4);
      if(delR < delR2Mu)delR2Mu = delR;
    }

    for(size_t k = 0; k < el->size(); k++){
      int iele = (*el)[k];
      double delR = DeltaR(vT[itau].p4, vE[iele].p4);
      if(delR < delR2Ele)delR2Ele = delR;
    }

    if(delR2Mu > DR && delR2Ele > DR)t_new->push_back(itau);
  }
}


void ObjectSelector::JetCleaning(const vector<MyJet> & vJ, const vector<MyMuon> & vM, const vector<MyElectron> & vE, const vector<MyTau> & vT, vector<int> * j_old, vector<int> * j_new, vector<int> * mu, vector<int> * el, vector<int> * tau, double DR, bool cleanTaus){

  for(size_t i = 0; i < j_old->size(); i++){
    int ijet = (*j_old)[i];

    double delR2Mu = 5.0, delR2Ele = 5.0, delR2Tau = 5.0;
    
    for(size_t j = 0; j < mu->size(); j++){
      int imu = (*mu)[j];
      double delR = DeltaR(vJ[ijet].p4, vM[imu].p4);
      if(delR < delR2Mu)delR2Mu = delR;
    }

    for(size_t k = 0; k < el->size(); k++){
      int iele = (*el)[k];
      double delR = DeltaR(vJ[ijet].p4, vE[iele].p4);
      if(delR < delR2Ele)delR2Ele = delR;
    }
    if(cleanTaus){
      for(size_t l = 0; l < tau->size(); l++){
      int itau = (*tau)[l];
      double delR = DeltaR(vJ[ijet].p4, vT[itau].p4);
      if(delR < delR2Tau)delR2Tau = delR;
      }
    }

    if(delR2Mu > DR && delR2Ele > DR && delR2Tau > DR)j_new->push_back(ijet);
  }
}

double ObjectSelector::DeltaR(MyLorentzVector aV, MyLorentzVector bV){
  
  double deta = TMath::Abs(aV.eta() - bV.eta());
  double dphi = TMath::Abs(aV.phi() - bV.phi());
  if(dphi > M_PI) dphi = 2*M_PI - dphi;

  double delR = sqrt(deta*deta + dphi*dphi);

  return delR;
}

