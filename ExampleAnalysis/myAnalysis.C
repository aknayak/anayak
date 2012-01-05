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
#include "TFile.h"
#include "TTree.h"

#include "selectEvent.C"
//#include "MiniTree/Selection/interface/Reader.h"
//#include "GetEvent.h"
//#include "GetEvent.C"

//#include "MiniTree/Selection/interface/MyEvent.h"

using namespace std;

int main()
{
  selectEvent t;
  cout<<" from mc file directly "<<endl;
  t.checkEvent("mc_tau.root");
  cout<<" from list of files "<<endl;
  t.checkEvent1("my_file_list");
  return 0; 
}

