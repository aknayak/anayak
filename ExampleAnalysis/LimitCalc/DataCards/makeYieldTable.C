#include <iostream>
#include <fstream>
#include <iomanip>
#include <algorithm>

void makeCutFlowTable()
{

  TFile *wh = new TFile("wjet_selection.root");
  TH1F* h_wh = ((TH1F*)wh->Get("base/cutflow") )->Clone("h_wh");

  TFile *ttbar = new TFile("wjet_selection.root");
  TH1F* h_ttbar = ((TH1F*)ttbar->Get("base/cutflow") )->Clone("h_ttbar");
  
  TFile *wjet = new TFile("wjet_selection.root");
  TH1F* h_wjet = ((TH1F*)wjet->Get("base/cutflow") )->Clone("h_wjet");
  
  h_ttbar->Sumw2();
  h_wjet->Sumw2();
  
  TH1F* h_TotalBkg = h_wjet->Clone("h_TotalBkg");
  h_TotalBkg->Reset();
  h_TotalBkg->Add(h_ttbar);
  h_TotalBkg->Add(h_wjet);


  ofstream outFile;
  outFile.open("cutflow.tex");
  
  outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl; 
  outFile<<"\\pagestyle{empty}"<<endl; 
  outFile<<"\\usepackage{epsfig}"<<endl; 
  outFile<<"\\usepackage{amsmath}"<<endl; 
  outFile<<"\\usepackage{array}"<<endl; 
  outFile<<"\\usepackage{multirow}"<<endl; 
  outFile<<"\\usepackage[cm]{fullpage}"<<endl; 
  outFile<<"\\textheight = 8.in"<<endl; 
  outFile<<"\\textwidth 7.0in"<<endl; 
  outFile<<"\\begin{document}"<<endl; 
  outFile<<"\\begin{center}"<<endl; 
  outFile<<"\\begin{LARGE}"<<endl; 
  outFile<<"\\begin{tabular}{ | c| c| c| c| c| c|}"<<endl; 
  outFile<<"\\multicolumn{6}{c}{ } \\\\"<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<<"\\multicolumn{1}{| c|}{ } & \\multicolumn{1}{ c|}{ $N_{muon}=1$ } & \\multicolumn{1}{ c|}{ $N_{jets}\\ge 4$ } & \\multicolumn{1}{ c|}{ $\\not\\!\\!E_T \\ge 30GeV$ } & \\multicolumn{1}{ c |}{ $\\ge$ 1btag } & \\multicolumn{1}{ c |}{ KinFit} \\\\ "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;
  outFile<<"HW, $M_{H}=120~GeV/c^{2}$"<<" & "<<h_wh->GetBinContent(2)<<" & "<<h_wh->GetBinContent(3)<<" & "<<h_wh->GetBinContent(4)<<" & "<<h_wh->GetBinContent(5)<<" & "<<h_wh->GetBinContent(6)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"SM $t\\bar{t}$"<<" & "<<h_ttbar->GetBinContent(2)<<" & "<<h_ttbar->GetBinContent(3)<<" & "<<h_ttbar->GetBinContent(4)<<" & "<<h_ttbar->GetBinContent(5)<<" & "<<h_ttbar->GetBinContent(6)<<" \\\\ "<<endl;
  outFile<<"W+Jets"<<" & "<<h_wjet->GetBinContent(2)<<" & "<<h_wjet->GetBinContent(3)<<" & "<<h_wjet->GetBinContent(4)<<" & "<<h_wjet->GetBinContent(5)<<" & "<<h_wjet->GetBinContent(6)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<<"Total Bkg"<<" & "<<h_TotalBkg->GetBinContent(2)<<" $\\pm$ "<<h_TotalBkg->GetBinError(2)<<" & "<<h_TotalBkg->GetBinContent(3)<<" $\\pm$ "<<h_TotalBkg->GetBinError(3)<<" & "<<h_TotalBkg->GetBinContent(4)<<" $\\pm$ "<<h_TotalBkg->GetBinError(4)<<" & "<<h_TotalBkg->GetBinContent(5)<<" $\\pm$ "<<h_TotalBkg->GetBinError(5)<<" & "<<h_TotalBkg->GetBinContent(6)<<" $\\pm$ "<<h_TotalBkg->GetBinError(6)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<"\\hline "<<endl; 
  outFile<<"Data "<<" & "<<h_TotalBkg->GetBinContent(2)<<" & "<<h_TotalBkg->GetBinContent(3)<<" & "<<h_TotalBkg->GetBinContent(4)<<" & "<<h_TotalBkg->GetBinContent(5)<<" & "<<h_TotalBkg->GetBinContent(6)<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\hline "<<endl;  
  outFile<<"\\end{tabular}"<<endl;
  outFile<<"\\end{LARGE}"<<endl; 
  outFile<<"\\end{center}"<<endl; 
  outFile<<"\\end{document}"<<endl;

  outFile.close();
}


void makeSummaryTableWithSys() 
{ 
  
  TString inputFilePath("/afs/cern.ch/user/g/gkole/public/chargedHiggs/8TeV/2012ABC/uncertainty/kinfit4_v10"); 
  TFile *wh = new TFile(inputFilePath+"/signal_M_120_selection.root"); 
  TH1F* h_wh_base = ((TH1F*)wh->Get("base/cutflow") )->Clone("h_wh_base"); 
  TH1F* h_wh_JESPlus = ((TH1F*)wh->Get("JESPlus/cutflow") )->Clone("h_wh_JESPlus");
  TH1F* h_wh_JESMinus = ((TH1F*)wh->Get("JESMinus/cutflow") )->Clone("h_wh_JESMinus");
  TH1F* h_wh_JERPlus = ((TH1F*)wh->Get("JERPlus/cutflow") )->Clone("h_wh_JERPlus");
  TH1F* h_wh_JERMinus = ((TH1F*)wh->Get("JERMinus/cutflow") )->Clone("h_wh_JERMinus");
  TH1F* h_wh_METUCPlus = ((TH1F*)wh->Get("METUCPlus/cutflow") )->Clone("h_wh_METUCPlus");
  TH1F* h_wh_METUCMinus = ((TH1F*)wh->Get("METUCMinus/cutflow") )->Clone("h_wh_METUCMinus");
  TH1F* h_wh_bTagPlus = ((TH1F*)wh->Get("bTagPlus/cutflow") )->Clone("h_wh_bTagPlus");
  TH1F* h_wh_bTagMinus = ((TH1F*)wh->Get("bTagMinus/cutflow") )->Clone("h_wh_bTagMinus");
 
  TFile *ttbar = new TFile(inputFilePath+"/ttbar_selection.root"); 
  TH1F* h_ttbar_base = ((TH1F*)ttbar->Get("base/cutflow") )->Clone("h_ttbar_base");  
  TH1F* h_ttbar_JESPlus = ((TH1F*)ttbar->Get("JESPlus/cutflow") )->Clone("h_ttbar_JESPlus"); 
  TH1F* h_ttbar_JESMinus = ((TH1F*)ttbar->Get("JESMinus/cutflow") )->Clone("h_ttbar_JESMinus"); 
  TH1F* h_ttbar_JERPlus = ((TH1F*)ttbar->Get("JERPlus/cutflow") )->Clone("h_ttbar_JERPlus"); 
  TH1F* h_ttbar_JERMinus = ((TH1F*)ttbar->Get("JERMinus/cutflow") )->Clone("h_ttbar_JERMinus"); 
  TH1F* h_ttbar_METUCPlus = ((TH1F*)ttbar->Get("METUCPlus/cutflow") )->Clone("h_ttbar_METUCPlus"); 
  TH1F* h_ttbar_METUCMinus = ((TH1F*)ttbar->Get("METUCMinus/cutflow") )->Clone("h_ttbar_METUCMinus"); 
  TH1F* h_ttbar_bTagPlus = ((TH1F*)ttbar->Get("bTagPlus/cutflow") )->Clone("h_ttbar_bTagPlus"); 
  TH1F* h_ttbar_bTagMinus = ((TH1F*)ttbar->Get("bTagMinus/cutflow") )->Clone("h_ttbar_bTagMinus"); 

   
  TFile *wjet = new TFile(inputFilePath+"/wjet_selection.root"); 
  TH1F* h_wjet_base = ((TH1F*)wjet->Get("base/cutflow") )->Clone("h_wjet_base");  
  TH1F* h_wjet_JESPlus = ((TH1F*)wjet->Get("JESPlus/cutflow") )->Clone("h_wjet_JESPlus"); 
  TH1F* h_wjet_JESMinus = ((TH1F*)wjet->Get("JESMinus/cutflow") )->Clone("h_wjet_JESMinus"); 
  TH1F* h_wjet_JERPlus = ((TH1F*)wjet->Get("JERPlus/cutflow") )->Clone("h_wjet_JERPlus"); 
  TH1F* h_wjet_JERMinus = ((TH1F*)wjet->Get("JERMinus/cutflow") )->Clone("h_wjet_JERMinus"); 
  TH1F* h_wjet_METUCPlus = ((TH1F*)wjet->Get("METUCPlus/cutflow") )->Clone("h_wjet_METUCPlus"); 
  TH1F* h_wjet_METUCMinus = ((TH1F*)wjet->Get("METUCMinus/cutflow") )->Clone("h_wjet_METUCMinus"); 
  TH1F* h_wjet_bTagPlus = ((TH1F*)wjet->Get("bTagPlus/cutflow") )->Clone("h_wjet_bTagPlus"); 
  TH1F* h_wjet_bTagMinus = ((TH1F*)wjet->Get("bTagMinus/cutflow") )->Clone("h_wjet_bTagMinus"); 

  TFile *zjet = new TFile(inputFilePath+"/zjet_selection.root");  
  TH1F* h_zjet_base = ((TH1F*)zjet->Get("base/cutflow") )->Clone("h_zjet_base");   
  TH1F* h_zjet_JESPlus = ((TH1F*)zjet->Get("JESPlus/cutflow") )->Clone("h_zjet_JESPlus");  
  TH1F* h_zjet_JESMinus = ((TH1F*)zjet->Get("JESMinus/cutflow") )->Clone("h_zjet_JESMinus");  
  TH1F* h_zjet_JERPlus = ((TH1F*)zjet->Get("JERPlus/cutflow") )->Clone("h_zjet_JERPlus");  
  TH1F* h_zjet_JERMinus = ((TH1F*)zjet->Get("JERMinus/cutflow") )->Clone("h_zjet_JERMinus");  
  TH1F* h_zjet_METUCPlus = ((TH1F*)zjet->Get("METUCPlus/cutflow") )->Clone("h_zjet_METUCPlus");  
  TH1F* h_zjet_METUCMinus = ((TH1F*)zjet->Get("METUCMinus/cutflow") )->Clone("h_zjet_METUCMinus");  
  TH1F* h_zjet_bTagPlus = ((TH1F*)zjet->Get("bTagPlus/cutflow") )->Clone("h_zjet_bTagPlus");  
  TH1F* h_zjet_bTagMinus = ((TH1F*)zjet->Get("bTagMinus/cutflow") )->Clone("h_zjet_bTagMinus");

  TFile *stop = new TFile(inputFilePath+"/singletop_selection.root");  
  TH1F* h_stop_base = ((TH1F*)stop->Get("base/cutflow") )->Clone("h_stop_base");   
  TH1F* h_stop_JESPlus = ((TH1F*)stop->Get("JESPlus/cutflow") )->Clone("h_stop_JESPlus");  
  TH1F* h_stop_JESMinus = ((TH1F*)stop->Get("JESMinus/cutflow") )->Clone("h_stop_JESMinus");  
  TH1F* h_stop_JERPlus = ((TH1F*)stop->Get("JERPlus/cutflow") )->Clone("h_stop_JERPlus");  
  TH1F* h_stop_JERMinus = ((TH1F*)stop->Get("JERMinus/cutflow") )->Clone("h_stop_JERMinus");  
  TH1F* h_stop_METUCPlus = ((TH1F*)stop->Get("METUCPlus/cutflow") )->Clone("h_stop_METUCPlus");  
  TH1F* h_stop_METUCMinus = ((TH1F*)stop->Get("METUCMinus/cutflow") )->Clone("h_stop_METUCMinus");  
  TH1F* h_stop_bTagPlus = ((TH1F*)stop->Get("bTagPlus/cutflow") )->Clone("h_stop_bTagPlus");  
  TH1F* h_stop_bTagMinus = ((TH1F*)stop->Get("bTagMinus/cutflow") )->Clone("h_stop_bTagMinus");
   
  TFile *vv = new TFile(inputFilePath+"/diboson_selection.root");  
  TH1F* h_vv_base = ((TH1F*)vv->Get("base/cutflow") )->Clone("h_vv_base");   
  TH1F* h_vv_JESPlus = ((TH1F*)vv->Get("JESPlus/cutflow") )->Clone("h_vv_JESPlus");  
  TH1F* h_vv_JESMinus = ((TH1F*)vv->Get("JESMinus/cutflow") )->Clone("h_vv_JESMinus");  
  TH1F* h_vv_JERPlus = ((TH1F*)vv->Get("JERPlus/cutflow") )->Clone("h_vv_JERPlus");  
  TH1F* h_vv_JERMinus = ((TH1F*)vv->Get("JERMinus/cutflow") )->Clone("h_vv_JERMinus");  
  TH1F* h_vv_METUCPlus = ((TH1F*)vv->Get("METUCPlus/cutflow") )->Clone("h_vv_METUCPlus");  
  TH1F* h_vv_METUCMinus = ((TH1F*)vv->Get("METUCMinus/cutflow") )->Clone("h_vv_METUCMinus");  
  TH1F* h_vv_bTagPlus = ((TH1F*)vv->Get("bTagPlus/cutflow") )->Clone("h_vv_bTagPlus");  
  TH1F* h_vv_bTagMinus = ((TH1F*)vv->Get("bTagMinus/cutflow") )->Clone("h_vv_bTagMinus");
  
  TFile *qcd = new TFile(inputFilePath+"/qcd_selection.root");  
  TH1F* h_qcd_base = ((TH1F*)qcd->Get("base/cutflow") )->Clone("h_qcd_base");   
  TH1F* h_qcd_JESPlus = ((TH1F*)qcd->Get("JESPlus/cutflow") )->Clone("h_qcd_JESPlus");  
  TH1F* h_qcd_JESMinus = ((TH1F*)qcd->Get("JESMinus/cutflow") )->Clone("h_qcd_JESMinus");  
  TH1F* h_qcd_JERPlus = ((TH1F*)qcd->Get("JERPlus/cutflow") )->Clone("h_qcd_JERPlus");  
  TH1F* h_qcd_JERMinus = ((TH1F*)qcd->Get("JERMinus/cutflow") )->Clone("h_qcd_JERMinus");  
  TH1F* h_qcd_METUCPlus = ((TH1F*)qcd->Get("METUCPlus/cutflow") )->Clone("h_qcd_METUCPlus");  
  TH1F* h_qcd_METUCMinus = ((TH1F*)qcd->Get("METUCMinus/cutflow") )->Clone("h_qcd_METUCMinus");  
  TH1F* h_qcd_bTagPlus = ((TH1F*)qcd->Get("bTagPlus/cutflow") )->Clone("h_qcd_bTagPlus");  
  TH1F* h_qcd_bTagMinus = ((TH1F*)qcd->Get("bTagMinus/cutflow") )->Clone("h_qcd_bTagMinus");

  h_ttbar_base->Sumw2(); 
  h_wjet_base->Sumw2(); 
  h_zjet_base->Sumw2();
  h_stop_base->Sumw2();
  h_vv_base->Sumw2();
  h_qcd_base->Sumw2();

  TH1F* h_TotalBkg = h_wjet_base->Clone("h_TotalBkg"); 
  h_TotalBkg->Reset(); 
  h_TotalBkg->Add(h_ttbar_base); 
  h_TotalBkg->Add(h_wjet_base); 
  h_TotalBkg->Add(h_wjet_base);
  h_TotalBkg->Add(h_zjet_base);
  h_TotalBkg->Add(h_stop_base);
  h_TotalBkg->Add(h_vv_base);
  h_TotalBkg->Add(h_qcd_base);

  ofstream outFile; 
  outFile.open("summaryWithSyst.tex"); 
   
  outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl;  
  outFile<<"\\pagestyle{empty}"<<endl;  
  outFile<<"\\usepackage{epsfig}"<<endl;  
  outFile<<"\\usepackage{amsmath}"<<endl;  
  outFile<<"\\usepackage{array}"<<endl;  
  outFile<<"\\usepackage{multirow}"<<endl;  
  outFile<<"\\usepackage[cm]{fullpage}"<<endl;  
  outFile<<"\\textheight = 8.in"<<endl;  
  outFile<<"\\textwidth 7.0in"<<endl;  
  outFile<<"\\begin{document}"<<endl;  
  outFile<<"\\begin{center}"<<endl;  
  outFile<<"\\begin{LARGE}"<<endl;  
  outFile<<"\\begin{tabular}{ | c| c| }"<<endl;  
  outFile<<"\\multicolumn{2}{c}{ } \\\\"<<endl;  
  outFile<<"\\hline "<<endl; 
  outFile<<"\\multicolumn{1}{|c|}{Source} & \\multicolumn{1}{|c|}{N$_{\\rm events}$ $\\pm$ MC stat $\\pm$ JES/MET scale $\\pm$ bTag } \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"\\hline "<<endl; 
  outFile<<"HW, $M_{H}=120~GeV/c^{2}$"<<" & "<<h_wh_base->GetBinContent(5)<<" $\\pm$ "<<h_wh_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_wh_JESPlus->GetBinContent(5) - h_wh_base->GetBinContent(5)), fabs(h_wh_base->GetBinContent(5) - h_wh_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wh_JERPlus->GetBinContent(5) - h_wh_base->GetBinContent(5)), fabs(h_wh_base->GetBinContent(5) - h_wh_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wh_METUCPlus->GetBinContent(5) - h_wh_base->GetBinContent(5)), fabs(h_wh_base->GetBinContent(5) - h_wh_METUCMinus->GetBinContent(5))), 2)) <<" $\\pm$ "<<max(fabs(h_wh_bTagPlus->GetBinContent(5) - h_wh_base->GetBinContent(5)), fabs(h_wh_base->GetBinContent(5) - h_wh_bTagMinus->GetBinContent(5))) <<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;  

  outFile<<"SM $t\\bar{t}$"<<" & "<<h_ttbar_base->GetBinContent(5)<<" $\\pm$ "<<h_ttbar_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_ttbar_JESPlus->GetBinContent(5) - h_ttbar_base->GetBinContent(5)), fabs(h_ttbar_base->GetBinContent(5) - h_ttbar_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_ttbar_JERPlus->GetBinContent(5) - h_ttbar_base->GetBinContent(5)), fabs(h_ttbar_base->GetBinContent(5) - h_ttbar_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_ttbar_METUCPlus->GetBinContent(5) - h_ttbar_base->GetBinContent(5)), fabs(h_ttbar_base->GetBinContent(5) - h_ttbar_METUCMinus->GetBinContent(5))), 2)) <<" $\\pm$ "<<max(fabs(h_ttbar_bTagPlus->GetBinContent(5) - h_ttbar_base->GetBinContent(5)), fabs(h_ttbar_base->GetBinContent(5) - h_ttbar_bTagMinus->GetBinContent(5))) <<" \\\\ "<<endl;

  outFile<<"W+Jets"<<" & "<<h_wjet_base->GetBinContent(5)<<" $\\pm$ "<<h_wjet_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_wjet_JESPlus->GetBinContent(5) - h_wjet_base->GetBinContent(5)), fabs(h_wjet_base->GetBinContent(5) - h_wjet_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wjet_JERPlus->GetBinContent(5) - h_wjet_base->GetBinContent(5)), fabs(h_wjet_base->GetBinContent(5) - h_wjet_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wjet_METUCPlus->GetBinContent(5) - h_wjet_base->GetBinContent(5)), fabs(h_wjet_base->GetBinContent(5) - h_wjet_METUCMinus->GetBinContent(5))), 2)) <<" $\\pm$ "<<max(fabs(h_wjet_bTagPlus->GetBinContent(5) - h_wjet_base->GetBinContent(5)), fabs(h_wjet_base->GetBinContent(5) - h_wjet_bTagMinus->GetBinContent(5))) <<" \\\\ "<<endl;  

  outFile<<"Z+Jets"<<" & "<<h_zjet_base->GetBinContent(5)<<" $\\pm$ "<<h_zjet_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_zjet_JESPlus->GetBinContent(5) - h_zjet_base->GetBinContent(5)), fabs(h_zjet_base->GetBinContent(5) - h_zjet_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_zjet_JERPlus->GetBinContent(5) - h_zjet_base->GetBinContent(5)), fabs(h_zjet_base->GetBinContent(5) - h_zjet_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_zjet_METUCPlus->GetBinContent(5) - h_zjet_base->GetBinContent(5)),fabs(h_zjet_base->GetBinContent(5) - h_zjet_METUCMinus->GetBinContent(5))), 2)) <<" $\\pm$ "<<max(fabs(h_zjet_bTagPlus->GetBinContent(5) - h_zjet_base->GetBinContent(5)), fabs(h_zjet_base->GetBinContent(5) - h_zjet_bTagMinus->GetBinContent(5))) <<" \\\\ "<<endl;

  outFile<<"single Top"<<" & "<<h_stop_base->GetBinContent(5)<<" $\\pm$ "<<h_stop_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_stop_JESPlus->GetBinContent(5) - h_stop_base->GetBinContent(5)), fabs(h_stop_base->GetBinContent(5) - h_stop_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_stop_JERPlus->GetBinContent(5) - h_stop_base->GetBinContent(5)), fabs(h_stop_base->GetBinContent(5) - h_stop_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_stop_METUCPlus->GetBinContent(5) - h_stop_base->GetBinContent(5)),fabs(h_stop_base->GetBinContent(5) - h_stop_METUCMinus->GetBinContent(5))), 2)) <<" $\\pm$ "<<max(fabs(h_stop_bTagPlus->GetBinContent(5) - h_stop_base->GetBinContent(5)), fabs(h_stop_base->GetBinContent(5) - h_stop_bTagMinus->GetBinContent(5))) <<" \\\\ "<<endl; 
  
  outFile<<"diboson"<<" & "<<h_vv_base->GetBinContent(5)<<" $\\pm$ "<<h_vv_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_vv_JESPlus->GetBinContent(5) - h_vv_base->GetBinContent(5)), fabs(h_vv_base->GetBinContent(5) - h_vv_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_vv_JERPlus->GetBinContent(5) - h_vv_base->GetBinContent(5)), fabs(h_vv_base->GetBinContent(5) - h_vv_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_vv_METUCPlus->GetBinContent(5) - h_vv_base->GetBinContent(5)),fabs(h_vv_base->GetBinContent(5) - h_vv_METUCMinus->GetBinContent(5))), 2)) <<" $\\pm$ "<<max(fabs(h_vv_bTagPlus->GetBinContent(5) - h_vv_base->GetBinContent(5)), fabs(h_vv_base->GetBinContent(5) - h_vv_bTagMinus->GetBinContent(5))) <<" \\\\ "<<endl;

  outFile<<"QCD"<<" & "<<h_qcd_base->GetBinContent(5)<<" $\\pm$ "<<h_qcd_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_qcd_JESPlus->GetBinContent(5) - h_qcd_base->GetBinContent(5)), fabs(h_qcd_base->GetBinContent(5) - h_qcd_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_qcd_JERPlus->GetBinContent(5) - h_qcd_base->GetBinContent(5)), fabs(h_qcd_base->GetBinContent(5) - h_qcd_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_qcd_METUCPlus->GetBinContent(5) - h_qcd_base->GetBinContent(5)),fabs(h_qcd_base->GetBinContent(5) - h_qcd_METUCMinus->GetBinContent(5))), 2)) <<" $\\pm$ "<<max(fabs(h_qcd_bTagPlus->GetBinContent(5) - h_qcd_base->GetBinContent(5)), fabs(h_qcd_base->GetBinContent(5) - h_qcd_bTagMinus->GetBinContent(5))) <<" \\\\ "<<endl;

  outFile<<"\\hline "<<endl;  
  outFile<<"\\hline "<<endl; 
  outFile<<"Total Bkg"<<" & "<<h_TotalBkg->GetBinContent(5)<<" $\\pm$ "<<h_TotalBkg->GetBinError(5)<<" $\\pm$ "<<" -- "<<" $\\pm$ "<<" -- "<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\hline "<<endl;  
  outFile<<"Data "<<" & "<<h_TotalBkg->GetBinContent(5)<<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;    
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"\\end{center}"<<endl;  
  outFile<<"\\end{document}"<<endl; 
 
  outFile.close(); 
} 

void makeSummaryTable()  
{  

  //declare fixed uncertainties on efficiency and scale factors
  float effUnWH = 0;
  float effUnTTbar = 0;
  float effUnWJets = 0;
  
  TFile *wh = new TFile("wjet_selection.root");  
  TH1F* h_wh_base = ((TH1F*)wh->Get("base/cutflow") )->Clone("h_wh_base");  
  TH1F* h_wh_JESPlus = ((TH1F*)wh->Get("JESPlus/cutflow") )->Clone("h_wh_JESPlus"); 
  TH1F* h_wh_JESMinus = ((TH1F*)wh->Get("JESMinus/cutflow") )->Clone("h_wh_JESMinus"); 
  TH1F* h_wh_JERPlus = ((TH1F*)wh->Get("JERPlus/cutflow") )->Clone("h_wh_JERPlus"); 
  TH1F* h_wh_JERMinus = ((TH1F*)wh->Get("JERMinus/cutflow") )->Clone("h_wh_JERMinus"); 
  TH1F* h_wh_METUCPlus = ((TH1F*)wh->Get("METUCPlus/cutflow") )->Clone("h_wh_METUCPlus"); 
  TH1F* h_wh_METUCMinus = ((TH1F*)wh->Get("METUCMinus/cutflow") )->Clone("h_wh_METUCMinus"); 
  TH1F* h_wh_bTagPlus = ((TH1F*)wh->Get("bTagPlus/cutflow") )->Clone("h_wh_bTagPlus"); 
  TH1F* h_wh_bTagMinus = ((TH1F*)wh->Get("bTagMinus/cutflow") )->Clone("h_wh_bTagMinus"); 
  
  TFile *ttbar = new TFile("wjet_selection.root");  
  TH1F* h_ttbar_base = ((TH1F*)ttbar->Get("base/cutflow") )->Clone("h_ttbar_base");   
  TH1F* h_ttbar_JESPlus = ((TH1F*)ttbar->Get("JESPlus/cutflow") )->Clone("h_ttbar_JESPlus");  
  TH1F* h_ttbar_JESMinus = ((TH1F*)ttbar->Get("JESMinus/cutflow") )->Clone("h_ttbar_JESMinus");  
  TH1F* h_ttbar_JERPlus = ((TH1F*)ttbar->Get("JERPlus/cutflow") )->Clone("h_ttbar_JERPlus");  
  TH1F* h_ttbar_JERMinus = ((TH1F*)ttbar->Get("JERMinus/cutflow") )->Clone("h_ttbar_JERMinus");  
  TH1F* h_ttbar_METUCPlus = ((TH1F*)ttbar->Get("METUCPlus/cutflow") )->Clone("h_ttbar_METUCPlus");  
  TH1F* h_ttbar_METUCMinus = ((TH1F*)ttbar->Get("METUCMinus/cutflow") )->Clone("h_ttbar_METUCMinus");  
  TH1F* h_ttbar_bTagPlus = ((TH1F*)ttbar->Get("bTagPlus/cutflow") )->Clone("h_ttbar_bTagPlus");  
  TH1F* h_ttbar_bTagMinus = ((TH1F*)ttbar->Get("bTagMinus/cutflow") )->Clone("h_ttbar_bTagMinus");  
 
  TFile *wjet = new TFile("wjet_selection.root");  
  TH1F* h_wjet_base = ((TH1F*)wjet->Get("base/cutflow") )->Clone("h_wjet_base");   
  TH1F* h_wjet_JESPlus = ((TH1F*)wjet->Get("JESPlus/cutflow") )->Clone("h_wjet_JESPlus");  
  TH1F* h_wjet_JESMinus = ((TH1F*)wjet->Get("JESMinus/cutflow") )->Clone("h_wjet_JESMinus");  
  TH1F* h_wjet_JERPlus = ((TH1F*)wjet->Get("JERPlus/cutflow") )->Clone("h_wjet_JERPlus");  
  TH1F* h_wjet_JERMinus = ((TH1F*)wjet->Get("JERMinus/cutflow") )->Clone("h_wjet_JERMinus");  
  TH1F* h_wjet_METUCPlus = ((TH1F*)wjet->Get("METUCPlus/cutflow") )->Clone("h_wjet_METUCPlus");  
  TH1F* h_wjet_METUCMinus = ((TH1F*)wjet->Get("METUCMinus/cutflow") )->Clone("h_wjet_METUCMinus");  
  TH1F* h_wjet_bTagPlus = ((TH1F*)wjet->Get("bTagPlus/cutflow") )->Clone("h_wjet_bTagPlus");  
  TH1F* h_wjet_bTagMinus = ((TH1F*)wjet->Get("bTagMinus/cutflow") )->Clone("h_wjet_bTagMinus");  
 
    
  h_ttbar_base->Sumw2();  
  h_wjet_base->Sumw2();  
    
  TH1F* h_TotalBkg = h_wjet->Clone("h_TotalBkg");  
  h_TotalBkg->Reset();  
  h_TotalBkg->Add(h_ttbar_base);  
  h_TotalBkg->Add(h_wjet_base);  
 
  ofstream outFile;  
  outFile.open("summaryWithTotalSyst.tex");  
    
  outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl;   
  outFile<<"\\pagestyle{empty}"<<endl;   
  outFile<<"\\usepackage{epsfig}"<<endl;   
  outFile<<"\\usepackage{amsmath}"<<endl;   
  outFile<<"\\usepackage{array}"<<endl;   
  outFile<<"\\usepackage{multirow}"<<endl;   
  outFile<<"\\usepackage[cm]{fullpage}"<<endl;   
  outFile<<"\\textheight = 8.in"<<endl;   
  outFile<<"\\textwidth 7.0in"<<endl;   
  outFile<<"\\begin{document}"<<endl;   
  outFile<<"\\begin{center}"<<endl;   
  outFile<<"\\begin{LARGE}"<<endl;   
  outFile<<"\\begin{tabular}{ | c| c| }"<<endl;   
  outFile<<"\\multicolumn{2}{c}{ } \\\\"<<endl;   
  outFile<<"\\hline "<<endl;  
  outFile<<"\\multicolumn{1}{|c|}{Source} & \\multicolumn{1}{|c|}{N$_{\rm events}$ $\\pm$ MC stat $\\pm$ syst} \\\\"<<endl;
  outFile<<"\\hline "<<endl; 
  outFile<<"\\hline "<<endl; 
  outFile<<"HW, $M_{H}=120~GeV/c^{2}$"<<" & "<<h_wh_base->GetBinContent(5)<<" $\\pm$ "<<h_wh_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_wh_JESPlus->GetBinContent(5) - h_wh_base->GetBinContent(5)), fabs(h_wh_base->GetBinContent(5) - h_wh_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wh_JERPlus->GetBinContent(5) - h_wh_base->GetBinContent(5)), fabs(h_wh_base->GetBinContent(5) - h_wh_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wh_METUCPlus->GetBinContent(5) - h_wh_base->GetBinContent(5)), fabs(h_wh_base->GetBinContent(5) - h_wh_METUCMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wh_bTagPlus->GetBinContent(5) - h_wh_base->GetBinContent(5)), fabs(h_wh_base->GetBinContent(5) - h_wh_bTagMinus->GetBinContent(5))), 2) + pow(h_wh_base->GetBinContent(5)*effUncWH, 2)) <<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;  
  outFile<<"SM $t\\bar{t}$"<<" & "<<h_ttbar_base->GetBinContent(5)<<" $\\pm$ "<<h_ttbar_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_ttbar_JESPlus->GetBinContent(5) - h_ttbar_base->GetBinContent(5)), fabs(h_ttbar_base->GetBinContent(5) - h_ttbar_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_ttbar_JERPlus->GetBinContent(5) - h_ttbar_base->GetBinContent(5)), fabs(h_ttbar_base->GetBinContent(5) - h_ttbar_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_ttbar_METUCPlus->GetBinContent(5) - h_ttbar_base->GetBinContent(5)), fabs(h_ttbar_base->GetBinContent(5) - h_ttbar_METUCMinus->GetBinContent(5))), 2) + pow(max(fabs(h_ttbar_bTagPlus->GetBinContent(5) - h_ttbar_base->GetBinContent(5)), fabs(h_ttbar_base->GetBinContent(5) - h_ttbar_bTagMinus->GetBinContent(5))), 2) + pow(h_ttbar_base->GetBinContent(5)*effUncTTbar, 2)) <<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;  
  outFile<<"W+Jets"<<" & "<<h_wjet_base->GetBinContent(5)<<" $\\pm$ "<<h_wjet_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_wjet_JESPlus->GetBinContent(5) - h_wjet_base->GetBinContent(5)), fabs(h_wjet_base->GetBinContent(5) - h_wjet_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wjet_JERPlus->GetBinContent(5) - h_wjet_base->GetBinContent(5)), fabs(h_wjet_base->GetBinContent(5) - h_wjet_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wjet_METUCPlus->GetBinContent(5) - h_wjet_base->GetBinContent(5)), fabs(h_wjet_base->GetBinContent(5) - h_wjet_METUCMinus->GetBinContent(5))), 2) + pow(max(fabs(h_wjet_bTagPlus->GetBinContent(5) - h_wjet_base->GetBinContent(5)), fabs(h_wjet_base->GetBinContent(5) - h_wjet_bTagMinus->GetBinContent(5))), 2) + pow(h_wjet_base->GetBinContent(5)*effUncWJets, 2)) <<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;
  outFile<<"\\hline "<<endl;  
  outFile<<"Total Bkg"<<" & "<<h_TotalBkg->GetBinContent(5)<<" $\\pm$ "<<h_TotalBkg->GetBinError(5)<<" $\\pm$ "<<" -- "<<" \\\\ "<<endl;
  outFile<<"\\hline "<<endl;   
  outFile<<"\\hline "<<endl;  
  outFile<<"Data "<<" & "<<h_TotalBkg->GetBinContent(5)<<" \\\\ "<<endl; 
  outFile<<"\\hline "<<endl;    
  outFile<<"\\hline "<<endl;   
  outFile<<"\\end{tabular}"<<endl; 
  outFile<<"\\end{LARGE}"<<endl;  
  outFile<<"\\end{center}"<<endl;  
  outFile<<"\\end{document}"<<endl; 
 
  outFile.close(); 
} 

  
void makeSummaryTableForSignal()  
{  
   
  TString inputFilePath("/afs/cern.ch/user/g/gkole/public/chargedHiggs/8TeV/2012ABC/uncertainty/kinfit4_v10/");  
  vector<int>massPoints;
  massPoints.push_back(80);
  massPoints.push_back(100); 
  massPoints.push_back(120); 
  massPoints.push_back(140);
  massPoints.push_back(150); 
  massPoints.push_back(155); 
  massPoints.push_back(160); 

  vector<TString>inputFiles;
  for(int i = 0; i < massPoints.size(); i++){
    inputFiles.push_back(Form("signal_M_%d_selection.root", massPoints[i])); 
    //cout<<inputFiles[i]<<endl;
  }

  ofstream outFile;  
  outFile.open("summaryTableForSignal.tex");  
    
  outFile<<"\\documentclass[landscape,letterpaper]{article}"<<endl;   
  outFile<<"\\pagestyle{empty}"<<endl;   
  outFile<<"\\usepackage{epsfig}"<<endl;   
  outFile<<"\\usepackage{amsmath}"<<endl;   
  outFile<<"\\usepackage{array}"<<endl;   
  outFile<<"\\usepackage{multirow}"<<endl;   
  outFile<<"\\usepackage[cm]{fullpage}"<<endl;   
  outFile<<"\\textheight = 8.in"<<endl;   
  outFile<<"\\textwidth 7.0in"<<endl;   
  outFile<<"\\begin{document}"<<endl;   
  outFile<<"\\begin{center}"<<endl;   
  outFile<<"\\begin{LARGE}"<<endl;   
  outFile<<"\\begin{tabular}{ | c| c| }"<<endl;   
  outFile<<"\\multicolumn{2}{c}{ } \\\\"<<endl;   
  outFile<<"\\hline "<<endl;  
  outFile<<"\\multicolumn{1}{|c|}{Source} & \\multicolumn{1}{|c|}{N$_{\\rm events}$ $\\pm$ MC stat $\\pm$ JES/MET scale $\\pm$ bTag } \\\\"<<endl; 
  outFile<<"\\hline "<<endl;  
  outFile<<"\\hline "<<endl;  

  for(int i = 0; i < massPoints.size(); i++){

    TFile *f = new TFile(inputFilePath+inputFiles[i]);  
    TH1F* h_base = ((TH1F*)f->Get("base/cutflow") )->Clone("h_base");  
    TH1F* h_JESPlus = ((TH1F*)f->Get("JESPlus/cutflow") )->Clone("h_JESPlus"); 
    TH1F* h_JESMinus = ((TH1F*)f->Get("JESMinus/cutflow") )->Clone("h_JESMinus"); 
    TH1F* h_JERPlus = ((TH1F*)f->Get("JERPlus/cutflow") )->Clone("h_JERPlus"); 
    TH1F* h_JERMinus = ((TH1F*)f->Get("JERMinus/cutflow") )->Clone("h_JERMinus"); 
    TH1F* h_METUCPlus = ((TH1F*)f->Get("METUCPlus/cutflow") )->Clone("h_METUCPlus"); 
    TH1F* h_METUCMinus = ((TH1F*)f->Get("METUCMinus/cutflow") )->Clone("h_METUCMinus"); 
    TH1F* h_bTagPlus = ((TH1F*)f->Get("bTagPlus/cutflow") )->Clone("h_bTagPlus"); 
    TH1F* h_bTagMinus = ((TH1F*)f->Get("bTagMinus/cutflow") )->Clone("h_bTagMinus"); 

    outFile<<"HW, $M_{H}="<<massPoints[i]<<" GeV/c^{2}$"<<" & "<<h_base->GetBinContent(5)<<" $\\pm$ "<<h_base->GetBinError(5)<<" $\\pm$ "<< sqrt(pow(max(fabs(h_JESPlus->GetBinContent(5) - h_base->GetBinContent(5)), fabs(h_base->GetBinContent(5) - h_JESMinus->GetBinContent(5))), 2) + pow(max(fabs(h_JERPlus->GetBinContent(5) - h_base->GetBinContent(5)), fabs(h_base->GetBinContent(5) - h_JERMinus->GetBinContent(5))), 2) + pow(max(fabs(h_METUCPlus->GetBinContent(5) - h_base->GetBinContent(5)), fabs(h_base->GetBinContent(5) - h_METUCMinus->GetBinContent(5))), 2)) <<" $\\pm$ "<<max(fabs(h_bTagPlus->GetBinContent(5) - h_base->GetBinContent(5)), fabs(h_base->GetBinContent(5) - h_bTagMinus->GetBinContent(5))) <<" \\\\ "<<endl;  

  }

  outFile<<"\\hline "<<endl;   
  outFile<<"\\hline "<<endl;  
  outFile<<"\\end{tabular}"<<endl;  
  outFile<<"\\end{LARGE}"<<endl;   
  outFile<<"\\end{center}"<<endl;   
  outFile<<"\\end{document}"<<endl;  

  outFile.close();
}
