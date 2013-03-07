#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
#include <iomanip>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <algorithm> // std::max

void makeHPlusDataCard()
{

  double scale_factor = 1.76877120000000012e+01;
  TFile* fout = new TFile("HplusShapes.root", "RECREATE");

  TString inputFilePath("/afs/cern.ch/user/g/gkole/public/chargedHiggs/8TeV/2012ABC/uncertainty/processedFiles/");
  //TString inputFilePath("./processedFiles/");
  
  TFile* fttbar = TFile::Open(inputFilePath+"ttbar_selection.root", "READ");
  if(fttbar == 0) return;
  if(fttbar->IsZombie()){fttbar->Close(); return;}

  fout->cd();
  TH1F* ttbar = (fttbar->Get("base/mjj_kfit"))->Clone("ttbar");
  ttbar->SetName("ttbar");
  ttbar->Scale(scale_factor);
  ttbar->Write("ttbar");
  TH1F* ttbar_JESUp = (fttbar->Get("JESPlus/mjj_kfit"))->Clone("ttbar_JESUp"); 
  ttbar_JESUp->SetName("ttbar_JESUp");
  ttbar_JESUp->Scale(scale_factor); 
  ttbar_JESUp->Write("ttbar_JESUp");
  TH1F* ttbar_JESDown = (fttbar->Get("JESMinus/mjj_kfit"))->Clone("ttbar_JESDown");  
  ttbar_JESDown->SetName("ttbar_JESDown");
  ttbar_JESDown->Scale(scale_factor);  
  ttbar_JESDown->Write("ttbar_JESUp");
  TH1F* ttbar_bTagUp = (fttbar->Get("bTagPlus/mjj_kfit"))->Clone("ttbar_bTagUp");  
  ttbar_bTagUp->Scale(scale_factor);  
  TH1F* ttbar_bTagDown = (fttbar->Get("bTagMinus/mjj_kfit"))->Clone("ttbar_bTagDown");   
  ttbar_bTagDown->Scale(scale_factor);   

  
  TH1F* ttll = ttbar->Clone("ttll"); 
  ttll->Scale(0);
  //ttll->Scale(scale_factor); 
  ttll->Write(); 
  TH1F* ttll_JESUp = ttbar_JESUp->Clone("ttll_JESUp");  
  ttll_JESUp->Scale(0);
  //ttll_JESUp->Scale(scale_factor);  
  ttll_JESUp->Write(); 
  TH1F* ttll_JESDown = ttbar_JESDown->Clone("ttll_JESDown");   
  ttll_JESDown->Scale(0);
  //ttll_JESDown->Scale(scale_factor);   
  ttll_JESDown->Write(); 

  TFile* fwjet = TFile::Open(inputFilePath+"wjet_selection.root", "READ"); 
  if(fwjet == 0) return; 
  if(fwjet->IsZombie()){fwjet->Close(); return;} 

  fout->cd();
  TH1F* wjet = (fwjet->Get("base/mjj_kfit"))->Clone("wjet"); 
  wjet->Scale(scale_factor); 
  wjet->Write(); 
  TH1F* wjet_JESUp = (fwjet->Get("JESPlus/mjj_kfit"))->Clone("wjet_JESUp");  
  wjet_JESUp->Scale(scale_factor);  
  wjet_JESUp->Write(); 
  TH1F* wjet_JESDown = (fwjet->Get("JESMinus/mjj_kfit"))->Clone("wjet_JESDown");   
  wjet_JESDown->Scale(scale_factor);   
  wjet_JESDown->Write(); 
  TH1F* wjet_bTagUp = (fwjet->Get("bTagPlus/mjj_kfit"))->Clone("wjet_bTagUp");   
  wjet_bTagUp->Scale(scale_factor);   
  TH1F* wjet_bTagDown = (fwjet->Get("bTagMinus/mjj_kfit"))->Clone("wjet_bTagDown");    
  wjet_bTagDown->Scale(scale_factor);    

  TFile* fzjet = TFile::Open(inputFilePath+"zjet_selection.root", "READ"); 
  if(fzjet == 0) return; 
  if(fzjet->IsZombie()){fzjet->Close(); return;} 

  fout->cd();
  TH1F* zjet = (fzjet->Get("base/mjj_kfit"))->Clone("zjet"); 
  zjet->Scale(scale_factor); 
  zjet->Write(); 
  TH1F* zjet_JESUp = (fzjet->Get("JESPlus/mjj_kfit"))->Clone("zjet_JESUp");  
  zjet_JESUp->Scale(scale_factor);  
  zjet_JESUp->Write(); 
  TH1F* zjet_JESDown = (fzjet->Get("JESMinus/mjj_kfit"))->Clone("zjet_JESDown");   
  zjet_JESDown->Scale(scale_factor);   
  zjet_JESDown->Write(); 
  TH1F* zjet_bTagUp = (fzjet->Get("bTagPlus/mjj_kfit"))->Clone("zjet_bTagUp");   
  zjet_bTagUp->Scale(scale_factor);   
  TH1F* zjet_bTagDown = (fzjet->Get("bTagMinus/mjj_kfit"))->Clone("zjet_bTagDown");    
  zjet_bTagDown->Scale(scale_factor);    


  TFile* fstop = TFile::Open(inputFilePath+"singletop_selection.root", "READ"); 
  if(fstop == 0) return; 
  if(fstop->IsZombie()){fstop->Close(); return;} 
 
  fout->cd();
  TH1F* stop = (fstop->Get("base/mjj_kfit"))->Clone("stop"); 
  stop->Scale(scale_factor); 
  stop->Write(); 
  TH1F* stop_JESUp = (fstop->Get("JESPlus/mjj_kfit"))->Clone("stop_JESUp");  
  stop_JESUp->Scale(scale_factor);  
  stop_JESUp->Write(); 
  TH1F* stop_JESDown = (fstop->Get("JESMinus/mjj_kfit"))->Clone("stop_JESDown");   
  stop_JESDown->Scale(scale_factor);   
  stop_JESDown->Write(); 
  TH1F* stop_bTagUp = (fstop->Get("bTagPlus/mjj_kfit"))->Clone("stop_bTagUp");   
  stop_bTagUp->Scale(scale_factor);   
  TH1F* stop_bTagDown = (fstop->Get("bTagMinus/mjj_kfit"))->Clone("stop_bTagDown");    
  stop_bTagDown->Scale(scale_factor);    


  TFile* fdiboson = TFile::Open(inputFilePath+"diboson_selection.root", "READ"); 
  if(fdiboson == 0) return; 
  if(fdiboson->IsZombie()){fdiboson->Close(); return;} 

  fout->cd();
  TH1F* diboson = (fdiboson->Get("base/mjj_kfit"))->Clone("diboson"); 
  diboson->Scale(scale_factor); 
  diboson->Write(); 
  TH1F* diboson_JESUp = (fdiboson->Get("JESPlus/mjj_kfit"))->Clone("diboson_JESUp");  
  diboson_JESUp->Scale(scale_factor);  
  diboson_JESUp->Write(); 
  TH1F* diboson_JESDown = (fdiboson->Get("JESMinus/mjj_kfit"))->Clone("diboson_JESDown");   
  diboson_JESDown->Scale(scale_factor);   
  diboson_JESDown->Write(); 
  TH1F* diboson_bTagUp = (fdiboson->Get("bTagPlus/mjj_kfit"))->Clone("diboson_bTagUp");   
  diboson_bTagUp->Scale(scale_factor);   
  TH1F* diboson_bTagDown = (fdiboson->Get("bTagMinus/mjj_kfit"))->Clone("diboson_bTagDown");    
  diboson_bTagDown->Scale(scale_factor);    


  TFile* fqcd = TFile::Open(inputFilePath+"qcd_selection.root", "READ"); 
  if(fqcd == 0) return; 
  if(fqcd->IsZombie()){fqcd->Close(); return;} 
 
  fout->cd();
  TH1F* qcd = (fqcd->Get("base/mjj_kfit"))->Clone("qcd"); 
  qcd->Scale(scale_factor); 
  qcd->Write(); 
  TH1F* qcd_JESUp = (fqcd->Get("JESPlus/mjj_kfit"))->Clone("qcd_JESUp");  
  qcd_JESUp->Scale(scale_factor);  
  qcd_JESUp->Write(); 
  TH1F* qcd_JESDown = (fqcd->Get("JESMinus/mjj_kfit"))->Clone("qcd_JESDown");   
  qcd_JESDown->Scale(scale_factor);   
  qcd_JESDown->Write(); 
  TH1F* qcd_bTagUp = (fqcd->Get("bTagPlus/mjj_kfit"))->Clone("qcd_bTagUp");   
  qcd_bTagUp->Scale(scale_factor);   
  TH1F* qcd_bTagDown = (fqcd->Get("bTagMinus/mjj_kfit"))->Clone("qcd_bTagDown");    
  qcd_bTagDown->Scale(scale_factor);    


  TFile* fdata = TFile::Open(inputFilePath+"data_selection.root", "READ");
  if(fdata == 0) return;
  if(fdata->IsZombie()){fdata->Close(); return;}
  fout->cd();
  TH1F* data_obs = (fdata->Get("base/mjj_kfit"))->Clone("data_obs");
  data_obs->Write("data_obs");
  
  
  vector<int>massPoints;
  massPoints.push_back(80);
  massPoints.push_back(100);
  massPoints.push_back(120);
  massPoints.push_back(140);
  massPoints.push_back(150);
  massPoints.push_back(155);
  massPoints.push_back(160);

  for(int i = 0; i < massPoints.size(); i++){
    cout<<" mass point "<<massPoints[i]<<endl;
    TFile* fwh = TFile::Open(inputFilePath+Form("wh_M_%d_selection.root", massPoints[i]), "READ");
    if(fwh == 0) return;
    if(fwh->IsZombie()){fwh->Close(); return;}
 
    fout->cd();
    TH1F* wh = (TH1F*)(fwh->Get("base/mjj_kfit"))->Clone(); 
    wh->SetName(Form("WH%d",massPoints[i]));
    wh->Scale(scale_factor); 
    wh->Write(Form("WH%d",massPoints[i])); 
    TH1F* wh_JESUp = (TH1F*)(fwh->Get("JESPlus/mjj_kfit"))->Clone();  
    wh_JESUp->SetName(Form("WH%d_JESUp",massPoints[i]));
    wh_JESUp->Scale(scale_factor);  
    wh_JESUp->Write(Form("WH%d_JESUp",massPoints[i])); 
    TH1F* wh_JESDown = (TH1F*)(fwh->Get("JESMinus/mjj_kfit"))->Clone();
    wh_JESDown->SetName(Form("WH%d_JESDown",massPoints[i]));
    wh_JESDown->Scale(scale_factor);   
    wh_JESDown->Write(Form("WH%d_JESDown",massPoints[i])); 
    TH1F* wh_bTagUp = (TH1F*)(fwh->Get("bTagPlus/mjj_kfit"))->Clone(Form("WH%d_bTagUp",massPoints[i]));   
    wh_bTagUp->Scale(scale_factor);   
    TH1F* wh_bTagDown = (TH1F*)(fwh->Get("bTagMinus/mjj_kfit"))->Clone(Form("WH%d_bTagDown",massPoints[i])); 
    wh_bTagDown->Scale(scale_factor);    

    //Temporary arrangement for HH, not fully correct
    TH1F* hh = (TH1F*)wh->Clone();
    hh->Scale(1./hh->Integral());
    hh->SetName(Form("HH%d",massPoints[i])); 
    hh->Write(Form("HH%d",massPoints[i]));  
    TH1F* hh_JESUp = (TH1F*)wh_JESUp->Clone();
    hh_JESUp->Scale(1./hh_JESUp->Integral());
    hh_JESUp->SetName(Form("HH%d_JESUp",massPoints[i])); 
    hh_JESUp->Write(Form("HH%d_JESUp",massPoints[i]));  
    TH1F* hh_JESDown = (TH1F*)wh_JESDown->Clone();
    hh_JESDown->Scale(1./hh_JESDown->Integral());
    hh_JESDown->SetName(Form("HH%d_JESDown",massPoints[i])); 
    hh_JESDown->Write(Form("HH%d_JESDown",massPoints[i]));  
    TH1F* hh_bTagUp = (TH1F*)wh_bTagUp->Clone(Form("HH%d_bTagUp",massPoints[i])); 
    hh_bTagUp->Scale(1./hh_bTagUp->Integral()); 
    TH1F* hh_bTagDown = (TH1F*)wh_bTagDown->Clone(Form("HH%d_bTagDown",massPoints[i])); 
    hh_bTagDown->Scale(1./hh_bTagDown->Integral()); 

    /*
    //This should be used when files are available
    TFile* fhh = TFile::Open(inputFilePath+Form("hh_M_%d_selection.root", massPoints[i]), "READ"); 
    if(fhh == 0) return; 
    if(fhh->IsZombie()){fhh->Close(); return;} 
    
    fout->cd(); 
    TH1F* hh = (fhh->Get("base/mjj_kfit"))->Clone();  
    hh->SetName(Form("HH%d",massPoints[i])); 
    hh->Scale(scale_factor);  
    hh->Write(Form("HH%d",massPoints[i]));  
    TH1F* hh_JESUp = (fhh->Get("JESPlus/mjj_kfit"))->Clone();   
    hh_JESUp->SetName(Form("HH%d_JESUp",massPoints[i])); 
    hh_JESUp->Scale(scale_factor);   
    hh_JESUp->Write(Form("HH%d_JESUp",massPoints[i]));  
    TH1F* hh_JESDown = (fhh->Get("JESMinus/mjj_kfit"))->Clone(); 
    hh_JESDown->SetName(Form("HH%d_JESDown",massPoints[i])); 
    hh_JESDown->Scale(scale_factor);    
    hh_JESDown->Write(Form("HH%d_JESDown",massPoints[i]));  
    TH1F* hh_bTagUp = (fhh->Get("bTagPlus/mjj_kfit"))->Clone(Form("HH%d_bTagUp",massPoints[i]));    
    hh_bTagUp->Scale(scale_factor);    
    TH1F* hh_bTagDown = (fhh->Get("bTagMinus/mjj_kfit"))->Clone(Form("HH%d_bTagDown",massPoints[i]));  
    hh_bTagDown->Scale(scale_factor);     
    */

    //Make Data card for each mass
    ifstream in;
    char* c = new char[1000];
    in.open("template/datacard_csbar.txt");
    
    ofstream out(Form("datacard_csbar_mH%d.txt", massPoints[i]));
    out.precision(8);

    time_t secs=time(0);
    tm *t=localtime(&secs);
    
    while (in.good())
      {
	in.getline(c,1000,'\n');
	if (in.good()){
	  string line(c);
	  
	  if(line.find("Date")!=string::npos){
	    string day = string(Form("%d",t->tm_mday));
	    string month = string(Form("%d",t->tm_mon+1));
	    string year = string(Form("%d",t->tm_year+1900));
	    line.replace( line.find("XXX") , 3 , day+"/"+month+"/"+year);
	    out << line << endl;
	  }
	  else if(line.find("Description")!=string::npos){
	    line.replace( line.find("YYY") , 3 , string(Form("%d", massPoints[i])) );
	    line.replace( line.find("ZZZ") , 3 , string(Form("%f", scale_factor)) ); 
	    out << line << endl;
	  }
	  else if(line.find("shapes")!=string::npos){
	    line.replace( line.find("XXX") , 3 , string("HplusShapes")  );
	    out << line << endl;
	  }
	  else if(line.find("Observation")!=string::npos){
	    line.replace( line.find("XXX") , 3 , string(Form("%.0f", data_obs->Integral()) )  );
	    out << line << endl;
	  }
	  else if(line.find("process")!=string::npos && line.find("WH")!=string::npos){
	    line.replace( line.find("XXX") , 3 , string(Form("%d", massPoints[i])) );
	    line.replace( line.find("YYY") , 3 , string(Form("%d", massPoints[i])) );
	    out << line << endl;
	  }
	  else if(line.find("rate")!=string::npos){
	    string rate = "rate      ";  
	    string space = "      ";
	    out << rate ;
	    out << space << hh->Integral()
		<< space << wh->Integral()
		<< space << ttbar->Integral()
		<< space << ttll->Integral()
		<< space << wjet->Integral()
		<< space << zjet->Integral()
		<< space << qcd->Integral()
		<< space << stop->Integral()
		<< space << diboson->Integral()
		<< space << "Projected event rates"<<endl;
	  }
	  else if(line.find("CMS_eff_b")!=string::npos){
	    float bTagUnc_hh = (hh->Integral() > 0) ? 1 + (max(fabs(hh_bTagUp->Integral() - hh->Integral()), fabs(hh->Integral() - hh_bTagDown->Integral()))/hh->Integral()) : 1.00;
	    line.replace( line.find("XXXX") , 4 , string(Form("%.2f", bTagUnc_hh)) );
	    	    
	    float bTagUnc_wh = (wh->Integral() > 0) ? 1 + (max(fabs(wh_bTagUp->Integral() - wh->Integral()), fabs(wh->Integral() - wh_bTagDown->Integral()))/wh->Integral()) : 1.00; 
            line.replace( line.find("YYYY") , 4 , string(Form("%.2f", bTagUnc_wh)) ); 

	    float bTagUnc_ttbar = (ttbar->Integral() > 0) ? 1 + (max(fabs(ttbar_bTagUp->Integral() - ttbar->Integral()), fabs(ttbar->Integral() - ttbar_bTagDown->Integral()))/ttbar->Integral()) : 1.00; 
            line.replace( line.find("ZZZZ") , 4 , string(Form("%.2f", bTagUnc_ttbar)) ); 

	    float bTagUnc_stop = (stop->Integral() > 0) ? 1 + (max(fabs(stop_bTagUp->Integral() - stop->Integral()), fabs(stop->Integral() - stop_bTagDown->Integral()))/stop->Integral()) : 1.00; 
            line.replace( line.find("KKKK") , 4 , string(Form("%.2f", bTagUnc_stop)) ); 

	    out << line << endl;
	  }
	  else if(line.find("CMS_mistag_b")!=string::npos){
	    float bTagUnc_wjet = (wjet->Integral() > 0) ? 1 + (max(fabs(wjet_bTagUp->Integral() - wjet->Integral()), fabs(wjet->Integral() - wjet_bTagDown->Integral()))/wjet->Integral()) : 1.00; 
            line.replace( line.find("XXXX") , 4 , string(Form("%.2f", bTagUnc_wjet)) ); 

	    float bTagUnc_zjet = (zjet->Integral() > 0) ? 1 + (max(fabs(zjet_bTagUp->Integral() - zjet->Integral()), fabs(zjet->Integral() - zjet_bTagDown->Integral()))/zjet->Integral()) : 1.00; 
            line.replace( line.find("YYYY") , 4 , string(Form("%.2f", bTagUnc_zjet)) ); 

	    float bTagUnc_qcd = (qcd->Integral() > 0) ? 1 + (max(fabs(qcd_bTagUp->Integral() - qcd->Integral()), fabs(qcd->Integral() - qcd_bTagDown->Integral()))/qcd->Integral()) : 1.00; 
            line.replace( line.find("ZZZZ") , 4 , string(Form("%.2f", bTagUnc_qcd)) ); 

	    float bTagUnc_diboson = (diboson->Integral() > 0) ? 1 + (max(fabs(diboson_bTagUp->Integral() - diboson->Integral()), fabs(diboson->Integral() - diboson_bTagDown->Integral()))/diboson->Integral()) : 1.00; 
            line.replace( line.find("KKKK") , 4 , string(Form("%.2f", bTagUnc_diboson)) ); 

	    out << line << endl;
	  }
	  else if(line.find("CMS_stat_hh")!=string::npos){
	    Double_t sError = 0;
	    Double_t norm = hh->IntegralAndError(1, hh->GetNbinsX(), sError);
	    float statUnc_hh = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00; 
            line.replace( line.find("XXXX") , 4 , string(Form("%.2f", statUnc_hh)) ); 
	    out << line << endl;
	  }
	  else if(line.find("CMS_stat_wh")!=string::npos){ 
            Double_t sError = 0; 
            Double_t norm = wh->IntegralAndError(1, wh->GetNbinsX(), sError); 
            float statUnc_wh = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;  
            line.replace( line.find("XXXX") , 4 , string(Form("%.2f", statUnc_wh)) );  
	    out << line << endl;
          } 
	  else if(line.find("CMS_stat_tt")!=string::npos){  
            Double_t sError = 0;  
            Double_t norm = ttbar->IntegralAndError(1, ttbar->GetNbinsX(), sError);  
            float statUnc_ttbar = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;   
            line.replace( line.find("XXXX") , 4 , string(Form("%.2f", statUnc_ttbar)) );   
	    out << line << endl;
          }  
	  else if(line.find("CMS_stat_wjet")!=string::npos){  
            Double_t sError = 0;  
            Double_t norm = wjet->IntegralAndError(1, wjet->GetNbinsX(), sError);  
            float statUnc_wjet = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;   
            line.replace( line.find("XXXX") , 4 , string(Form("%.2f", statUnc_wjet)) );   
	    out << line << endl;
          }  
	  else if(line.find("CMS_stat_zjet")!=string::npos){ 
            Double_t sError = 0; 
            Double_t norm = zjet->IntegralAndError(1, zjet->GetNbinsX(), sError); 
            float statUnc_zjet = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;  
            line.replace( line.find("XXXX") , 4 , string(Form("%.2f", statUnc_zjet)) );  
            out << line << endl; 
          } 
	  else if(line.find("CMS_stat_stop")!=string::npos){ 
            Double_t sError = 0; 
            Double_t norm = stop->IntegralAndError(1, stop->GetNbinsX(), sError); 
            float statUnc_stop = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;  
            line.replace( line.find("XXXX") , 4 , string(Form("%.2f", statUnc_stop)) );  
            out << line << endl; 
          } 
	  else if(line.find("CMS_stat_qcd")!=string::npos){  
            Double_t sError = 0;  
            Double_t norm = qcd->IntegralAndError(1, qcd->GetNbinsX(), sError);  
            float statUnc_qcd = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;   
            line.replace( line.find("XXXX") , 4 , string(Form("%.2f", statUnc_qcd)) );   
            out << line << endl;  
          }  
	  else if(line.find("CMS_stat_vv")!=string::npos){  
            Double_t sError = 0;  
            Double_t norm = diboson->IntegralAndError(1, diboson->GetNbinsX(), sError);  
            float statUnc_diboson = (norm > 0) ? 1 + (fabs(sError)/norm) : 1.00;   
            line.replace( line.find("XXXX") , 4 , string(Form("%.2f", statUnc_diboson)) );   
            out << line << endl;  
          }  
	  else{ //default without changes
	    out << line << endl;
	  }
	}
	else{
	  out << c << endl;
	}
      }
	
    out.close();
    in.close();
  }
  fout->Close();
}

