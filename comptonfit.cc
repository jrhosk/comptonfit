// Standard C++ libraries
#include <cstdlib>
#include <cstdio>
#include <cstring>
#include <iostream>

// Custom libraries
#include "ComptonConfig.hh"
#include "DataFile.hh"
#include "FontColor.hh"

// Boost libraries
#include "boost/filesystem.hpp"

// Root libraries
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TF1.h"
#include "TFitResultPtr.h"
#include "TVirtualFitter.h"
#include "TMath.h"
#include "Math/MinimizerOptions.h"

int main(int argc, char *argv[])
{

  ComptonConfig *compton = new ComptonConfig();
  compton->GetOptions(argv);
  if(compton->fGraphicsShow) compton->InitGraphicsEngine(argc, argv); 

  compton->GetConfigurationFiles("config/", ".config");
  compton->ReadConfiguration();
  compton->CalculateDriftDistance();

  DataFile *asym_data = new DataFile("data/asymmetry.dat");
  DataFile *yield_data = new DataFile("data/yield.dat");
  DataFile *background_data = new DataFile("data/background.dat");

  asym_data->ReadFile();
  yield_data->ReadFile();
  background_data->ReadFile();

  double initCE = 0.0;

  bool kFoundCE = false;
  bool polSign = false;

  for(int i = 0; i < (int)(asym_data->GetArraySize()); i++) {
    if( (yield_data->GetY().at(i)-background_data->GetY().at(i) < 100.0*background_data->GetYError().at(i) ) && i > 25) {

      initCE = asym_data->GetX().at(i-1); 
      std::cout << "Initial estimate of Compton edge (CE) : "<< initCE << std::endl;
      polSign = asym_data->GetY().at(i-1) > 0 ? 1 : 0;
      kFoundCE = true;
      break;
    }
  }

  if(kFoundCE) {
    if(polSign) std::cout << "Sign of polarization estimated to be positive.\n" << std::endl;
    else std::cout << "Sign of polarization estimated to be negative.\n" << std::endl;
  } else std::cout << red << "**** Alert: Could not determine compton edge. " << white << std::endl;

  double xCedge = compton->RhoToPositionMap(initCE);

  TCanvas *canvas = new TCanvas("canvas", "Asymmetry in run", 10, 10, 1500, 800);
  canvas->SetGridx(1);

  if(compton->fResiduals){
    canvas->Divide(1,2);
    canvas->cd(1);
    canvas->GetPad(1)->SetGridx(1);
  }
  else{
    canvas->cd();
    canvas->SetGridx(1);
  }

  std::vector<double> zero;
  zero.resize(asym_data->GetArraySize());
  std::fill(zero.begin(), zero.end(), 0);

  TGraphErrors *graph = new TGraphErrors(asym_data->GetArraySize(), (Double_t *)(asym_data->GetX().data()), 
					      (Double_t *)(asym_data->GetY().data()), (Double_t *)(zero.data()), (Double_t *)(asym_data->GetYError().data())); 

  graph->GetXaxis()->SetTitle("Strip number");
  graph->GetYaxis()->SetTitle("Asymmetry");
  graph->SetTitle();
  graph->SetMarkerStyle(kFullCircle);
  graph->SetLineColor(kRed);
  graph->SetMarkerColor(kRed);
  graph->SetMaximum(0.5);
  graph->SetMinimum(-0.5);
  graph->GetXaxis()->SetLimits(0,150);

  graph->Draw("AP");

  TF1 *polFit = new TF1("polFit", compton, &ComptonConfig::TheoreticalAsymFit, 1, initCE, 2);

  if(polSign) {                                                                                                                                       
    polFit->SetParameters(initCE, 1.0);   
    polFit->SetParLimits(0, 110.0, 124.0);
    polFit->SetParLimits(1, 0.90, 1.0);   
  } else { 
    polFit->SetParameters(initCE, -1.0);  
    polFit->SetParLimits(0, 110.0, 124.0);
    polFit->SetParLimits(1, -1.0, -0.90);
  }

  polFit->SetParNames("comptonEdge","polarization");
  polFit->SetLineColor(kBlack);

  std::cout << blue << "The maxdist used: " << xCedge << white << std::endl;
  TVirtualFitter::SetMaxIterations(10000);

  int numbIterations = 0;

  TFitResultPtr result = graph->Fit(polFit,"RS 0");

  std::cout << green << "the polarization fit status " << int(result) << white << std::endl;

  if(int(result) == 0) {
    std::cout << blue 
	      << "pol%: " << (polFit->GetParameter(1))*100 << " +/- " << (polFit->GetParError(1))*100 
	      << "\t CE: " << polFit->GetParameter(0) << " +/- " << polFit->GetParError(0) 
	      << white << std::endl;

    std::cout << blue << "Fit converged, hence applying MINOS for better errors" << white << std::endl;

    result = graph->Fit(polFit,"RES 0");
    polFit = graph->GetFunction("polFit");
    numbIterations++;

    if (int(result) !=0 ) {
      std::cout << red << "Failed MINOS though the initial fit was successful" << std::endl;
      std::cout << "Do not update the preceeding fit parameters" << white << std::endl;

    } else {
      std::cout << green << "Updated the fit parameters with the MINOS fit results" << white << std::endl;
      polFit = graph->GetFunction("polFit");
    }
  }

  int maxIterations = 4;
  int iterations = 0;

  while (int(result) != 0) {
    iterations++;

    std::cout << blue << "Repeating attempt to fit without MINOS, attempt# " << numbIterations << white << std::endl;
    ROOT::Math::MinimizerOptions::SetDefaultMinimizer("Minuit", "Scan");
    result = graph->Fit("polFit","RS 0");
    polFit = graph->GetFunction("polFit");

    if(iterations > maxIterations) {
      std::cout << red << "Maximum no.of iterations (" << maxIterations << ") attempted. Exiting." << white << std::endl;
      break;
    }
  }

  if((std::fabs((polFit->GetParameter(0)) - initCE) > 0.1) && (iterations < maxIterations)) { 
    if((polFit->GetParameter(0)) != 0.0) initCE = polFit->GetParameter(0); 

    xCedge = compton->RhoToPositionMap(initCE);
    result = graph->Fit("polFit","SR 0");

    std::cout << green << "Refitting # " << iterations << ", the polarization fit status is: " << int(result) << white << std::endl;
    if(int(result) == 0) {                                             
      std::cout << blue << "Updating the function parameters" << white << std::endl;
      iterations++;
      polFit  = graph->GetFunction("polFit");                        
    }
    std::cout << blue 
	      << "Polariztion: "<< 100*(polFit->GetParameter(1)) << " +/- " << 100*(polFit->GetParError(1)) 
	      << "\tCompton edge: " << polFit->GetParameter(0) << " +/- " << polFit->GetParError(0) 
	      << "\tChi^2/NDF: " << (polFit->GetChisquare())/(polFit->GetNDF()) 
	      << white << std::endl;
  }

  polFit->DrawCopy("same");

  if(compton->fResiduals){
    canvas->cd(2);
    canvas->GetPad(2)->SetGridx(2);

    std::vector <double> residuals;
    
    for(int i = 1; i < (asym_data->GetArraySize()); i++){
      residuals.push_back(asym_data->GetY()[i] - polFit->Eval(asym_data->GetX()[i]));
    }
    TGraphErrors *residual = new TGraphErrors(asym_data->GetArraySize(), (Double_t *)(asym_data->GetX().data()), 
					      (Double_t *)(residuals.data()), (Double_t *)(zero.data()), (Double_t *)(zero.data())); 

    residual->GetXaxis()->SetTitle("Strip number");
    residual->GetYaxis()->SetTitle("Residual");
    residual->SetTitle();
    residual->SetMarkerStyle(kFullCircle);
    residual->SetMarkerColor(kRed);
    residual->SetLineColor(kRed);
    residual->SetMaximum(0.05);
    residual->SetMinimum(-0.05);
    residual->GetXaxis()->SetLimits(0,150);
    residual->SetFillColor(38);
    residual->Draw("AB");

    canvas->Update();
    canvas->Modified();
  } 
  canvas->SaveAs("plots/asymmetry_fit.C");
  
  if(compton->fGraphicsShow){
    compton->RunGraphicsEngine(); 
  }
  return 0;
}
