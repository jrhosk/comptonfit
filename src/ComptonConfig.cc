#define compton_config_cxx
// C++ Libraries

#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <cmath>

using namespace std;

// Root libraries

#include "TAxis.h"
#include "TF1.h"
#include "TGraph.h"
#include "TString.h"
#include "TMath.h"
#include "TApplication.h"

// Boost libraries

#include "boost/filesystem.hpp"
#include "boost/math/constants/constants.hpp"

// Custom Libraries

#include "ComptonConfig.hh"
#include "FontColor.hh"
#include "PhysicalConstants.hh"

#ifdef compton_config_cxx

ComptonConfig::ComptonConfig()
{
  fLowerLimit = 1.0;
  fUpperLimit = 120.0;
  fGraphicsShow = false;
  fResiduals = false; 
  fLowerPolLimit = 0.8;
  fUpperPolLimit = 1.0;
  fLowerCELimit = 115.0;
  fUpperCELimit = 120.0;

}

ComptonConfig::~ComptonConfig()
{
  // Filler
}

double ComptonConfig::ComputeAlpha()
{

  return(1/(1+4*beam.laser_energy*beam.beam_energy/pow(pchep::electron_mass_GeV,2)));

}

double ComptonConfig::CrossSectionFit(double *thisStrip, double *parCx)
{
  double a_const = ComputeAlpha();

  double xStrip = xCedge - (parCx[1] - (*thisStrip))*(det.width+det.spacing)*parCx[0];
  double rhoStrip = (param[0]+ xStrip*param[1]+ xStrip*xStrip*param[2]+ xStrip*xStrip*xStrip*param[3]);
  double rhoPlus = 1.0-rhoStrip*(1.0+a_const);
  double rhoMinus = 1.0-rhoStrip*(1.0 - a_const);

  double dsdrho1 = rhoPlus/rhoMinus;

  return (parCx[2]*((rhoStrip*(1.0 - a_const)*rhoStrip*(1.0 - a_const)/rhoMinus)+1.0+dsdrho1*dsdrho1));
}

double ComptonConfig::TheoreticalAsymFit(double *thisStrip, double *par)
{
  
  double a_const = ComputeAlpha();
  double delta;
  double xStrip = xCedge - (par[0] -(*thisStrip))*(det.width+det.spacing); 
  double rhoStrip = param[0]+ xStrip*param[1]+ xStrip*xStrip*param[2]+ xStrip*xStrip*xStrip*param[3];


  double eGamma   = rhoStrip* kprimemax;
  double betaBar  = TMath::Sqrt(1.0 - (pchep::electron_mass_GeV/beam.beam_energy)*(pchep::electron_mass_GeV/beam.beam_energy));
  double betaCM   = (betaBar*beam.beam_energy-beam.laser_energy)/(beam.beam_energy+beam.laser_energy);   
  double sEq211   = pchep::electron_mass_GeV*pchep::electron_mass_GeV + 2*beam.laser_energy*beam.beam_energy*(1.0+betaCM);       
  double gammaCM  = (beam.beam_energy+beam.laser_energy)/TMath::Sqrt(sEq211);
  double eLaserCM = beam.laser_energy*gammaCM*(1.0+betaCM);
  double eBeamCM  = TMath::Sqrt(beam.laser_energy*beam.laser_energy + std::pow(pchep::electron_mass_GeV, 2));
  double eBetaCM  = eLaserCM/eBeamCM;   
  double costhcm  = (gammaCM*eLaserCM - eGamma)/(eLaserCM*gammaCM*betaCM);

  if(rhoStrip <= 1.0){
    delta=((pchep::fine_struct)/(boost::math::constants::pi<double>()))*(3.0*costhcm-1.0)/(4*(eBetaCM + costhcm));
  }
  else {
    delta = 0.0;
  }

  double radCor = (1.0+delta);

  double rhoPlus = 1.0-rhoStrip*(1.0+a_const);
  double rhoMinus = 1.0-rhoStrip*(1.0 - a_const); 
  double dsdrho1 = rhoPlus/rhoMinus;
  double dsdrho =((rhoStrip*(1.0 - a_const)*rhoStrip*(1.0 - a_const)/rhoMinus)+1.0+dsdrho1*dsdrho1); 

  // radCor = 1;

  return (radCor*(par[1]*(rhoPlus*(1.0-1.0/(rhoMinus*rhoMinus)))/dsdrho));

}

void ComptonConfig::CalculateDriftDistance()
{

  // Calculates the first and second drift distances

  geo.drift_first = std::abs(mag[0].z-mag[1].z) 
    - std::abs(0.5*mag[0].length*std::cos(mag[0].theta)) 
    - std::abs(0.5*mag[1].length*std::cos(mag[1].theta));     // drift distance between third adn fourth dipoles

  det.position = (mag[1].z - 0.5*mag[1].length*std::cos(mag[1].theta)) - det.z;
  geo.drift_second = det.position - det.ce_to_det_bot*std::sin(det.theta);

  std::cout << geo.drift_first << "\t" << geo.drift_second << "\t" << det.position << std::endl;

}

double ComptonConfig::RhoToPositionMap(double initCE)
{

  std::cout << "\nCalculating rho to position mapping .... " << initCE << "\n" << std::endl;

  bool debug = true;

  const int nPoints = 10000;

  double xPrime[nPoints] = {0.0};
  double rho[nPoints] = {0.0};
  double dsdx[nPoints] = {};
  double asym[nPoints] = {};
  double dsdx_0[nPoints] = {};

  ofstream QEDasym;

  double det_angle = det.theta*boost::math::constants::pi<double>()/180;                 //(radians)
  double a_const = ComputeAlpha();

  det.ce_to_det_bot = (initCE)*(det.width + det.spacing);                         // Used the normal method but had to change what I assume was strip spacing (not documented) 

  if(debug) {
    std::cout << blue << "Using the following parameters for rhoToX conversion" << std::endl;
    std::cout << Form("det.ce_to_det_bot: %f = (%f)*(%f + %f)\n", det.ce_to_det_bot, initCE, det.width, det.spacing);
    std::cout << "Beam energy: " << beam.beam_energy << " +/- " << 0.1 << std::endl;
    std::cout << "Field(s):    " << mag[0].dipole << " " << mag[1].dipole << std::endl;
    std::cout << Form("geo.drift_first: %f ", geo.drift_first) << std::endl;
    std::cout << Form("geo.drift_second: %f = %f - %f * %f", geo.drift_second, det.position, det.ce_to_det_bot, sin(det_angle)) << white << std::endl;
  }

  kprimemax = 4*a_const*beam.laser_energy*pow(beam.beam_energy/pchep::electron_mass_GeV, 2);       // eqn.16{max recoil photon energy} (GeV)                                                                                                  
  double kDummy = kprimemax;

  double r_dipole_1;
  double r_dipole_2;
  double th_dipole_1;
  double th_dipole_2;
  double h1 = 0;
  double h2 = 0;

  double thetabend_1 = asin(pchep::c_light_lol*mag[0].dipole*mag[0].length/beam.beam_energy);       // 10.131*pi/180 ! bend angle in Compton chicane (radians)
  double thetabend_2 = asin(pchep::c_light_lol*mag[1].dipole*mag[1].length/beam.beam_energy);       // Dont ask why he original author made light in that units lol

  double R_bend_1    = beam.beam_energy/(pchep::c_light_lol*mag[0].dipole);
  double R_bend_2    = beam.beam_energy/(pchep::c_light_lol*mag[1].dipole);

  QEDasym.open("data/QEDasymP.txt");
  ofstream test;
  test.open("test.dat");

  if(!QEDasym.is_open()){
    std::cout << red << "Cant open file: QEDasymP.txt" << white << std::endl;
    exit(1);
  }
  for (Int_t i = 0; i < nPoints; i++) {
    rho[i]   = (double)i/(nPoints-1);
    kDummy   = rho[i]*kprimemax;

    // Compute height change due to first dipole                                                                                                                                                   

    r_dipole_1 = (beam.beam_energy - kDummy + beam.laser_energy)/(pchep::c_light_lol*mag[0].dipole);
    th_dipole_1 = asin(pchep::c_light_lol*mag[0].dipole*mag[0].length/(beam.beam_energy+beam.laser_energy-kDummy));

    h1 = r_dipole_1*(1-cos(th_dipole_1)) - R_bend_1*(1-cos(thetabend_1)) + (geo.drift_first + mag[1].length)*(tan(th_dipole_1) - tan(thetabend_1));   // x-component for first dipole                            

    // Compute height change due to second dipole                                                                                                                                                  

    r_dipole_2  = (beam.beam_energy - kDummy + beam.laser_energy)/(pchep::c_light_lol*mag[1].dipole);
    th_dipole_2 = asin(pchep::c_light_lol*mag[1].dipole*mag[1].length/(beam.beam_energy+beam.laser_energy-kDummy));

    h2 = r_dipole_2*(1-cos(th_dipole_2)) - R_bend_2*(1-cos(thetabend_2)) + (geo.drift_second)*(tan(th_dipole_1 + th_dipole_2) - tan(thetabend_1 + thetabend_2));   // x-component for second dipole      

    xPrime[i] = (h2 + h1)*cos(thetabend_2)/cos(det_angle-thetabend_2);

    dsdx_0[i] = ((1.0-rho[i]*(1.0+a_const))/(1.0-rho[i]*(1.0-a_const)));
    dsdx[i]   = 2*boost::math::constants::pi<double>()*pow(pchep::classic_e_radius, 2)/100.0*a_const*(rho[i]*rho[i]*(1-a_const)*(1-a_const)/(1-rho[i]*(1.0-a_const))+1.0+(dsdx_0[i] *dsdx_0[i]));
    asym[i]   = 2*boost::math::constants::pi<double>()*pow(pchep::classic_e_radius, 2)/100.0*a_const/dsdx[i]*(1-rho[i]*(1+a_const))*(1.0-1.0/(pow((1.0-rho[i]*(1.0-a_const)),2))) ;
    test << rho[i] << " " << asym[i] << std::endl; 

    if(QEDasym.is_open()) {
      QEDasym << xPrime[i] << "\t" << rho[i] << "\t" << asym[i] << "\t" << dsdx[i] << std::endl;
    } else std::cout << "**Alert: couldn't open file to write QEDasym**\n" << std::endl;
  }

  QEDasym.close();
  test.close();
  xCedge = xPrime[nPoints-1]; // 'xCedge' is used in determining QED asym, hence this should be evaluated before calling the function to fit theoretical asym                                      

  if(debug) {
  }

  TGraph *grtheory = new TGraph("data/QEDasymP.txt", "%lg %lg");

  grtheory->GetXaxis()->SetTitle("dist from compton scattered electrons(m)");
  grtheory->GetYaxis()->SetTitle("#rho");
  grtheory->GetYaxis()->CenterTitle();
  grtheory->SetTitle("#rho to x");
  grtheory->SetMarkerStyle(20);
  grtheory->SetLineColor(2);
  grtheory->SetMarkerColor(2);

  TF1 *fn0 = new TF1("fn0","pol3");
  grtheory->Fit("fn0","0");            //,"0","goff");                                                                                                                                             
  fn0->GetParameters(param);
  ofstream checkfile;

  checkfile.open("data/scheckfileP.txt");
  if(checkfile.is_open()) {
    checkfile << param[0] << "\t" << param[1] << "\t" << param[2] << "\t" << param[3] << std::endl;
  } else std::cout << "**Alert: couldn't open file to write QED fit parameters**\n" << std::endl;
  checkfile.close();
  return xCedge;
}

void ComptonConfig::GetOptions(char **options){

  int i = 0;
  
  std::string flag;

 while(options[i] != NULL){
   flag = options[i];

   if(flag.compare("--config") == 0){
     std::string opt(options[i+1]);
     fConfigDirectory = opt;
     flag.clear();
     std::cout << blue << "Loading exteranl configuration files from:\t" 
	       << fConfigDirectory 
	       << white << std::endl;
   }
    if(flag.compare("--graphics") == 0){
      flag.clear();
      fGraphicsShow = true;
    }
    if(flag.compare("--residuals") == 0){
      flag.clear();
      fResiduals = true;
    }
    if(flag.compare("--fit_range") == 0){
      flag.clear();
      std::string opt(options[i+1]);
      int index = opt.find_first_of(":");
      fLowerLimit = atoi(opt.substr(0, index).c_str());
      fUpperLimit = atoi(opt.substr(index + 1, opt.length()-index).c_str());
    }
    if(flag.compare("--pol_range") == 0){
      flag.clear();
      std::string opt(options[i+1]);
      int index = opt.find_first_of(":");
      fLowerPolLimit = atof(opt.substr(0, index).c_str());
      fUpperPolLimit = atof(opt.substr(index + 1, opt.length()-index).c_str());
    }
    if(flag.compare("--compton_range") == 0){
      flag.clear();
      std::string opt(options[i+1]);
      int index = opt.find_first_of(":");
      fLowerCELimit = atof(opt.substr(0, index).c_str());
      fUpperCELimit = atof(opt.substr(index + 1, opt.length()-index).c_str());
    }

   if(flag.compare("--help") == 0){
     printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
     printf("Usage: ./comptonfit <options>\n");
     printf("         --graphics \tGraphical output.\n");
     printf("         --residuals \tCalculate and plot residuals of asymmetry fit.\n");
     printf("         --fit_range \tDefine the range of strips over which the fit will be done.\n");
     printf("                     \tex. --fit_range 1:120 \n");
     printf("         --pol_range \tDefine the range of polarization over which the fit will be done.\n");
     printf("                     \tex. --pol_range 0.8:1.0 \n");
     printf("         --compton_range \tDefine the range of compton edge over which the fit will be done.\n");
     printf("                     \tex. --compton_range 115:120 \n");
     printf("\n~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n");
     exit(0);
   }
   i++;
 }
}

void ComptonConfig::GetConfigurationFiles(const boost::filesystem::path& root, const std::string ext)
{

  if(!boost::filesystem::exists(root) || !boost::filesystem::is_directory(root)) return;

  boost::filesystem::recursive_directory_iterator it(root);
  boost::filesystem::recursive_directory_iterator endit;

  while(it != endit)
    {
      if(boost::filesystem::is_regular_file(*it) && it->path().extension() == ext) file_list.push_back(it->path().filename());
      ++it;

    }
}

void ComptonConfig::ReadConfiguration()
{

  std::fstream config;

  for(int i = 0; i < (int)file_list.size(); i++){
    std::cout << blue << "Opening file:\t " << "\t" << file_list[i] << white << std::endl;

    config.open(Form("config/%s", file_list[i].string().c_str()), std::ios_base::in);
    if(!config.is_open()){
      std::cout << red << "Error opening detector config file" << white << std::endl;
      exit(1);
    }
    // Determine the type of config file we are dealing with. Search for 'keyword' in string.

    if(file_list[i].string().find("magnet")!= std::string::npos){
      std::cout << green << "Found magnet file." << white << std::endl;
      ParseMagnetFile(config);
    }
    if(file_list[i].string().find("detector")!= std::string::npos){
      std::cout << green << "Found detector file." << white << std::endl;
      ParseDetectorFile(config);
    }
    if(file_list[i].string().find("parameter")!= std::string::npos){
      std::cout << green << "Found beam paramters file." << white << std::endl;
      ParseBeamParametersFile(config);
    }
    else {
      std::cout << red << "Unknow file type: skipping." << white << std::endl;
      continue;
    }
  }
  std::cout << blue << "Finished reading configuration files." << white << std::endl;
  SetMagnetOrder();

  return;
}

void ComptonConfig::ParseMagnetFile(std::fstream &config)
{

  std::string line;

  char *token;
				    
  mag.push_back(magnet());       // Push a new magnet definition and move on

  while(config.good()){
    getline(config, line);
    token = new char[line.size() + 1];
    strcpy(token, line.c_str());       // convert line to c string and copy to char array token
    token = strtok(token, " ,.");      // breaks token into tokens defined by delimiters

    while(token){
      if(strcmp("magnet", token) == 0){
	// mag.push_back(magnet());       // Push a new magnet definition and move on
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
      }
      if(strcmp("theta", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	mag.back().theta = atof(token);
      }
      if(strcmp("x", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	mag.back().x = atof(token);
      }
      if(strcmp("y", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	mag.back().y = atof(token);
      }
      if(strcmp("z", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	mag.back().z = atof(token);
      }
      if(strcmp("length", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	mag.back().length = atof(token);
      }
      if(strcmp("dipole", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	mag.back().dipole = atof(token);
      }
      else
	token = strtok(NULL, " ,");
    }
  }
  config.close();
  return;
}

void ComptonConfig::ParseDetectorFile(std::fstream &config)
{

  std::string line;

  char *token;
				    
   while(config.good()){
    getline(config, line);
    token = new char[line.size() + 1];
    strcpy(token, line.c_str());       // convert line to c string and copy to char array token
    token = strtok(token, " ,");      // breaks token into tokens defined by delimiters

    while(token){
       if(strcmp("theta", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	det.theta = atof(token);
      }
      if(strcmp("x", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	det.x = atof(token);
      }
      if(strcmp("y", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	det.y = atof(token);
      }
      if(strcmp("z", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	det.z = atof(token);
      }
      if(strcmp("width", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	det.width = atof(token);
      }
      if(strcmp("spacing", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	det.spacing = atof(token);
      }
      else
	token = strtok(NULL, " ,");
    }
  }

  config.close();
  return;
}

void ComptonConfig::ParseBeamParametersFile(std::fstream &config)
{

  std::string line;

  char *token;
				    
   while(config.good()){
    getline(config, line);
    token = new char[line.size() + 1];
    strcpy(token, line.c_str());       // convert line to c string and copy to char array token
    token = strtok(token, " ,");      // breaks token into tokens defined by delimiters

    while(token){
       if(strcmp("beam", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	beam.beam_energy = atof(token);
      }
      if(strcmp("laser", token) == 0){
	token = strtok(NULL, " ,");   // Here the extra strtok(NULL, " .,") keeps scanning for next token
	beam.laser_energy = atof(token);
      }
      else
	token = strtok(NULL, " ,");
    }
  }
  config.close();
  return;
}

void ComptonConfig::SetMagnetOrder()
{
  std::cout << "Setting dipole order correctly." << std::endl;

  if(mag[0].z > mag[1].z) return;
  else{

    magnet temp;
    temp.x = mag[0].x;
    mag[0].x = mag[1].x;
    mag[1].x = temp.x;
    
    temp.y   = mag[0].y;
    mag[0].y = mag[1].y;
    mag[1].y = temp.y;
    
    temp.z   = mag[0].z;
    mag[0].z = mag[1].z;
    mag[1].z = temp.z;

    temp.length   = mag[0].length;
    mag[0].length = mag[1].length;
    mag[1].length = temp.length;

    temp.dipole   = mag[0].dipole;
    mag[0].dipole = mag[1].dipole;
    mag[1].dipole = temp.dipole;

    temp.theta   = mag[0].theta;
    mag[0].theta = mag[1].theta;
    mag[1].theta = temp.theta;
  }
  return;
}

void ComptonConfig::InitGraphicsEngine(int Argc, char **Argv)
{
  std::cout << green << "<<<< Initialize Graphics Engine." << white << std::endl;
  app = new TApplication("App", &Argc, Argv);

}

void ComptonConfig::RunGraphicsEngine()
{
  std::cout << green << "<<<< Running Graphics Engine." << white << std::endl;
  app->Run();
}
#endif

