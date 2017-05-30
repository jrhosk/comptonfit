#ifndef compon_config_hh
#define compton_config_hh

#include <vector>
#include <string>

#include "TApplication.h"

#include "boost/filesystem.hpp"

class ComptonConfig {

private:

  bool kRacCor;

  std::string fConfigDirectory;

  std::vector <boost::filesystem::path> file_list;

public:

  // Structure

  struct detector {
    double width;             // width of each strip of electron detector
    double spacing;           // distance between each strip
    double ce_to_det_bot   ;  // distance from compton edge to detector bottom
    double theta;             // detector tilt in theta
    double position;          // detector position wrt last dipole

    double x;
    double y;
    double z;
 
  };

  struct geometry {
    double drift_first;   // first drift distance in the six magnoet chicane setup
    double drift_second;  // second drift distance in the six magnoet chicane setup

  };

  struct magnet {
    double x;             // x position of magnet 
    double y;             // y position of magnet 
    double z;             // z position of magnet 
    double theta;         // theta of magnet 
    double length;        // length of magnet 
    double dipole;
  };

  struct parameters { 

  double beam_energy;
  double laser_energy;

  };

  double fLowerLimit;
  double fUpperLimit;
  double fLowerPolLimit;
  double fUpperPolLimit;
  double fLowerCELimit;
  double fUpperCELimit;

  double xCedge;
  double kprimemax;
  double param[4];

  detector det;
  geometry geo;
  parameters beam;

  TApplication *app;

  bool fGraphicsShow;
  bool fResiduals;

  std::vector <magnet> mag;

  std::vector <double> sub;
  // Functions

  ComptonConfig();
  ~ComptonConfig();

  double CrossSectionFit(double *thisStrip, double *parCx);
  double TheoreticalAsymFit(double *thisStrip, double *par);
  double RhoToPositionMap(double initCE);
  double ComputeAlpha();

  void GetConfigurationFiles(const boost::filesystem::path& root, const std::string ext);
  void ReadConfiguration();
  void GetOptions(char **options);
  void CalculateDriftDistance();
  void ParseMagnetFile(std::fstream &file);
  void ParseDetectorFile(std::fstream &file);
  void ParseBeamParametersFile(std::fstream &file);
  void SetMagnetOrder();
  void InitGraphicsEngine(int, char** );
  void RunGraphicsEngine();

};
#endif
