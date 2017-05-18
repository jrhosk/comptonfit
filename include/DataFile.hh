#ifndef DataFile_h
#define DataFile_h

#include <fstream>
#include <ostream>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <string>
#include <cstring>

class DataFile : private std::fstream {

private:

//   std::ifstream file;

  std::vector <double> x;
  std::vector <double> y;
  std::vector <double> y_error;

public:

  double fUpperX;
  double fUpperY;
  double fLowerX;
  double fLowerY;
  double fMeanY;

  DataFile(const char* filename = 0, ios_base::openmode mode = ios_base::in) : std::fstream(filename, mode)
  {
    if( (good()) ){
      std::cout << "Opening file: " << filename<< std::endl;
    }
    else{
      std::cerr << "Error opening file: " << std::endl;
      exit(1);
    }
  };

  ~DataFile(){close();}

  void ReadFile();
  void SetExtrema();
  void ScaleEntries(double);
  void ScaleError(double);

  int GetArraySize() const;

  std::vector <double> GetX() const;
  std::vector <double> GetY() const;
  std::vector <double> GetYError() const;
};

#endif
