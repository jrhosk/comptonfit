#define DataFile_cxx
#include "DataFile.hh"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#ifdef DataFile_cxx

void DataFile::ReadFile()
{

  std::string line;

  char *token;

  while(!eof()){
    std::getline(*this, line);
    token = new char[line.size() + 1];
    strcpy(token, line.c_str());
    token = strtok(token, " ,\t");
    while(token){
      if(!token){
	x.push_back(0); 
	y.push_back(0); 
	y_error.push_back(0); 
	break;
      }
      x.push_back(atof(token));

      token = strtok(NULL, " ,\t");
      if(!token){
	y.push_back(0); 
	y_error.push_back(0);
	break;
      }
      y.push_back(atof(token));

      token = strtok(NULL, " ,\t");
      if(!token){
	y_error.push_back(0); 
	break;
      }
      y_error.push_back(atof(token));

      token = strtok(NULL, " ,\t");
    }
  }
  SetExtrema();

  return;
}

void DataFile::SetExtrema()
{
  fUpperX = x[0];
  fLowerX = x[0];
  fUpperY = y[0];
  fLowerY = y[0];

  for(int i = 0; i < (int)GetArraySize(); i++){
    if(x[i] > fUpperX) fUpperX = x[i];
    if(x[i] < fLowerX) fLowerX = x[i];
    if(y[i] > fUpperY) fUpperY = y[i];
    if(y[i] < fLowerY) fLowerY = y[i];

    fMeanY += fMeanY;
  }

  fMeanY /= GetArraySize();
  return;
}

std::vector <double> DataFile::GetX() const
{
  return x;
}
std::vector <double> DataFile::GetY() const
{
  return y;
}
std::vector <double> DataFile::GetYError() const
{
  return y_error;
}

int DataFile::GetArraySize() const
{
  
  if( ( x.size() != y.size() ) || (x.size() != y_error.size()) ){
    std::cerr << "Array sizes do not match!" << std::endl;
    exit(1);
  }
  else
    return(x.size());
  
}

void DataFile::ScaleEntries(double factor = 1.0)
{
  for(int i = 0; i < (int)GetArraySize(); i++){
    y[i] *= factor;
    y_error[i] *= factor;
  }
  return;
}

void DataFile::ScaleError(double factor = 1.0)
{
  for(int i = 0; i < (int)GetArraySize(); i++){
    y_error[i] *= factor;
  }
  return;
}

#endif
