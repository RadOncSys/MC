// Radiation Oncology Monte Carlo open source project
//
// Author: [2005-2017] Gennady Gorlachev (ggorlachev@roiss.ru) 
//---------------------------------------------------------------------------
#pragma once

#include <string>
#include <vector>

short GetTwoStringsFromLine(const std::string&, std::string&, std::string&, const char* pmask = nullptr);
int GetFloatArray(const std::string& line, double* x, int n = 0);
int GetFloatArray(const std::string&, std::vector<double>& a);
int GetIntArray(const std::string&, int*, int n = 0);
int GetIntArray(const std::string& line, std::vector<int>& a);
