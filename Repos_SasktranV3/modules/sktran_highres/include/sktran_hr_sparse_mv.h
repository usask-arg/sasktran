//#include "sktran_hr_internals.h"

void SKTRAN_HR_Sparse_mmv( std::vector<double>& mvals,
						   std::vector<int>&    mcolind,
						   std::vector<int>&    mrowptr,
						   std::vector<double>& v,
						   std::vector<double>& result,
						   int numthreads);
