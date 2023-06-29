#include <vector>
#include <algorithm>
using namespace std;

//Some standard mathematical functions
vector<double>				LinearInterpolate	( const vector<double>& x,			const vector<double>& y, const vector<double>& xx );
vector<double>				spline				( const vector<double>& x,			const vector<double>& y, const vector<double>& xx );
vector<double>				PolynomialRegression( const vector<double>& x,			const vector<double>& y, int n );
double						Median				( const vector<double>& x );

//Input/output vectors to/from text files
bool						LoadVectorFromFile	( string filename, vector<vector< double > >* v, size_t num_columns );
void						OutputVectorToText	( string filename, vector< double > &v );

//Crude (and slow) but convenient implementation of some mathematical functions with vectors
vector<vector<double> >		Transpose			( const vector<vector<double> >	&x );
vector<vector<double> >		Invert				( const vector<vector<double> >	&x );
vector<vector<double> >		Diagonal			( const vector<vector<double> >	&x );
vector<vector<double> >		Multiply			( const vector<vector<double> >	&x, const vector<vector<double> > &y );
vector<vector<double> >		Multiply			( const vector<vector<double> >	&x, double						  y );
vector<vector<double> >		Multiply			( double						 y, const vector<vector<double> > &x );
vector<double>				Multiply			( double						 y, const vector<double>		 &x );
vector<double>				Multiply			( const vector<vector<double> >	&x, const vector<double>		 &y );
double						Multiply			( const vector<double>			&x, const vector<double>		 &y );
vector<vector<double> >		Add 				( const vector<vector<double> >	&x, const vector<vector<double> > &y );
vector<double>				Add 				( const vector<double>			&x,	const vector<double>		 &y );