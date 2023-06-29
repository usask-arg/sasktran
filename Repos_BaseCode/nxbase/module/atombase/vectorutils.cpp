#include "nxbase_linearalgebra.h"
//#include "nxbase_math.h"
//#include "../math/linearalgebra/nxlapack.cpp"
//#include "vectorutils.h"
using namespace std;

std::vector<double> LinearInterpolate( const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xx )
{
	size_t numElements = xx.size();
	size_t elementCtr;
	std::vector<double> yy( numElements );
	double slope;

	
	for (elementCtr = 0; elementCtr < numElements; elementCtr++)
	{
		double xVal = xx[elementCtr];
		double yVal;
		if (xVal <= x[0])
		yVal = y[0];
		else	if (xVal >= x[x.size()-1])
				yVal = y[y.size()-1];
				else
				{
					int ctr=0;
					while (xVal > x[ctr])
					ctr++;
					
					slope = (y[ctr]-y[ctr-1])/(x[ctr]-x[ctr-1]);
					
					double offset = xVal-x[ctr-1];
					
					yVal = y[ctr-1]+slope*offset;
				}
		yy[elementCtr] = yVal;				
				
	}


	return(yy);
}

std::vector<double> spline( const std::vector<double>& x, const std::vector<double>& y, const std::vector<double>& xx )
{
	int numCubics = int(x.size()-1);

	/* Step 1: Calculate the h's */
	double* h = new double[numCubics];
	for( int i=0; i<numCubics; i++ )
		h[i] = x[i+1]-x[i];

	/* Steps 2+3: Calculate the d's and u's*/
	double* d = new double[numCubics];
	double* u = new double[numCubics-1];

	d[0] = (y[1]-y[0])/h[0];
	for( int i=1; i<numCubics; i++ )
	{
		d[i] = (y[i+1]-y[i])/h[i];
		u[i-1] = 6*(d[i]-d[i-1]);
	}

	/* Step 4: Construct System of Equations, Matrix A and Solutions m */
	// A is tri-diagonal. This fact will be exploited to solve the system.
	double* A = new double[numCubics-1];
	double* m = new double[numCubics-1];

	// First pass eliminates lower diagonal
	A[0] = 1.5*h[0]+2.0*h[1];
	m[0] = u[0]-3*d[0];
	for( int i=1; i<numCubics-2; i++ )
	{
		A[i] = 2*(h[i]+h[i+1])-h[i]*h[i]/A[i-1];
		m[i] = u[i]-m[i-1]*h[i]/A[i-1];
	}
	A[numCubics-2] = 2*h[numCubics-2]+1.5*h[numCubics-1]-h[numCubics-2]*h[numCubics-2]/A[numCubics-3];
	m[numCubics-2] = u[numCubics-2]+3*d[numCubics-1]-m[numCubics-3]*h[numCubics-2]/A[numCubics-3];

	// Second pass eliminates upper diagonal
	m[numCubics-2] /= A[numCubics-2];
	for( int i=numCubics-3; i>=0; i-- )
	{
		m[i] = ( m[i] - m[i+1]*h[i+1] ) / A[i];
	}

	/* Step 5: Calculate polynomial coefficients. */
	// a*x^3 + b*x^2 + c*x + d
	// d's are just the y's and need not be recalculated.
	std::vector<double> a(numCubics);
	std::vector<double> b(numCubics);
	std::vector<double> c(numCubics);

	double m_first = 3.0*d[0]/h[0] - m[0]/2.0;
	a[0] = (m[0]-m_first)/(6*h[0]);
	b[0] = m_first/2.0;
	c[0] = d[0] - h[0]*(2*m_first+m[0])/6.0;
	for( int i=1; i<numCubics-1; i++ )
	{
		a[i] = (m[i]-m[i-1])/(6*h[i]);
		b[i] = m[i-1]/2.0;
		c[i] = d[i] - h[i]*(2*m[i-1]+m[i])/6.0;
	}
	double m_last = -3.0*d[numCubics-1]/h[numCubics-1] - m[numCubics-2]/2.0;
	a[numCubics-1] = (m_last-m[numCubics-2])/(6*h[numCubics-1]);
	b[numCubics-1] = m[numCubics-2]/2.0;
	c[numCubics-1] = d[numCubics-1] - h[numCubics-1]*(2*m[numCubics-2]+m_last)/6.0;

	delete[] h;
	delete[] d;
	delete[] u;
	delete[] A;
	delete[] m;

	/* Step through all xx's and find yy's */
	std::vector<double> yy( xx.size() );
	std::vector<double>::const_iterator itr_x, itr_y, itr_xx, itr_xx_lower, itr_xx_upper;
	std::vector<double>::iterator itr_yy, itr_yy_lower, itr_yy_upper;
	std::vector<double>::const_iterator itr_a, itr_b, itr_c;
	std::vector<double>::difference_type diff_yx = y.begin() - x.begin();
	std::vector<double>::difference_type diff_ax = a.begin() - x.begin();
	std::vector<double>::difference_type diff_bx = b.begin() - x.begin();
	std::vector<double>::difference_type diff_cx = c.begin() - x.begin();

	itr_xx = itr_xx_lower = std::lower_bound( xx.begin(), xx.end(), x.at(0) );
	if( itr_xx==xx.begin() )
	{
		itr_x = std::lower_bound( x.begin(), x.end(), xx.at(0) );
	}
	else
		itr_x = x.begin();

	while( itr_xx<xx.end() )
	{
		if( itr_x > x.end()-2 )
			break;
		itr_y	= itr_x + diff_yx;
		itr_a	= itr_x + diff_ax;
		itr_b	= itr_x + diff_bx;
		itr_c	= itr_x + diff_cx;
		itr_yy	= itr_xx - xx.begin() + yy.begin();
		double dx = (itr_xx[0]-itr_x[0]);
		itr_yy[0] = ((itr_a[0]*dx+itr_b[0])*dx+itr_c[0])*dx+itr_y[0];
		itr_xx_upper = ++itr_xx;
		while( (itr_x<x.end()-1) && (itr_xx[0] > itr_x[1]) )
			itr_x++;
	}

	itr_yy_lower = itr_xx_lower - xx.begin() + yy.begin();
	std::fill( yy.begin(), itr_yy_lower, y.front() );
	itr_yy_upper = itr_xx_upper - xx.begin() + yy.begin();
	std::fill( itr_yy_upper, yy.end(), y.back() );

	return(yy);
}

vector<double> PolynomialRegression( vector< double > &x, vector<double> &y, int n )
{
	// HAS NOT BEEN DEBUGGED (like at all)
	// [Inverse(Xt*X)]*[X*y]

	LapackMatrix		Xl;
	LapackMatrix		XTl;
	LapackMatrix		yl;
	LapackMatrix		al;
	LapackFactorLU		solver;

	size_t						i,j;
	vector< double >			a;
	vector<vector< double > >	X;

	//create X matrix
	X.resize( x.size() );
	for ( j=0; j<x.size(); j++ )
	{
		X[j].resize( n+1 );
		for( i=0; i< (size_t)n+1; i++ )
			X[j][i] = pow( x[j], i );
	}

	Xl.CreateFromVectorVector(X);
	XTl.CreateFromVectorVector( Transpose(X) );
	yl.CreateFromColumnVector( y );

	solver.Solve(  XTl*Xl, XTl*yl, &al );		//solve coefficients

	for( i=0; i< (size_t)n+1; i++ )
		a.push_back(al(i+1,1));

	return a;
}

double Median( const vector< double > &x )
{
	double				median;
	vector< double >	t;
	t = x;
	size_t size = t.size();

	sort( t.begin(), t.end() );

	if (size % 2 == 0 )
		median = ( t[size/2-1]+t[size/2] )/2;
	else
		median = t[size/2];

	return median;

}

bool LoadVectorFromFile( string filename, vector<vector< double > >* v, size_t num_columns )
{
	bool	ok=true;
	size_t	i;
	double	temp;
	//Instantiate read file object
	std::ifstream inputFile(filename.c_str());

	//Check if the file is open
	if(inputFile.fail() == 1)
		ok = false;
	else
	{
		v->clear();
		v->resize(num_columns);
		while(!inputFile.eof()) 
		{
			//Read file into class
			for( i=0; i<num_columns; i++ )
			{
				inputFile >> temp;
				v->at(i).push_back(temp);
			}
		}
		//Close file
		inputFile.close();
	}
	return ok;
}

void OutputVectorToText( string filename, vector< double >& v )
{
	size_t j;
	ofstream output(filename.c_str(), std::ofstream::app);
	for(j=0;j<v.size();j++)
		output<<v[j]<<" ";

	output<<endl;
	output.close();
}

vector<vector<double> > Transpose( const vector<vector<double> > &x )
{
	size_t i,j;
	vector<vector<double> > xt;

	xt.resize(x[0].size() );
	for (i=0; i<x[0].size(); i++)
		xt[i].resize( x.size(),0.0);
	
	for(i=0; i<x.size(); i++)
		for(j=0; j<x[i].size(); j++)
			xt[j][i] = x[i][j];

	return xt;
}

vector<vector<double> > Diagonal( const vector<vector<double> > &x )
{
	size_t					i;
	vector<vector<double> >	d;

	d.resize( x.size() );
	for(i=0; i<x.size(); i++)
	{
		d[i].resize( x[i].size(), 0.0 );
		d[i][i] = x[i][i];
	}

	return d;
}

vector<vector<double> > Invert( const vector<vector<double> > &x )
{
	size_t							i,j;
	LapackMatrix				X;
	LapackMatrix				Xi;
	vector<vector< double > >	xi;

	if ( x.size()==x[0].size() )	//check that matrix is square
	{
		xi.resize( x.size() );
		for ( i=0; i<x.size(); i++)
			xi[i].resize(x.size() );

		X.CreateFromVectorVector( x );
		Xi = X.LUInverse();
	}

	for ( i=0; i<x.size(); i++ )
		for ( j=0; j<x.size(); j++ )
			xi[i][j] = Xi.At( i+1, j+1 );

	return xi;
}

vector<vector<double> > Multiply( const vector<vector<double> > &x, const vector<vector<double> > &y )
{
	size_t							row, col, inner;
	vector< vector< double > >	product;

	if ( x.size()==y[0].size()  &&  x[0].size() == y.size() )	//check that matrices can be multiplied
	{
		product.resize(x.size());
		for ( row=0; row < x.size(); row++ )
			product[row].resize( y[0].size(), 0.0 );

		for ( row = 0; row < x.size(); row++ )
		{
			for ( col = 0; col < y[0].size(); col++ ) 
			{
				// Multiply the row of A by the column of B to get the row, column of product.
				for ( inner = 0; inner < x[0].size(); inner++ ) 
				{
					product[row][col] += x[row][inner] * y[inner][col];
				}
			}
		}
	}
	return product;
}
vector<double>	Multiply( const vector<vector<double> > &x, const vector<double>	&y )
{
	size_t							row, inner;
	vector< double >			product;

	if ( x[0].size() == y.size() )	//check that matrices can be multiplied
	{
		product.resize(x.size(), 0.0);
		for ( row = 0; row < x.size(); row++ )
		{
			// Multiply the row of A by the column of B to get the row, column of product.
			for ( inner = 0; inner < x[0].size(); inner++ ) 
			{
				product[row] += x[row][inner] * y[inner];
			}
		}
	}
	return product;
}
double	Multiply( const vector<double>	&x, const vector<double> &y )
{
	size_t							inner;
	double						product=0;

	if ( x.size() == y.size() )	//check that matrices can be multiplied
	{
		// Multiply the row of A by the column of B to get the product.
		for ( inner = 0; inner < x.size(); inner++ ) 
		{
			product += x[inner] * y[inner];
		}
	}
	return product;
}
vector<vector<double> >	Multiply( const vector<vector<double> >	&x, double  y )
{
	size_t							row, col;
	vector< vector< double > >	product;

		product.resize(x.size());
		for ( row=0; row < x.size(); row++ )
			product[row].resize( x[0].size(), 0.0 );

		for ( row = 0; row < x.size(); row++ )
		{
			for ( col = 0; col < x[0].size(); col++ ) 
			{
				product[row][col] = x[row][col] * y;
			}
		}

	return product;
}
vector<vector<double> >	Multiply( double y, const vector<vector<double> >	&x )
{
	size_t							row, col;
	vector< vector< double > >	product;

		product.resize(x.size());
		for ( row=0; row < x.size(); row++ )
			product[row].resize( x[0].size(), 0.0 );

		for ( row = 0; row < x.size(); row++ )
		{
			for ( col = 0; col < x[0].size(); col++ ) 
			{
				product[row][col] = x[row][col] * y;
			}
		}

	return product;
}
vector<double>	Multiply( double y, const vector<double>	&x )
{
	size_t							row;
	vector<  double >	product;

	product.resize(x.size(), 0.0);

	for ( row = 0; row < x.size(); row++ )
		product[row] = x[row] * y;
	
	return product;
}
vector<vector<double> > Add ( const vector<vector<double> >	&x, const vector<vector<double> > &y )
{
	size_t							row, col;
	vector< vector< double > >	result;

	if ( x.size() == y.size() && x[0].size() == y[0].size() )
	{
		result.resize(x.size());
		for ( row=0; row < x.size(); row++ )
			result[row].resize( x[0].size(), 0.0 );

		for ( row = 0; row < x.size(); row++ )
		{
			for ( col = 0; col < x[0].size(); col++ ) 
			{
				result[row][col] = x[row][col] + y[row][col];
			}
		}
	}
	return result;
}
vector<double>	Add( const vector<double>	&x,	const vector<double> &y )
{
	size_t					row;
	vector<  double >	result;

	if ( x.size() == y.size() )
	{
		result.resize(x.size(), 0.0);
		for ( row = 0; row < x.size(); row++ )
			result[row] = x[row] + y[row];
	}
	return result;
}