/* @doc curveFitting Math SpaceAndTime */

/****************************************************************************
*                              arg-cpp-library
*           Atmospheric Research Group, University of Saskatchewan
****************************************************************************/
/*<--------------------------------------------------------------------------
 *'						class nxVector2D
 *	Implements 2-D vectors as X,Y
 *>------------------------------------------------------------------------*/

class nxVector2D
{
	private:
		double		x;
		double		y;

	public:
						nxVector2D	()								  { x = 0,  y = 0;}
						nxVector2D	( double ax, double ay)			  { x = ax; y = ay;}
		double			X			() const						  { return x;}
		double			Y			() const						  { return y;}
		void			SetX		( double ax )					  { x = ax;}
		void			SetY		( double ay )					  { y = ay;}
		void			SetTo		( double ax, double ay )		  { x = ax; y = ay;}
		nxVector2D		Abs			()							const { return nxVector2D( fabs(x), fabs(y));}
		double			Magnitude	()							const { return sqrt(x*x+y*y);}
		double			DotProduct	( const nxVector2D& other ) const { return x*other.x + y*other.y;}
		double			DistanceTo	( const nxVector2D& p2)			  { return (*this-p2).Magnitude();}
		nxVector2D		operator+	( const nxVector2D& other ) const { return nxVector2D( x + other.x,  y + other.y); }
		nxVector2D		operator+	( double constant)			const { return nxVector2D( x + constant, y + constant); }
		nxVector2D		operator-	( const nxVector2D& other ) const { return nxVector2D( x - other.x,  y - other.y); }
		nxVector2D		operator-	( double constant)			const { return nxVector2D( x - constant, y - constant); }
		nxVector2D		operator/	( double denom )			const { return nxVector2D( x/denom,      y/denom);}
		nxVector2D		operator/	( const nxVector2D& other )	const { return nxVector2D( x/other.x,    y/other.y);}
		nxVector2D		operator*	( double mult     )			const { return nxVector2D( x*mult,       y*mult);}
		nxVector2D		operator*	( const nxVector2D& other )	const { return nxVector2D( x*other.x,    y*other.y);}
		
		void			Normalize()
		{
			double m = Magnitude();
			if (m > 0.0)
			{
				x = x/m;
				y = y/m;
			}
		}

		nxVector2D		NormalizedXY() const
		{
			double m = Magnitude();
			double ax = 0;
			double ay = 0;
			if (m > 0.0)
			{
				ax = x/m;
				ay = y/m;
			}
			return nxVector2D(ax,ay);
		}
		nxBOOL			OnlyOneComponentIsZero()	{ return ((x == 0.0) && (y != 0.0)) || ((y == 0.0 && x != 0.0));}
       	nxBOOL			OnlyOneCompononentLessThan( double thres) { return ((x < thres) && (y >= thres))  || (( x>= thres) && (y < thres));}
		nxBOOL			BothComponentsGreaterThan ( double thres) { return ((x >= thres)  && (y >= thres));}
};

nxVector2D	operator*( double mult,     const nxVector2D& other ); //{ return other*mult;}
nxVector2D	operator+( double constant, const nxVector2D& other ); //{ return other+constant; }
nxVector2D	operator-( double constant, const nxVector2D& other ); //{ return other-constant; }

