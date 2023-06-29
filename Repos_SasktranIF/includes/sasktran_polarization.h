

/*-----------------------------------------------------------------------------
 *					ISKBasisDirection		 2015- 9- 28*/
 /**	The three directions that specify the basis for the ray and its polarization
  *	Directions are specified in ECEF coordinates. The basis must be constructed
  *	such that m_theta cross m_phi is equal to m_propagation
 **/
 /*---------------------------------------------------------------------------*/

class ISKBasisDirection
{
	private:
		nxVector			m_propagation;		// 0 The ray propagation direction
		nxVector			m_theta;			// 1 The theta direction
		nxVector			m_phi;				// 2 The phi direction

	public:
		void				Assign		(const nxVector& prop, const nxVector& theta, const nxVector& phi);
		const nxVector&		Propagation	() const { return m_propagation; }		// 0 The ray propagation direction
		const nxVector&		Theta		() const { return m_theta; }				// 1 The theta direction
		const nxVector&		Phi			() const { return m_phi; }				// 2 The phi direction
		ISKBasisDirection&	operator=   (const ISKBasisDirection& other) { this->Assign(other.Propagation(), other.Theta(), other.Phi()); return *this; }
};


/*-----------------------------------------------------------------------------
 *					IQUV		 2015- 9- 29*/
 /** **/
 /*---------------------------------------------------------------------------*/

struct IQUV
{
	double	I;
	double	Q;
	double	U;
	double	V;
};

/*-----------------------------------------------------------------------------
 *					ISKStokesVector		 2015- 9- 29*/
 /** **/
 /*---------------------------------------------------------------------------*/

class ISKStokesVector
{
	private:
		IQUV							m_stokes;				// The stokes vector [I, Q, U, V]
		ISKBasisDirection				m_basis;				// Coordinate basis the stokes vector is defined in.


	public:
										ISKStokesVector			();
										ISKStokesVector			(const IQUV& stokes, const ISKBasisDirection& new_basis);
		void							Assign					(const IQUV& stokes, const ISKBasisDirection& new_basis);
		const IQUV&						Stokes					() const { return m_stokes; }
		const ISKBasisDirection&		Basis					() const { return m_basis; }
		void							to_new_basis			(const ISKBasisDirection& new_basis);
		void							to_new_basis			(const nxVector& prop, const nxVector& theta, const nxVector& phi);
		const nxVector&					propagation_direction	() const { return m_basis.Propagation(); }			//!< Returns the propagation direction of the current basis in ECEF coordinates
		const nxVector&					theta_direction			() const { return m_basis.Theta(); }					//!< Returns the theta direction of the current basis in ECEF coordinates
		const nxVector&					phi_direction			() const { return m_basis.Phi(); }					//!< Returns the phi direction of the current basis in ECEF coordinates
		double							I						() const { return m_stokes.I; }
		double							Q						() const { return m_stokes.Q; }
		double							U						() const { return m_stokes.U; }
		double							V						() const { return m_stokes.V; }
};


