

class ISKMie_Stub_Wiscombe : public ISKMie_Stub {
private:
	sk_MieSphericalWiscombeWrapper m_mie;

	void MakeVectorSetFunctions();
	void MakeScalarSetFunctions();

public:
	ISKMie_Stub_Wiscombe();
	virtual ~ISKMie_Stub_Wiscombe() {};

	virtual bool Calculate(double lambda, double radius, double refrac_real, double refrac_imag) override final;

	virtual double								Qext() override final;
	virtual double								Qsca() override final;
	virtual double								Qabs() override final;
	virtual double								Cext() override final;
	virtual double								Csca() override final;
	virtual double								Cabs() override final;
	virtual nx1dArray< std::complex<double> >* S1() override final;
	virtual nx1dArray< std::complex<double> >* S2() override final;
	virtual nx2dArray<double>* PMom() override final;

	virtual std::complex<double>				SForward() override final;
	virtual std::complex<double>				SBackward() override final;
	virtual std::complex<double>				TForward(int i) override final;
	virtual std::complex<double>				TBackward(int i) override final;
	virtual double								Spike() override final;
};