#pragma once

template<typename T>
class SKTRAN_HR_WF_TVector
{
    private:
		std::vector<double> m_loweraltitudes;
		std::vector<double> m_upperaltitudes;
		bool				m_issorted;

	public:
		std::vector<T> m_storage;

		void push_back(T item) { 
			m_storage.push_back(item); 
			m_loweraltitudes.push_back(item.PerturbationAltitudeLower());
			m_upperaltitudes.push_back(item.PerturbationAltitudeUpper()); 
			m_issorted = std::is_sorted(std::begin(m_loweraltitudes), std::end(m_loweraltitudes)) && std::is_sorted(std::begin(m_upperaltitudes), std::end(m_upperaltitudes));
		}

		std::pair<size_t, size_t> bounding_indicies(double altitude) const;

};

class SKTRAN_HR_WF_Store
{
	private:
		SKTRAN_HR_WF_TVector<SKTRAN_HR_Perturbation_Absorption_Box> m_boxstore;

		std::vector<double>	m_boxstorealts;
		std::vector<double> m_boxstorewidths;

		SKTRAN_HR_WF_TVector<SKTRAN_HR_Perturbation_Absorption_Linear> m_linearstore;


	public:
		void add_item(SKTRAN_HR_Perturbation_Absorption_Box item)
		{
			m_boxstore.push_back(item);
			m_boxstorealts.push_back(item.PerturbationCenterAltitude());
			m_boxstorewidths.push_back(item.PerturbationAltitudeWidth() / 2);
		}

		void add_item(SKTRAN_HR_Perturbation_Absorption_Linear item)
		{
			m_linearstore.push_back(item);
		}

		size_t StoreSize() const { return m_boxstore.m_storage.size() + m_linearstore.m_storage.size(); }

		void ExtinctionPerturbation(const HELIODETIC_POINT& location, std::vector<double>& value) const;
		void AddExtinctionPerturbation(const HELIODETIC_POINT & location, std::vector<double>& value, double ds) const;
		void PerturbationLocation(const SKTRAN_CoordinateTransform_V2& coords, std::vector<HELIODETIC_POINT>& data) const;
		void PerturbationAltitudeWidth(std::vector<double>& data) const;

		void AddGeometryToRayTracer(std::shared_ptr<const SKTRAN_CoordinateTransform_V2>& coords, SKTRAN_RayTracer_Straight_Generic& raytracer) const;

		void SetPertVal(double value);

		const SKTRAN_HR_Perturbation_Base* RawAccess(size_t idx) const;
		SKTRAN_HR_Perturbation_Base* RawAccess(size_t idx);
};