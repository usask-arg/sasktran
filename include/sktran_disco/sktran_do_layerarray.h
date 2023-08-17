#pragma once
#include "sktran_disco/sktran_do.h"
#include <sasktran2/atmosphere/atmosphere.h>

namespace sasktran_disco
{
	template<Propagating dir, int NSTOKES, int CNSTR=-1>
	class OpticalLayerArrayIterator;

	template <int NSTOKES, int CNSTR=-1>
	using OpticalLayerArrayROP = ReadOnlyProperties<
		BasicProperties<NSTOKES>,
		SolarProperties<NSTOKES>,
		UserSpecProperties,
		TestProperties
	>;

	// A vector of OpticalLayers. Provides additional helpful functions.
    template <int NSTOKES, int CNSTR=-1>
	class OpticalLayerArray:
		public OpticalLayerArrayROP<NSTOKES>,
		public AzimuthDependencyCascade
	{
	public:
		OpticalLayerArray(const PersistentConfiguration<NSTOKES, CNSTR>& config,
			int wavelidx,
			const std::vector<LineOfSight>& los,
			std::unique_ptr<BRDF_Base> brdf,
			const sasktran_disco_lowlevel::Atmosphere& atmosphere,
			const sasktran_disco_lowlevel::WeightingFunctions* weightingfunctions,
			const GeometryLayerArray<NSTOKES, CNSTR>& geometry_layers,
            const ThreadData<NSTOKES, CNSTR>& thread_data
		);

        OpticalLayerArray(const PersistentConfiguration<NSTOKES, CNSTR>& config,
                          const std::vector<LineOfSight>& los,
                          const GeometryLayerArray<NSTOKES, CNSTR>& geometry_layers,
                          const ThreadData<NSTOKES, CNSTR>& thread_data
        );

        void set_optical(int wavelidx,
                         std::unique_ptr<BRDF_Base> brdf,
                         const sasktran_disco_lowlevel::Atmosphere& atmosphere,
                         const sasktran_disco_lowlevel::WeightingFunctions* weightingfunctions
                         );

        OpticalLayerArray(const PersistentConfiguration<NSTOKES, CNSTR>& config,
                          int wavelidx,
                          const std::vector<LineOfSight>& los,
                          std::unique_ptr<BRDF_Base> brdf,
                          const GeometryLayerArray<NSTOKES, CNSTR>& geometry_layers,
                          const sasktran2::atmosphere::Atmosphere<NSTOKES>& atmosphere,
                          const sasktran2::Config& sk_config
        );

		// Configures a the array for a test. This just skips pulling physical 
		// atmospheric properties.
		void configureTest(const PersistentConfiguration<NSTOKES, CNSTR>& config, const std::vector<testing::TestLayer<NSTOKES>>& testlayers);

		// Computes the reflected intensity for the given order of the azimuth
		// expansion.
		//
		// Performs Spurr2001 (C.6) and (C.21) and stores the results in m_reflected
		void computeReflectedIntensities(AEOrder m, const sasktran_disco::LineOfSight& los);

		// Returns the mutable layer at index lidx.
		inline OpticalLayer<NSTOKES, CNSTR>& operator[](uint lidx) {
			return *m_layers[lidx];
		}
		// Returns the non-mutable layer at index lidx.
		inline const OpticalLayer<NSTOKES, CNSTR>& operator[](uint lidx) const {
			return *m_layers[lidx];
		}
		// Returns a mutable reference to the bottom layer
		inline OpticalLayer<NSTOKES, CNSTR>& bottom() {
			return *m_layers.back();
		}
		// Returns a non-mutable reference to the bottom layer
		inline const OpticalLayer<NSTOKES, CNSTR>& bottom() const
		{
			return *m_layers.back();
		}
		// Returns the mutable layer at index lidx.
		inline OpticalLayer<NSTOKES, CNSTR>& layer(uint lidx)
		{
			return *m_layers[lidx];
		}
		// Returns the non-mutable layer at index lidx.
		inline const OpticalLayer<NSTOKES, CNSTR>& layer(uint lidx) const
		{
			return *m_layers[lidx];
		}
		// Returns a mutable reference to the top layer
		inline OpticalLayer<NSTOKES, CNSTR>& top()
		{
			return *m_layers[0];
		}
		// Returns a non-mutable reference to the bottom layer
		inline const OpticalLayer<NSTOKES, CNSTR>& top() const
		{
			return *m_layers[0];
		}

		inline double directIntensityTOA() const {
			return m_direct_toa;
		}

		inline const Albedo& albedo(AEOrder m) const {
			return m_albedo[m];
		}

		// Returns the number of layers.
		inline uint numLayers() const {
			return this->M_NLYR;
		}

		inline const InputDerivatives<NSTOKES>& inputDerivatives() const {
			return m_input_derivatives;
		}

		// Returns the total reflected intensity in the direction of the 
		// corresponding los.
		inline const Radiance<NSTOKES>& reflectedIntensity(AEOrder m, const LineOfSight& los) {
			if(!m_reflection_computed[m][los.unsorted_index]) {
				computeReflectedIntensities(m, los);
			}
			return m_ground_reflection[m][los.unsorted_index];
		}

		// Returns an iterator which will iterator in the specified direction 
		// to the given terminal optical depth.
		template<Propagating dir>
		inline OpticalLayerArrayIterator<dir, NSTOKES, CNSTR> iteratorTo(double terminal) {
			return OpticalLayerArrayIterator<dir, NSTOKES, CNSTR>(*this, terminal);
		}

		// Returns an iterator which will iterate across all layers in the 
		// specified direction.
		template<Propagating dir>
		inline OpticalLayerArrayIterator<dir, NSTOKES, CNSTR> iteratorAcross() {
			if(dir == Propagating::UP) {
				return iteratorTo<dir>(0);
			} else {
				return iteratorTo<dir>(bottom().opticalDepth(Location::FLOOR));
			}
		}

		inline OpticalLayer<NSTOKES, CNSTR>* layerAt(double optical_depth) {
			for(LayerIndex p = 0; p < this->M_NLYR; ++p) {
				if(layer(p).opticalDepth(Location::FLOOR) >= optical_depth) return &layer(p);
			}
			return nullptr;
		}

        inline const OpticalLayer<NSTOKES, CNSTR>* layerAtAltitude(double altitude) const {
            // Binary search since this can be slow and called many times

            int lower = 0;
            int upper = this->M_NLYR - 1;

            int retindex = 0;
            while(true) {
                if(lower == upper) {
                    retindex = upper;
                    break;
                }
                if(upper - lower == 1) {
                    retindex = layer(lower).altitude(Location::FLOOR) <= altitude ? lower : upper;
                    break;
                }
                int checkIndex = (lower + upper) / 2;

                if(layer(checkIndex).altitude(Location::FLOOR) > altitude) {
                    lower = checkIndex;
                } else {
                    upper = checkIndex;
                }
            }
            return &layer(retindex);
        }

        // TODO: This is just temporary
        void copyLegendre(VectorDim1<sasktran_disco::LegendreCoefficient<NSTOKES>>& container, const Eigen::Matrix<double, Eigen::Dynamic, 6>& legendre);

		double triangleFragmentArea(double h_lower, double h_upper, double h_center, double w_above, double w_below) {
			// area above center
			double left = (h_center + w_above) - h_upper;
			left = (std::min)(left, w_above);
			left = (std::max)(left, 0.0);

			double right = (h_center + w_above) - h_lower;
			right = (std::min)(right, w_above);
			right = (std::max)(right, 0.0);

			double top_area;
			if (w_above > 0) {
				top_area = 0.5 * 1.0 / w_above * (right*right - left * left);
			}
			else {
				top_area = 0;
			}


			// aread below center
			left = h_upper - (h_center - w_below);
			left = (std::min)(left, w_below);
			left = (std::max)(left, 0.0);

			right = h_lower - (h_center - w_below);
			right = (std::min)(right, w_below);
			right = (std::max)(right, 0.0);

			double bot_area;
			if (w_below > 0) {
				bot_area = 0.5 * 1 / w_below * (left*left - right * right);
			}
			else {
				bot_area = 0.0;
			}

			return top_area + bot_area;
		}

		inline const OpticalLayer<NSTOKES, CNSTR>* layerAt(double optical_depth) const {
			for(LayerIndex p = 0; p < this->M_NLYR; ++p) {
				if(layer(p).opticalDepth(Location::FLOOR) >= optical_depth) return &layer(p);
			}
			return nullptr;
		}
		
		inline double opticalDepthAt(double altitude) const {
            const auto& layer = layerAtAltitude(altitude);

            if(altitude > layer->altitude(Location::CEILING)) {
                return 0.0;
            } else {
                double f = 1 - (layer->altitude(Location::CEILING) - altitude) / (layer->altitude(Location::CEILING) - layer->altitude(Location::FLOOR));

                return layer->opticalDepth(Location::FLOOR) - layer->opticalDepth(Location::INSIDE) * f;
            }
		}

		inline double altitudeAt(double optical_depth) const {
            // TODO: Calculate directly from the layer quantities
            return 0.0;
		}

	protected:
		void configureTransmission();
        void assignLegendreDerivative(std::vector<LegendreCoefficient<NSTOKES>>& d_legendre,
                                      const Eigen::Matrix<double, Eigen::Dynamic, 6>& species_legendre,
                                      const std::vector<LegendreCoefficient<NSTOKES>>& layer_legendre,
                                      double species_ssa,
                                      double layer_ssa,
                                      double thickness
                                      ) const;
        void assignLegendreDerivativeTest(std::vector<LegendreCoefficient<NSTOKES>>& d_legendre,
                                          const std::vector<double>& legendre) const;
	protected:
		// Members
		VectorDim1<std::unique_ptr<OpticalLayer<NSTOKES, CNSTR>>>	m_layers;
		double m_direct_toa;
        size_t                                      m_num_los;

		InputDerivatives<NSTOKES>&					m_input_derivatives;
        size_t                                      m_wavel_index;

		VectorDim2<Radiance<NSTOKES>>				m_ground_reflection;
		VectorDim2<bool>							m_reflection_computed;
		
		const Eigen::VectorXd*						m_cell_depths;

		// Pseudo-spherical beam transmittances
		Eigen::MatrixXd								m_chapman_factors;

		const PersistentConfiguration<NSTOKES, CNSTR>&		m_config;

		AlbedoExpansion			                    m_albedo;

		bool										m_include_direct_bounce;
		
		double										m_surfaceemission;

	};

	// STL random access iterator with additional layer of functionality for 
	// iterating over a range of a OpticalLayerArray.
	template<Propagating dir, int NSTOKES, int CNSTR>
	class OpticalLayerArrayIterator
	{
	public:
		// Constructors
		OpticalLayerArrayIterator() noexcept :
			m_layers(nullptr),
			m_current_layer_index(-1),
			m_terminal_optical_depth(std::nan("1"))
		{
			// empty
		}
		OpticalLayerArrayIterator(const OpticalLayerArrayIterator& other) noexcept = default;
		OpticalLayerArrayIterator(OpticalLayerArrayIterator&& other) noexcept = default;
		OpticalLayerArrayIterator& operator=(OpticalLayerArrayIterator&& other) noexcept = default;
		OpticalLayerArrayIterator& operator=(const OpticalLayerArrayIterator& other) noexcept = default;

		// STL type info
		typedef int difference_type;
		typedef OpticalLayer<NSTOKES, CNSTR> value_type;
		typedef OpticalLayer<NSTOKES, CNSTR>* pointer;
		typedef OpticalLayer<NSTOKES, CNSTR>& reference;
		typedef std::random_access_iterator_tag iterator_category;

		// Comparison operators
		bool operator==(const OpticalLayerArrayIterator& other) const {
			return m_layers == other.m_layers
				&& m_current_layer_index == other.m_current_layer_index
				&& m_terminal_optical_depth == other.m_terminal_optical_depth;
		}
		inline bool operator!=(const OpticalLayerArrayIterator& other) const {
			return *this == other;
		}
		bool operator<(const OpticalLayerArrayIterator& other) const {
			if(isDownwelling()) {
				return m_current_layer_index < other.m_current_layer_index;
			} else {
				return m_current_layer_index > other.m_current_layer_index;
			}
		}
		bool operator<=(const OpticalLayerArrayIterator& other) const {
			return *this < other || *this == other;
		}
		bool operator>(const OpticalLayerArrayIterator& other) const {
			return !(*this <= other);
		}
		bool operator>=(const OpticalLayerArrayIterator& other) const {
			return !(*this < other);
		}

		// Dereference operators
		OpticalLayer<NSTOKES, CNSTR>& operator*() {						// needs to be removed for const
			return layer();
		}
		const OpticalLayer<NSTOKES, CNSTR>& operator*() const {
			return layer();
		}
		OpticalLayer<NSTOKES, CNSTR>& operator->() {					// needs to be removed for const
			return layer();
		}
		const OpticalLayer<NSTOKES, CNSTR>& operator->() const {
			return layer();
		}
		inline OpticalLayer<NSTOKES, CNSTR>& operator[](int offset) {
			return m_layers->layer(m_current_layer_index + offset);
		}
		inline const OpticalLayer<NSTOKES, CNSTR>& operator[](int offset) const {
			return m_layers->layer(m_current_layer_index + offset);
		}

		// Increment/decrement operators
		inline OpticalLayerArrayIterator& operator++() {
			if(isDownwelling()) ++m_current_layer_index;
			else --m_current_layer_index;
			return *this;
		}
		inline OpticalLayerArrayIterator& operator++(int) {
			if(isDownwelling()) ++m_current_layer_index;
			else --m_current_layer_index;
			return *this;
		}
		inline OpticalLayerArrayIterator& operator--() {
			if(isDownwelling()) --m_current_layer_index;
			else ++m_current_layer_index;
			return *this;
		}
		inline OpticalLayerArrayIterator& operator--(int) {
			if(isDownwelling()) --m_current_layer_index;
			else ++m_current_layer_index;
			return *this;

		}
		inline OpticalLayerArrayIterator& operator+=(int offset) {
			if(isDownwelling()) m_current_layer_index += offset;
			else m_current_layer_index -= offset;
			return *this;
		}
		inline OpticalLayerArrayIterator& operator-=(int offset) {
			return operator+=(-offset);
		}

		// Arithmetic operators
		inline OpticalLayerArrayIterator operator+(int offset) const {
			return OpticalLayerArrayIterator(*this) += offset;
		}
		inline OpticalLayerArrayIterator operator-(int offset) const {
			return  OpticalLayerArrayIterator(*this) -= offset;
		}
		inline int operator+(const OpticalLayerArrayIterator& other) const {
			return m_current_layer_index + other.m_current_layer_index;
		}
		inline int operator-(const OpticalLayerArrayIterator& other) const {
			return m_current_layer_index - other.m_current_layer_index;
		}

	public: // Interface
		// Constructs a iterator over all layers.
		OpticalLayerArrayIterator(OpticalLayerArray<NSTOKES, CNSTR>& layers) {
			m_layers = &layers;
			if(isDownwelling()) {
				m_current_layer_index = 0;
			} else {
				m_current_layer_index = layers.numLayers() - 1;
			}
			if(isDownwelling()) {
				m_terminal_optical_depth = layers.bottom().opticalDepth(Location::FLOOR);
			} else {
				m_terminal_optical_depth = layers.top().opticalDepth(Location::CEILING);
			}
		}
		// Constructs an iterator which terminates at terminal_depth.
		OpticalLayerArrayIterator(OpticalLayerArray<NSTOKES, CNSTR>& layers, double terminal_depth):
			OpticalLayerArrayIterator(layers)
		{
			m_terminal_optical_depth = terminal_depth;
		}

		// Returns non-mutable OpticalLayer reference
		inline const OpticalLayer<NSTOKES, CNSTR>& layer() const {
			return m_layers->layer(m_current_layer_index);
		}
		// Returns mutable OpticalLayer reference
		inline OpticalLayer<NSTOKES, CNSTR>& layer() {
			return m_layers->layer(m_current_layer_index);
		}
		// Returns true if this iterator is still accessing a valid range of the 
		// OpticalLayerArray false otherwise. It is important to not perform any
		// operations an an iterator when isValid() is false. This indicates 
		// that the iterator has expired.
		inline bool isValid() const {
			bool ok = true;
			ok &= m_current_layer_index < static_cast<int>(m_layers->numLayers()) && m_current_layer_index >= 0;
			if(!ok) return ok;
			if(isUpwelling()) {
				ok &= layer().opticalDepth(Location::FLOOR) > m_terminal_optical_depth;
			} else {
				ok &= layer().opticalDepth(Location::CEILING) < m_terminal_optical_depth;
			}
			return ok;
		}

		// Return the entry point to a layer (upwelling is FLOOR, down welling 
		// is CEILING)
		inline Location entry() const {
			return isUpwelling() ? Location::FLOOR : Location::CEILING;
		}
		// Returns the exit point to a layer (see entry())
		inline Location exit() const {
			if(m_terminal_optical_depth > layer().opticalDepth(Location::FLOOR) && m_terminal_optical_depth < layer().opticalDepth(Location::CEILING)) {
				return Location::INSIDE;
			} else {
				return isUpwelling() ? Location::CEILING : Location::FLOOR;
			}
		}
		// Returns the optical depth at entry to the current layer
		inline double entryOpticalDepth() const {
			return isUpwelling() ?
				layer().opticalDepth(Location::FLOOR) : layer().opticalDepth(Location::CEILING);
		}
		// Returns the optical depth at the exit of the current layer
		inline double exitOpticalDepth() const {
			if(exit() == Location::INSIDE) return m_terminal_optical_depth;
			else {
				return isUpwelling() ?
					layer().opticalDepth(Location::CEILING) : layer().opticalDepth(Location::FLOOR);
			}
		}
		// Returns true if the iterator is upwelling false otherwise.
		constexpr inline bool isUpwelling() const {
			return dir == Propagating::UP;
		}
		// Returns true if the iterator is down welling false otherwise.
		constexpr inline bool isDownwelling() const {
			return dir == Propagating::DOWN;
		}
		// Returns the attenuation factor for this layer at the given coszentih
		inline void attenuationFactor(double coszenith, double losod, double losaltitude, Dual<double>& factors) const {
			// TODO: Layer dual?
			size_t layerStart = m_layers->inputDerivatives().layerStartIndex(layer().index());
			size_t numDeriv = m_layers->inputDerivatives().numDerivativeLayer(layer().index());

			double exitod = std::max(exitOpticalDepth(), losod);

			factors.value = exp(-std::abs((entryOpticalDepth() - exitod) / coszenith));

			// We do not have cross derivatives
            // But we store as a full dual for some reason so we have to zero everything
            factors.deriv.setZero();
			for (uint i = 0; i < numDeriv; ++i) {
				const auto& deriv = m_layers->inputDerivatives().layerDerivatives()[i + layerStart];
				
				for (uint j = 0; j < deriv.group_and_triangle_fraction.size(); ++j) {
					double pertfactor = deriv.group_and_triangle_fraction[j].second * deriv.extinctions[j];
					if (exitod != exitOpticalDepth()) {
						double altitude_floor = m_layers->altitudeAt(entryOpticalDepth());

						double lower_triangle = m_layers->triangleFragmentArea(altitude_floor, losaltitude,
							std::get<0>(deriv.alt_and_widths[j]),
							std::get<1>(deriv.alt_and_widths[j]),
							std::get<2>(deriv.alt_and_widths[j]));
                        // Account for rounding errors
                        if(losaltitude <= altitude_floor) {
                            lower_triangle = 0.0;
                        }
                        // Can rarely happen for odd cases, but means that the area was 0
                        if(lower_triangle != lower_triangle) {
                            lower_triangle = 0.0;
                        }
                        // Account for extinction in cm2
						pertfactor = lower_triangle * deriv.extinctions[j] * 100;
					}
					double added_od =  pertfactor * deriv.d_optical_depth;

					factors.deriv[deriv.group_and_triangle_fraction[j].first] = -1.0*factors.value * added_od / coszenith;
				}
			}

		}
		// Returns mutable pointer to the layer
		inline OpticalLayer<NSTOKES, CNSTR>* ptr() {
			return &layer();
		}

		// Returns non-mutable pointer to the layer
		inline const OpticalLayer<NSTOKES, CNSTR>* ptr() const {
			return &layer();
		}

	private:
		double m_terminal_optical_depth;
		int m_current_layer_index;
		OpticalLayerArray<NSTOKES, CNSTR>* m_layers;

	};
}
