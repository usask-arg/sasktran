#pragma once

#include <sasktran2/internal_common.h>
#include <sasktran2/dual.h>

namespace sasktran2 {

    /** Essentially a pure virtual void class to interface with SWIG, removing the NSTOKES template.
     *
     */
    class OutputInterface {
    public:
        virtual ~OutputInterface() {}
    };


    /** Base class for the output container.  Provides storage for the output values, as well as defines what quantities
     *  are actually output.  The Sasktran2 engine both solves the radiative transfer equation to get the full radiance
     *  field, \f$I(location, solid angle, wavelength)\f$ and integrates this along the specified LOS vectors to get
     *  \f$(I_{los_i}(wavelength)\f$.  As well as corresponding derivatives.
     *
     *  The user output is going to be a function of these two parameters.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES>
    class Output : public OutputInterface {
    protected:
        int m_nlos;
        int m_nwavel;
        int m_nderiv;
    public:
        Output() {};
        virtual ~Output() {}

        /** Method called by the Sasktran2 engine which specifies the native number of lines of sight, wavelength batches,
         *  and derivatives.  These are then internally stored for derived classes to access.
         *
         * @param nlos
         * @param nwavel
         * @param nderiv
         */
        virtual void resize(int nlos, int nwavel, int nderiv) {
            m_nlos = nlos;
            m_nwavel = nwavel;
            m_nderiv = nderiv;
        }

        /** Method the Sasktran2 engine calls for each integrated line of sight/wavelength
         *
         * @param radiance The final calculated radiance and corresponding derivatives
         * @param losidx The index of this line of sight
         * @param wavelidx The index of this wavelength
         */
        virtual void assign(const sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& radiance, int losidx, int wavelidx) = 0;

        /**
         *
         * @return The engine number of lines of sight
         */
        int num_los() const { return m_nlos; }

        /**
         *
         * @return The engine number of wavelengths
         */
        int num_wavel() const { return m_nwavel; }

        /**
         *
         * @return The engine number of derivatives
         */
        int num_deriv() const { return m_nderiv; }
    };

    /** An idealized output container where only the line of sight radiances are stored, and are stored on the native
     *  model calculation grid.
     *
     * @tparam NSTOKES
     */
    template <int NSTOKES>
    class OutputIdealDense : public Output<NSTOKES> {
    private:
        sasktran2::Dual<double, sasktran2::dualstorage::dense> m_radiance; /**< Internal storage */
    public:
        OutputIdealDense() {};

        void resize(int nlos, int nwavel, int nderiv);

        void assign(const sasktran2::Dual<double, sasktran2::dualstorage::dense, NSTOKES>& radiance, int losidx, int wavelidx);

        /**
         *
         * @return The stored radiance container
         */
        sasktran2::Dual<double, sasktran2::dualstorage::dense>& radiance() { return m_radiance; }
    };
}