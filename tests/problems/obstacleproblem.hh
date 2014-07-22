/*
  Copyright (C) 2008-2013 by Andreas Lauser

  This file is part of the Open Porous Media project (OPM).

  OPM is free software: you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation, either version 2 of the License, or
  (at your option) any later version.

  OPM is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with OPM.  If not, see <http://www.gnu.org/licenses/>.
*/
/*!
 * \file
 *
 * \copydoc Ewoms::ObstacleProblem
 */
#ifndef EWOMS_OBSTACLE_PROBLEM_HH
#define EWOMS_OBSTACLE_PROBLEM_HH

#include <ewoms/models/ncp/ncpproperties.hh>

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedBrooksCorey.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/heatconduction/Somerton.hpp>

#include <dune/grid/yaspgrid.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <string>
#include <iostream>

namespace Ewoms {
template <class TypeTag>
class ObstacleProblem;
}

namespace Opm {
namespace Properties {
NEW_TYPE_TAG(ObstacleBaseProblem);

// Set the grid type
SET_TYPE_PROP(ObstacleBaseProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ObstacleBaseProblem, Problem, Ewoms::ObstacleProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(ObstacleBaseProblem, FluidSystem,
              Opm::FluidSystems::H2ON2<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Set the material Law
SET_PROP(ObstacleBaseProblem, MaterialLaw)
{
private:
    // define the material law
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef Opm::TwoPhaseMaterialTraits<Scalar,
                                        /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                        /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx>
    MaterialTraits;

    typedef Opm::LinearMaterial<MaterialTraits> EffMaterialLaw;

public:
    typedef Opm::EffToAbsLaw<EffMaterialLaw> type;
};

// Set the heat conduction law
SET_PROP(ObstacleBaseProblem, HeatConductionLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // define the material law parameterized by absolute saturations
    typedef Opm::Somerton<FluidSystem, Scalar> type;
};

// Enable gravity
SET_BOOL_PROP(ObstacleBaseProblem, EnableGravity, true);

// The default for the end time of the simulation
SET_SCALAR_PROP(ObstacleBaseProblem, EndTime, 1e4);

// The default for the initial time step size of the simulation
SET_SCALAR_PROP(ObstacleBaseProblem, InitialTimeStepSize, 250);

// The default DGF file to load
SET_STRING_PROP(ObstacleBaseProblem, GridFile, "./data/obstacle_24x16.dgf");
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup TestProblems
 *
 * \brief Problem where liquid water is first stopped by a
 *        low-permeability lens and then seeps though it.
 *
 * Liquid water is injected by using of a free-flow condition on the
 * lower right of the domain. This water level then raises until
 * hydrostatic pressure is reached. On the left of the domain, a
 * rectangular obstacle with \f$10^3\f$ lower permeability than the
 * rest of the domain first stops the for a while until it seeps
 * through it.
 *
 * The domain is sized 60m times 40m and consists of two media, a
 * moderately permeable soil (\f$ K_0=10e-12 m^2\f$) and an obstacle
 * at \f$[10; 20]m \times [0; 35]m \f$ with a lower permeablility of
 * \f$ K_1=K_0/1000\f$.
 *
 * Initially the whole domain is filled by nitrogen, the temperature
 * is \f$20^\circ C\f$ for the whole domain. The gas pressure is
 * initially 1 bar, at the inlet of the liquid water on the right side
 * it is 2 bar.
 *
 * The boundary is no-flow except on the lower 10 meters of the left
 * and the right boundary where a free flow condition is assumed.
 */
template <class TypeTag>
class ObstacleProblem : public GET_PROP_TYPE(TypeTag, BaseProblem)
{
    typedef typename GET_PROP_TYPE(TypeTag, BaseProblem) ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, EqVector) EqVector;
    typedef typename GET_PROP_TYPE(TypeTag, RateVector) RateVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryRateVector) BoundaryRateVector;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLawParams) MaterialLawParams;
    typedef typename GET_PROP_TYPE(TypeTag, HeatConductionLaw) HeatConductionLaw;
    typedef typename HeatConductionLaw::Params HeatConductionLawParams;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        gasPhaseIdx = FluidSystem::gasPhaseIdx,
        liquidPhaseIdx = FluidSystem::liquidPhaseIdx,
        H2OIdx = FluidSystem::H2OIdx,
        N2Idx = FluidSystem::N2Idx
    };

    typedef Dune::FieldVector<typename GridView::ctype, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, numPhases> PhaseVector;
    typedef Dune::FieldMatrix<Scalar, dimWorld, dimWorld> DimMatrix;
    typedef typename GET_PROP_TYPE(TypeTag, Simulator) Simulator;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    ObstacleProblem(Simulator &simulator)
        : ParentType(simulator)
    {
        eps_ = 1e-6;
        temperature_ = 273.15 + 25; // -> 25°C

        // initialize the tables of the fluid system
        Scalar Tmin = temperature_ - 1.0;
        Scalar Tmax = temperature_ + 1.0;
        int nT = 3;

        Scalar pmin = 1.0e5 * 0.75;
        Scalar pmax = 2.0e5 * 1.25;
        int np = 1000;

        FluidSystem::init(Tmin, Tmax, nT, pmin, pmax, np);

        // intrinsic permeabilities
        coarseK_ = this->toDimMatrix_(1e-12);
        fineK_ = this->toDimMatrix_(1e-15);

        // the porosity
        finePorosity_ = 0.3;
        coarsePorosity_ = 0.3;

        // residual saturations
        fineMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.0);
        fineMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);
        coarseMaterialParams_.setResidualSaturation(liquidPhaseIdx, 0.0);
        coarseMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);

        // parameters for the linear law, i.e. minimum and maximum
        // pressures
        fineMaterialParams_.setPcMinSat(liquidPhaseIdx, 0.0);
        fineMaterialParams_.setPcMaxSat(liquidPhaseIdx, 0.0);
        coarseMaterialParams_.setPcMinSat(liquidPhaseIdx, 0.0);
        coarseMaterialParams_.setPcMaxSat(liquidPhaseIdx, 0.0);

        /*
        // entry pressures for Brooks-Corey
        fineMaterialParams_.setEntryPressure(5e3);
        coarseMaterialParams_.setEntryPressure(1e3);

        // Brooks-Corey shape parameters
        fineMaterialParams_.setLambda(2);
        coarseMaterialParams_.setLambda(2);
        */

        fineMaterialParams_.finalize();
        coarseMaterialParams_.finalize();

        // parameters for the somerton law of heat conduction
        computeHeatCondParams_(fineHeatCondParams_, finePorosity_);
        computeHeatCondParams_(coarseHeatCondParams_, coarsePorosity_);

        initFluidStates_();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        this->model().checkConservativeness();

        // Calculate storage terms of the individual phases
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            PrimaryVariables phaseStorage;
            this->model().globalPhaseStorage(phaseStorage, phaseIdx);

            if (this->gridView().comm().rank() == 0) {
                std::cout << "Storage in " << FluidSystem::phaseName(phaseIdx)
                          << "Phase: [" << phaseStorage << "]"
                          << "\n"  << std::flush;
            }
        }

        // Calculate total storage terms
        EqVector storage;
        this->model().globalStorage(storage);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage total: [" << storage << "]"
                      << "\n"  << std::flush;
        }
#endif // NDEBUG
    }

    /*!
     * \name Problem parameters
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::name
     */
    std::string name() const
    {
        std::ostringstream oss;
        oss << "obstacle"
            << "_" << Model::name();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     *
     * This problem simply assumes a constant temperature.
     */
    template <class Context>
    Scalar temperature(const Context &context, int spaceIdx, int timeIdx) const
    { return temperature_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix &intrinsicPermeability(const Context &context, int spaceIdx,
                                           int timeIdx) const
    {
        if (isFineMaterial_(context.pos(spaceIdx, timeIdx)))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams &materialLawParams(const Context &context,
                                               int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::heatCapacitySolid
     *
     * For this problem, we assume that the solid phase of the porous
     * medium is granite.
     */
    template <class Context>
    Scalar heatCapacitySolid(const Context &context, int spaceIdx,
                             int timeIdx) const
    {
        return 790     // specific heat capacity of granite [J / (kg K)]
               * 2700; // density of granite [kg/m^3]
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::heatConductionParams
     */
    template <class Context>
    const HeatConductionLawParams &
    heatConductionParams(const Context &context, int spaceIdx, int timeIdx) const
    {
        const GlobalPosition &pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineHeatCondParams_;
        return coarseHeatCondParams_;
    }

    //! \}

    /*!
     * \name Boundary conditions
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::boundary
     */
    template <class Context>
    void boundary(BoundaryRateVector &values, const Context &context,
                  int spaceIdx, int timeIdx) const
    {
        const auto &pos = context.pos(spaceIdx, timeIdx);

        if (onInlet_(pos))
            values.setFreeFlow(context, spaceIdx, timeIdx, inletFluidState_);
        else if (onOutlet_(pos))
            values.setFreeFlow(context, spaceIdx, timeIdx, outletFluidState_);
        else
            values.setNoFlow();
    }

    //! \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables &values, const Context &context, int spaceIdx,
                 int timeIdx) const
    {
        const auto &matParams = materialLawParams(context, spaceIdx, timeIdx);
        values.assignMassConservative(outletFluidState_, matParams);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector &rate, const Context &context, int spaceIdx,
                int timeIdx) const
    { rate = 0.0; }

    //! \}

private:
    /*!
     * \brief Returns whether a given global position is in the
     *        fine-permeability region or not.
     */
    bool isFineMaterial_(const GlobalPosition &pos) const
    { return 10 <= pos[0] && pos[0] <= 20 && 0 <= pos[1] && pos[1] <= 35; }

    bool onInlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x >= 60 - eps_ && y <= 10;
    }

    bool onOutlet_(const GlobalPosition &globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];
        return x < eps_ && y <= 10;
    }

    void initFluidStates_()
    {
        initFluidState_(inletFluidState_, coarseMaterialParams_,
                        /*isInlet=*/true);
        initFluidState_(outletFluidState_, coarseMaterialParams_,
                        /*isInlet=*/false);
    }

    template <class FluidState>
    void initFluidState_(FluidState &fs, const MaterialLawParams &matParams,
                         bool isInlet)
    {
        int refPhaseIdx;
        int otherPhaseIdx;

        // set the fluid temperatures
        fs.setTemperature(temperature_);

        if (isInlet) {
            // only liquid on inlet
            refPhaseIdx = liquidPhaseIdx;
            otherPhaseIdx = gasPhaseIdx;

            // set liquid saturation
            fs.setSaturation(liquidPhaseIdx, 1.0);

            // set pressure of the liquid phase
            fs.setPressure(liquidPhaseIdx, 2e5);

            // set the liquid composition to pure water
            fs.setMoleFraction(liquidPhaseIdx, N2Idx, 0.0);
            fs.setMoleFraction(liquidPhaseIdx, H2OIdx, 1.0);
        }
        else {
            // elsewhere, only gas
            refPhaseIdx = gasPhaseIdx;
            otherPhaseIdx = liquidPhaseIdx;

            // set gas saturation
            fs.setSaturation(gasPhaseIdx, 1.0);

            // set pressure of the gas phase
            fs.setPressure(gasPhaseIdx, 1e5);

            // set the gas composition to 99% nitrogen and 1% steam
            fs.setMoleFraction(gasPhaseIdx, N2Idx, 0.99);
            fs.setMoleFraction(gasPhaseIdx, H2OIdx, 0.01);
        };

        // set the other saturation
        fs.setSaturation(otherPhaseIdx, 1.0 - fs.saturation(refPhaseIdx));

        // calulate the capillary pressure
        PhaseVector pC;
        MaterialLaw::capillaryPressures(pC, matParams, fs);
        fs.setPressure(otherPhaseIdx, fs.pressure(refPhaseIdx)
                                      + (pC[otherPhaseIdx] - pC[refPhaseIdx]));

        // make the fluid state consistent with local thermodynamic
        // equilibrium
        typedef Opm::ComputeFromReferencePhase<Scalar, FluidSystem>
        ComputeFromReferencePhase;

        typename FluidSystem::ParameterCache paramCache;
        ComputeFromReferencePhase::solve(fs, paramCache, refPhaseIdx,
                                         /*setViscosity=*/false,
                                         /*setEnthalpy=*/false);
    }

    void computeHeatCondParams_(HeatConductionLawParams &params, Scalar poro)
    {
        Scalar lambdaWater = 0.6;
        Scalar lambdaGranite = 2.8;

        Scalar lambdaWet = std::pow(lambdaGranite, (1 - poro))
                           * std::pow(lambdaWater, poro);
        Scalar lambdaDry = std::pow(lambdaGranite, (1 - poro));

        params.setFullySaturatedLambda(gasPhaseIdx, lambdaDry);
        params.setFullySaturatedLambda(liquidPhaseIdx, lambdaWet);
        params.setVacuumLambda(lambdaDry);
    }

    DimMatrix coarseK_;
    DimMatrix fineK_;

    Scalar coarsePorosity_;
    Scalar finePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    HeatConductionLawParams fineHeatCondParams_;
    HeatConductionLawParams coarseHeatCondParams_;

    Opm::CompositionalFluidState<Scalar, FluidSystem> inletFluidState_;
    Opm::CompositionalFluidState<Scalar, FluidSystem> outletFluidState_;

    Scalar temperature_;
    Scalar eps_;
};
} // namespace Ewoms

#endif
