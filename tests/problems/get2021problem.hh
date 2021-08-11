// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
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

  Consult the COPYING file in the top-level source directory of this
  module for the precise wording of the license and the list of
  copyright holders.
*/
/*!
 * \file
 *
 * \copydoc Opm::Get2021Problem
 */
#ifndef EWOMS_CO2_GET2021_PROBLEM_HH
#define EWOMS_CO2_GET2021_PROBLEM_HH

#include <opm/models/immiscible/immisciblemodel.hh>
#include <opm/simulators/linalg/parallelamgbackend.hh>

#include <opm/material/fluidsystems/H2ON2FluidSystem.hpp>
#include <opm/material/fluidsystems/BrineCO2FluidSystem.hpp>
#include <opm/material/fluidstates/CompositionalFluidState.hpp>
#include <opm/material/fluidstates/ImmiscibleFluidState.hpp>
#include <opm/material/constraintsolvers/ComputeFromReferencePhase.hpp>
#include <opm/material/fluidmatrixinteractions/LinearMaterial.hpp>
#include <opm/material/fluidmatrixinteractions/RegularizedDynamicWa.hpp>
#include <opm/material/fluidmatrixinteractions/EffToAbsLaw.hpp>
#include <opm/material/fluidmatrixinteractions/MaterialTraits.hpp>
#include <opm/material/thermal/SomertonThermalConductionLaw.hpp>
#include <opm/material/thermal/ConstantSolidHeatCapLaw.hpp>
#include <opm/material/binarycoefficients/Brine_CO2.hpp>
#include <opm/material/common/UniformTabulated2DFunction.hpp>
#include <opm/material/common/Unused.hpp>

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <sstream>
#include <iostream>
#include <string>

namespace Opm {
//! \cond SKIP_THIS
template <class TypeTag>
class Get2021Problem;

namespace Get2021 {
#include <opm/material/components/co2tables.inc>
}
//! \endcond
}

namespace Opm::Properties {

namespace TTag {
struct Get2021BaseProblem {};
}

// declare the CO2 wa problem specific property tags
template<class TypeTag, class MyTypeTag>
struct FluidSystemPressureLow { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemPressureHigh { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemNumPressure { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemTemperatureLow { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemTemperatureHigh { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FluidSystemNumTemperature { using type = UndefinedProperty; };

template<class TypeTag, class MyTypeTag>
struct MaxDepth { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct Temperature { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct SimulationName { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct KC { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct K { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct PhiC { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct Phi { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct CiC { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct CfC { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct S0w { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct X0n { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct FineLayerBottom { using type = UndefinedProperty; };
template<class TypeTag, class MyTypeTag>
struct InjRate { using type = UndefinedProperty; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::Get2021BaseProblem> { using type = Dune::YaspGrid<3>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::Get2021BaseProblem>
{ using type = Opm::Get2021Problem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::Get2021BaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using CO2Tables = Opm::Get2021::CO2Tables;

public:
    using type = Opm::BrineCO2FluidSystem<Scalar, CO2Tables>;
};

// Set the material Law
template<class TypeTag>
struct MaterialLaw<TypeTag, TTag::Get2021BaseProblem>
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Traits = Opm::TwoPhaseMaterialTraits<Scalar,
                                               /*wettingPhaseIdx=*/FluidSystem::liquidPhaseIdx,
                                               /*nonWettingPhaseIdx=*/FluidSystem::gasPhaseIdx>;

    // define the material law which is parameterized by effective
    // saturations
    using EffMaterialLaw = Opm::RegularizedDynamicWa<Traits>;

public:
    // define the material law parameterized by absolute saturations
    using type = Opm::EffToAbsLaw<EffMaterialLaw>;
};

// Set the thermal conduction law
template<class TypeTag>
struct ThermalConductionLaw<TypeTag, TTag::Get2021BaseProblem>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    // define the material law parameterized by absolute saturations
    using type = Opm::SomertonThermalConductionLaw<FluidSystem, Scalar>;
};

// set the energy storage law for the solid phase
template<class TypeTag>
struct SolidEnergyLaw<TypeTag, TTag::Get2021BaseProblem>
{ using type = Opm::ConstantSolidHeatCapLaw<GetPropType<TypeTag, Properties::Scalar>>; };

// Use the algebraic multi-grid linear solver for this problem
template<class TypeTag>
struct LinearSolverSplice<TypeTag, TTag::Get2021BaseProblem> { using type = TTag::ParallelAmgLinearSolver; };

// Write the Newton convergence behavior to disk?
template<class TypeTag>
struct NewtonWriteConvergence<TypeTag, TTag::Get2021BaseProblem> { static constexpr bool value = false; };

// Enable molecular diffusion for this problem
template<class TypeTag>
struct EnableDiffusion<TypeTag, TTag::Get2021BaseProblem> { static constexpr bool value = true; };

// Enable gravity
template<class TypeTag>
struct EnableGravity<TypeTag, TTag::Get2021BaseProblem> { static constexpr bool value = true; };

// set the defaults for the problem specific properties
template<class TypeTag>
struct FluidSystemPressureLow<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 3e7;
};
template<class TypeTag>
struct FluidSystemPressureHigh<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 4e7;
};
template<class TypeTag>
struct FluidSystemNumPressure<TypeTag, TTag::Get2021BaseProblem> { static constexpr int value = 100; };
template<class TypeTag>
struct FluidSystemTemperatureLow<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 290;
};
template<class TypeTag>
struct FluidSystemTemperatureHigh<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 500;
};
template<class TypeTag>
struct FluidSystemNumTemperature<TypeTag, TTag::Get2021BaseProblem> { static constexpr int value = 100; };

template<class TypeTag>
struct MaxDepth<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 2500;
};
template<class TypeTag>
struct KC<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1E-16;
};
template<class TypeTag>
struct K<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1E-10;
};
template<class TypeTag>
struct PhiC<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.2;
};
template<class TypeTag>
struct Phi<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.2;
};
template<class TypeTag>
struct CiC<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e4;
};
template<class TypeTag>
struct CfC<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 1e2;
};
template<class TypeTag>
struct S0w<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.5;
};
template<class TypeTag>
struct X0n<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 5e-3;
};
template<class TypeTag>
struct FineLayerBottom<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 70.;
};
template<class TypeTag>
struct InjRate<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = -1e-4;
};
template<class TypeTag>
struct Temperature<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 293.15;
};
template<class TypeTag>
struct SimulationName<TypeTag, TTag::Get2021BaseProblem> { static constexpr auto value = "wa"; };

// The default for the end time of the simulation
template<class TypeTag>
struct EndTime<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 365*86400;
};

// The default for the initial time step size of the simulation
template<class TypeTag>
struct InitialTimeStepSize<TypeTag, TTag::Get2021BaseProblem>
{
    using type = GetPropType<TypeTag, Scalar>;
    static constexpr type value = 0.01;
};

// The default DGF file to load
template<class TypeTag>
struct GridFile<TypeTag, TTag::Get2021BaseProblem> { static constexpr auto value = "get2021.dgf"; };

} // namespace Opm::Properties

namespace Opm {
/*!
 * \ingroup TestProblems
 *
 * \brief Problem where \f$CO_2\f$ is injected under a low permeable
 *        layer at a depth of 2700m.
 *
 * The domain is sized 4000m times 1m times 100m and consists of two layers, one
 * which is permeable (\f$K = 10^{-10}\;m^2\f$) for \f$ z <
 * 70\; m\f$ and one with a lower intrinsic permeablility (\f$
 * K=10^{-16}\;m^2\f$) in the rest of the domain.
 *
 * \f$CO_2\f$ gets injected by means of a forced-flow boundary
 * condition into water-filled aquifer, at the lower-left boundary and
 * migrates upwards due to buoyancy. It accumulates and eventually
 * enters the caprock by means of wettability alteration.
 *
 * The boundary conditions applied by this problem are no-flow
 * conditions on the top bottom and right boundaries and a free-flow
 * boundary condition on the right. For the free-flow condition,
 * hydrostatic pressure is assumed.
 */
template <class TypeTag>
class Get2021Problem : public GetPropType<TypeTag, Properties::BaseProblem>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseProblem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    enum { dim = GridView::dimension };
    enum { dimWorld = GridView::dimensionworld };

    // copy some indices for convenience
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    enum { numPhases = FluidSystem::numPhases };
    enum { gasPhaseIdx = FluidSystem::gasPhaseIdx };
    enum { liquidPhaseIdx = FluidSystem::liquidPhaseIdx };
    enum { CO2Idx = FluidSystem::CO2Idx };
    enum { BrineIdx = FluidSystem::BrineIdx };
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { contiCO2EqIdx = conti0EqIdx + CO2Idx };

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using BoundaryRateVector = GetPropType<TypeTag, Properties::BoundaryRateVector>;
    using MaterialLaw = GetPropType<TypeTag, Properties::MaterialLaw>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using MaterialLawParams = GetPropType<TypeTag, Properties::MaterialLawParams>;
    using ThermalConductionLaw = GetPropType<TypeTag, Properties::ThermalConductionLaw>;
    using SolidEnergyLawParams = GetPropType<TypeTag, Properties::SolidEnergyLawParams>;
    using ThermalConductionLawParams = typename ThermalConductionLaw::Params;

    using Toolbox = Opm::MathToolbox<Evaluation>;
    using CoordScalar = typename GridView::ctype;
    using GlobalPosition = Dune::FieldVector<CoordScalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:
    /*!
     * \copydoc Doxygen::defaultProblemConstructor
     */
    Get2021Problem(Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
     * \copydoc FvBaseProblem::finishInit
     */
    void finishInit()
    {
        ParentType::finishInit();

        eps_ = 1e-6;

        temperatureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow);
        temperatureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh);
        nTemperature_ = EWOMS_GET_PARAM(TypeTag, unsigned, FluidSystemNumTemperature);

        pressureLow_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureLow);
        pressureHigh_ = EWOMS_GET_PARAM(TypeTag, Scalar, FluidSystemPressureHigh);
        nPressure_ = EWOMS_GET_PARAM(TypeTag, unsigned, FluidSystemNumPressure);

        maxDepth_ = EWOMS_GET_PARAM(TypeTag, Scalar, MaxDepth);
        temperature_ = EWOMS_GET_PARAM(TypeTag, Scalar, Temperature);
        kC_ = EWOMS_GET_PARAM(TypeTag, Scalar, KC);
        k_ = EWOMS_GET_PARAM(TypeTag, Scalar, K);
        phiC_ = EWOMS_GET_PARAM(TypeTag, Scalar, PhiC);
        phi_ = EWOMS_GET_PARAM(TypeTag, Scalar, Phi);
        ciC_ = EWOMS_GET_PARAM(TypeTag, Scalar, CiC);
        cfC_ = EWOMS_GET_PARAM(TypeTag, Scalar, CfC);
        ci_ = EWOMS_GET_PARAM(TypeTag, Scalar, Ci);
        cf_ = EWOMS_GET_PARAM(TypeTag, Scalar, Cf);
        beta_ = EWOMS_GET_PARAM(TypeTag, Scalar, Beta);
        eta_ = EWOMS_GET_PARAM(TypeTag, Scalar, Eta);
        ei_ = EWOMS_GET_PARAM(TypeTag, Scalar, Ei);
        ef_ = EWOMS_GET_PARAM(TypeTag, Scalar, Ef);
        lambda_ = EWOMS_GET_PARAM(TypeTag, Scalar, Lambda);
        llambda_ = EWOMS_GET_PARAM(TypeTag, Scalar, Llambda);
        srw_ = EWOMS_GET_PARAM(TypeTag, Scalar, Srw);
        srn_ = EWOMS_GET_PARAM(TypeTag, Scalar, Srn);
        fineLayerBottom_ = EWOMS_GET_PARAM(TypeTag, Scalar, FineLayerBottom);
        injRate_ = EWOMS_GET_PARAM(TypeTag, Scalar, InjRate);

        // initialize the tables of the fluid system
        // FluidSystem::init();
        FluidSystem::init(/*Tmin=*/temperatureLow_,
                          /*Tmax=*/temperatureHigh_,
                          /*nT=*/nTemperature_,
                          /*pmin=*/pressureLow_,
                          /*pmax=*/pressureHigh_,
                          /*np=*/nPressure_);

        // intrinsic permeabilities
        fineK_ = this->toDimMatrix_(kC_);
        coarseK_ = this->toDimMatrix_(k_);

        // porosities
        finePorosity_ = phiC_;
        coarsePorosity_ = phi_;

        // residual saturations
        fineMaterialParams_.setResidualSaturation(liquidPhaseIdx, srw_);
        fineMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);
        coarseMaterialParams_.setResidualSaturation(liquidPhaseIdx, srw_);
        coarseMaterialParams_.setResidualSaturation(gasPhaseIdx, 0.0);

        // parameters for the DynamicWa saturation functions
        fineMaterialParams_.setEntryPressure(ciC_);
        coarseMaterialParams_.setEntryPressure(ci_);
        fineMaterialParams_.setFinalEntryPressure(cfC_);
        coarseMaterialParams_.setFinalEntryPressure(cf_);
        fineMaterialParams_.setLambda(lambda_);
        coarseMaterialParams_.setLambda(lambda_);
        fineMaterialParams_.setLlambda(llambda_);
        coarseMaterialParams_.setLlambda(llambda_);
        fineMaterialParams_.setBeta(beta_);
        coarseMaterialParams_.setBeta(beta_);
        fineMaterialParams_.setEta(eta_);
        coarseMaterialParams_.setEta(eta_);
        fineMaterialParams_.setEi(ei_);
        coarseMaterialParams_.setEi(ei_);
        fineMaterialParams_.setEf(ef_);
        coarseMaterialParams_.setEf(ef_);

        fineMaterialParams_.finalize();
        coarseMaterialParams_.finalize();

        // parameters for the somerton law thermal conduction
        computeThermalCondParams_(fineThermalCondParams_, finePorosity_);
        computeThermalCondParams_(coarseThermalCondParams_, coarsePorosity_);

        // assume constant heat capacity and granite
        solidEnergyLawParams_.setSolidHeatCapacity(790.0 // specific heat capacity of granite [J / (kg K)]
                                                   * 2700.0); // density of granite [kg/m^3]
        solidEnergyLawParams_.finalize();
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::registerParameters
     */
    static void registerParameters()
    {
        ParentType::registerParameters();

        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureLow,
                             "The lower temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemTemperatureHigh,
                             "The upper temperature [K] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumTemperature,
                             "The number of intervals between the lower and "
                             "upper temperature");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureLow,
                             "The lower pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FluidSystemPressureHigh,
                             "The upper pressure [Pa] for tabulation of the "
                             "fluid system");
        EWOMS_REGISTER_PARAM(TypeTag, unsigned, FluidSystemNumPressure,
                             "The number of intervals between the lower and "
                             "upper pressure");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, KC,
                             "The permeability of the caprock");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, K,
                             "The permeability of the aquifer");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, PhiC,
                             "The porosity of the caprock");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Phi,
                             "The porosity of the aquifer");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, CiC,
                             "The initial entry pressure of the caprock");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, CfC,
                             "The final entry pressure of the caprock");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, S0w,
                             "The initial wetting saturation in the aquifer");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, X0n,
                             "The initial co2 mole fraction in the brine");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, FineLayerBottom,
                             "The height of the aquifer");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, Temperature,
                             "The temperature [K] in the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, MaxDepth,
                             "The maximum depth [m] of the reservoir");
        EWOMS_REGISTER_PARAM(TypeTag, Scalar, InjRate,
                             "The co2 injection rate [kg/(m^3 s)]");
        EWOMS_REGISTER_PARAM(TypeTag, std::string, SimulationName,
                             "The name of the simulation used for the output "
                             "files");
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
        oss << EWOMS_GET_PARAM(TypeTag, std::string, SimulationName)
            << "_" << Model::name();
        if (getPropValue<TypeTag, Properties::EnableEnergy>())
            oss << "_ni";
        oss << "_" << Model::discretizationName();
        return oss.str();
    }

    /*!
     * \copydoc FvBaseProblem::endTimeStep
     */
    void endTimeStep()
    {
#ifndef NDEBUG
        Scalar tol = this->model().newtonMethod().tolerance()*1e5;
        this->model().checkConservativeness(tol);

        // Calculate storage terms
        PrimaryVariables storageL, storageG;
        this->model().globalPhaseStorage(storageL, /*phaseIdx=*/0);
        this->model().globalPhaseStorage(storageG, /*phaseIdx=*/1);

        // Write mass balance information for rank 0
        if (this->gridView().comm().rank() == 0) {
            std::cout << "Storage: liquid=[" << storageL << "]"
                      << " gas=[" << storageG << "]\n" << std::flush;
        }
#endif // NDEBUG
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::temperature
     */
    template <class Context>
    Scalar temperature(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        return temperature_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::intrinsicPermeability
     */
    template <class Context>
    const DimMatrix& intrinsicPermeability(const Context& context, unsigned spaceIdx,
                                           unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::porosity
     */
    template <class Context>
    Scalar porosity(const Context& context, unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return finePorosity_;
        return coarsePorosity_;
    }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::materialLawParams
     */
    template <class Context>
    const MaterialLawParams& materialLawParams(const Context& context,
                                               unsigned spaceIdx, unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineMaterialParams_;
        return coarseMaterialParams_;
    }

    /*!
     * \brief Return the parameters for the heat storage law of the rock
     *
     * In this case, we assume the rock-matrix to be granite.
     */
    template <class Context>
    const SolidEnergyLawParams&
    solidEnergyLawParams(const Context& context OPM_UNUSED,
                         unsigned spaceIdx OPM_UNUSED,
                         unsigned timeIdx OPM_UNUSED) const
    { return solidEnergyLawParams_; }

    /*!
     * \copydoc FvBaseMultiPhaseProblem::thermalConductionParams
     */
    template <class Context>
    const ThermalConductionLawParams &
    thermalConductionLawParams(const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        if (isFineMaterial_(pos))
            return fineThermalCondParams_;
        return coarseThermalCondParams_;
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
     void boundary(BoundaryRateVector& values, const Context& context,
                   unsigned spaceIdx, unsigned timeIdx) const
     {
         const auto& pos = context.pos(spaceIdx, timeIdx);
         if (onRightBoundary_(pos)) {
             Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
             initialFluidState_(fs, context, spaceIdx, timeIdx);
             fs.checkDefined();

             // impose an freeflow boundary condition
             values.setFreeFlow(context, spaceIdx, timeIdx, fs);
         }
         else if (onInlet_(pos)) {
             RateVector massRate(0.0);
             massRate[contiCO2EqIdx] = injRate_;

             using FluidState = Opm::ImmiscibleFluidState<Scalar, FluidSystem>;
             FluidState fs;
             fs.setSaturation(gasPhaseIdx, 1.0);
             const auto& pg =
                 context.intensiveQuantities(spaceIdx, timeIdx).fluidState().pressure(gasPhaseIdx);
             fs.setPressure(gasPhaseIdx, Toolbox::value(pg));
             fs.setTemperature(temperature(context, spaceIdx, timeIdx));

             typename FluidSystem::template ParameterCache<Scalar> paramCache;
             paramCache.updatePhase(fs, gasPhaseIdx);
             Scalar h = FluidSystem::template enthalpy<FluidState, Scalar>(fs, paramCache, gasPhaseIdx);

             // impose an forced inflow boundary condition for pure CO2
             values.setMassRate(massRate);
             values.setEnthalpyRate(massRate[contiCO2EqIdx] * h);
         }
         else
             // no flow on top and bottom
             values.setNoFlow();
     }

    // \}

    /*!
     * \name Volumetric terms
     */
    //! \{

    /*!
     * \copydoc FvBaseProblem::initial
     */
    template <class Context>
    void initial(PrimaryVariables& values, const Context& context, unsigned spaceIdx,
                 unsigned timeIdx) const
    {
        Opm::CompositionalFluidState<Scalar, FluidSystem> fs;
        initialFluidState_(fs, context, spaceIdx, timeIdx);

        // const auto& matParams = this->materialLawParams(context, spaceIdx,
        // timeIdx);
        // values.assignMassConservative(fs, matParams, /*inEquilibrium=*/true);
        values.assignNaive(fs);
    }

    /*!
     * \copydoc FvBaseProblem::source
     *
     * For this problem, the source term of all components is 0
     * everywhere.
     */
    template <class Context>
    void source(RateVector& rate,
                const Context& context OPM_UNUSED,
                unsigned spaceIdx OPM_UNUSED,
                unsigned timeIdx OPM_UNUSED) const
    { rate = Scalar(0.0); }

    //! \}

private:
    template <class Context, class FluidState>
    void initialFluidState_(FluidState& fs,
                            const Context& context,
                            unsigned spaceIdx,
                            unsigned timeIdx) const
    {
        const GlobalPosition& pos = context.pos(spaceIdx, timeIdx);
        unsigned globalIdx = context.globalSpaceIndex(0, /*timeIdx=*/0);
        Scalar densityL = FluidSystem::Brine::liquidDensity(temperature_, Scalar(1e5));
        Scalar depth = maxDepth_ - pos[dim - 1];

        //////
        // set temperature
        //////
        fs.setTemperature(temperature(context, spaceIdx, timeIdx));


        //////set saturations and pressure
        Scalar s0w_ = EWOMS_GET_PARAM(TypeTag, Scalar, S0w);
        if (pos[dim - 1] > fineLayerBottom_){
            fs.setSaturation(FluidSystem::liquidPhaseIdx, 1.0);
            fs.setSaturation(FluidSystem::gasPhaseIdx, 0.0);}
        else {
            fs.setSaturation(FluidSystem::liquidPhaseIdx, s0w_);
            fs.setSaturation(FluidSystem::gasPhaseIdx, 1.0 - s0w_);}

        Scalar pl = 1e5 - densityL * this->gravity()[dim - 1] * depth;

        Scalar pC[numPhases];
        const auto& matParams = this->materialLawParams(context, spaceIdx, timeIdx);
        MaterialLaw::capillaryPressures(pC, matParams, fs);

        fs.setPressure(liquidPhaseIdx, pl + (pC[liquidPhaseIdx] - pC[liquidPhaseIdx]));
        fs.setPressure(gasPhaseIdx, pl + (pC[gasPhaseIdx] - pC[liquidPhaseIdx]));

        //////
        // set composition of the liquid phase
        //////
        Scalar x0n_ = EWOMS_GET_PARAM(TypeTag, Scalar, X0n);
        fs.setMoleFraction(liquidPhaseIdx, CO2Idx, x0n_);
        fs.setMoleFraction(liquidPhaseIdx, BrineIdx,
                           1.0 - fs.moleFraction(liquidPhaseIdx, CO2Idx));

        typename FluidSystem::template ParameterCache<Scalar> paramCache;
        using CFRP = Opm::ComputeFromReferencePhase<Scalar, FluidSystem>;
        CFRP::solve(fs, paramCache,
                    /*refPhaseIdx=*/liquidPhaseIdx,
                    /*setViscosity=*/true,
                    /*setEnthalpy=*/true);
    }

    bool onLeftBoundary_(const GlobalPosition& pos) const
    { return pos[0] < eps_; }

    bool onRightBoundary_(const GlobalPosition& pos) const
    { return pos[0] > this->boundingBoxMax()[0] - eps_; }

    bool onInlet_(const GlobalPosition& pos) const
    { return onLeftBoundary_(pos) && (pos[dim - 1] < fineLayerBottom_); }

    void computeThermalCondParams_(ThermalConductionLawParams& params, Scalar poro)
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

    bool isFineMaterial_(const GlobalPosition& pos) const
    { return pos[dim - 1] > fineLayerBottom_; }

    DimMatrix fineK_;
    DimMatrix coarseK_;
    Scalar fineLayerBottom_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

    ThermalConductionLawParams fineThermalCondParams_;
    ThermalConductionLawParams coarseThermalCondParams_;
    SolidEnergyLawParams solidEnergyLawParams_;

    Scalar kC_;
    Scalar k_;
    Scalar phiC_;
    Scalar phi_;
    Scalar ciC_;
    Scalar cfC_;
    Scalar ci_;
    Scalar cf_;
    Scalar temperature_;
    Scalar maxDepth_;
    Scalar eps_;
    Scalar beta_;
    Scalar eta_;
    Scalar ei_;
    Scalar ef_;
    Scalar lambda_;
    Scalar llambda_;
    Scalar srw_;
    Scalar srn_;
    Scalar injRate_;

    unsigned nTemperature_;
    unsigned nPressure_;

    Scalar pressureLow_, pressureHigh_;
    Scalar temperatureLow_, temperatureHigh_;
};
} // namespace Opm

#endif
