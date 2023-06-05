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
* \brief Contains the classes required to extend the black-oil model with (simple) microbe interaction.
*/
#ifndef EWOMS_BLACK_OIL_MICROBES_MODULE_HH
#define EWOMS_BLACK_OIL_MICROBES_MODULE_HH

#include "blackoilproperties.hh"

#include <opm/models/blackoil/blackoilmicrobesparams.hh>
#include <opm/models/io/vtkblackoilmicrobesmodule.hh>

#if HAVE_ECL_INPUT
#include <opm/input/eclipse/EclipseState/EclipseState.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PermporoTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/PefactTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/TableManager.hpp>
#include <opm/input/eclipse/EclipseState/Tables/SimpleTable.hpp>
#include <opm/input/eclipse/EclipseState/Tables/Bactpara.hpp>
#endif

#include <dune/common/fvector.hh>

#include <cmath>
#include <stdexcept>
#include <string>

namespace Opm {
/*!
* \ingroup BlackOil
* \brief Contains the high level supplements required to extend the black oil
*        model with (simple) microbe interaction.
*/
template <class TypeTag, bool enableMicrobesV = getPropValue<TypeTag, Properties::EnableMicrobes>()>
class BlackOilMicrobesModule
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using IntensiveQuantities = GetPropType<TypeTag, Properties::IntensiveQuantities>;
    using ExtensiveQuantities = GetPropType<TypeTag, Properties::ExtensiveQuantities>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Model = GetPropType<TypeTag, Properties::Model>;
    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using EqVector = GetPropType<TypeTag, Properties::EqVector>;
    using RateVector = GetPropType<TypeTag, Properties::RateVector>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;

    using Toolbox = MathToolbox<Evaluation>;

    using TabulatedFunction = typename BlackOilMicrobesParams<Scalar>::TabulatedFunction;

    enum { gasCompIdx = FluidSystem::gasCompIdx };

    static constexpr unsigned conti0EqIdx = Indices::conti0EqIdx;
    static constexpr unsigned bacteriaConcentrationIdx = Indices::bacteriaConcentrationIdx;
    static constexpr unsigned contiBacteriaEqIdx = Indices::contiBacteriaEqIdx;
    static constexpr unsigned gasPhaseIdx = FluidSystem::gasPhaseIdx;
    static constexpr unsigned waterPhaseIdx = FluidSystem::waterPhaseIdx;

    static constexpr unsigned enableMicrobes = enableMicrobesV;

    static constexpr unsigned numEq = getPropValue<TypeTag, Properties::NumEq>();

public:

#if HAVE_ECL_INPUT
    /*!
    * \brief Initialize internal data needed for microbes module
    */
    static void initFromState(const EclipseState& eclState)
    {
        // Check if module is enabled but MICROBES is not present
        if (enableMicrobes && !eclState.runspec().microbes()) {
            throw std::runtime_error("Microbes module enabled at compile time, but the deck does not contain "
                                     "MICROBES!");
        }
        // Check for opposite of the above: module disabled but MICROBES is in deck
        else if (!enableMicrobes && eclState.runspec().microbes()) {
            throw std::runtime_error("Microbes module disabled at compile time, but deck contains MICROBES!");
        }
        
        // If MICROBES is not present, this module should not be active regardless
        if (!eclState.runspec().microbes())
            return;

        const auto& tableManager = eclState.getTableManager();
        unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
        const TableContainer& permporoTables = tableManager.getPermporoTables();
        params_.permporoTable_.resize(numSatRegions);
        for (size_t i = 0; i < permporoTables.size(); ++i) {
            const PermporoTable& permporoTable = permporoTables.getTable<PermporoTable>(i);
            params_.permporoTable_[i].setXYContainers(permporoTable.getPorosityChangeColumn(), permporoTable.getPermeabilityMultiplierColumn());
        }
        const TableContainer& pefactTables = tableManager.getPefactTables();
        if (!pefactTables.empty()) {
            params_.pefactTable_.resize(numSatRegions);
            for (size_t i = 0; i < pefactTables.size(); ++i) {
                const PefactTable& pefactTable = pefactTables.getTable<PefactTable>(i);
                params_.pefactTable_[i].setXYContainers(pefactTable.getPorosityChangeColumn(), pefactTable.getPcMultiplierColumn());
            }
        }
        const auto& bactpara = tableManager.getBactpara();
        if (!bactpara.empty()) {
            unsigned numSatRegions = tableManager.getTabdims().getNumSatTables();
            params_.biofilmDensity_.resize(numSatRegions);
            params_.maxGrowthRate_.resize(numSatRegions);
            params_.halfVelocityCoeff_.resize(numSatRegions);
            params_.yieldCoeff_.resize(numSatRegions);
            params_.decayCoeff_.resize(numSatRegions);
            for (size_t i = 0; i < bactpara.size(); ++i) {
                params_.biofilmDensity_[i] = bactpara[i].biofilm_density;
                params_.maxGrowthRate_[i] = bactpara[i].max_growth_rate;
                params_.halfVelocityCoeff_[i] = bactpara[i].half_velocity_coefficient;
                params_.yieldCoeff_[i] = bactpara[i].yield_coefficient;
                params_.decayCoeff_[i] = bactpara[i].decay_coefficient;
            }
        }
        else {
            throw std::runtime_error("BACTPARA must be specified in MICROBES runs\n");
        }
    }
#endif

    /*!
    * \brief Register all run-time parameters for the black-oil microbes module.
    */
    static void registerParameters()
    {
        if (!enableMicrobes)
            return;

        VtkBlackOilMicrobesModule<TypeTag>::registerParameters();
    }

    /*!
    * \brief Register all microbes specific VTK and ECL output modules.
    */
    static void registerOutputModules(Model& model,
                                      Simulator& simulator)
    {
        if (!enableMicrobes)
            return;

        model.addOutputModule(new VtkBlackOilMicrobesModule<TypeTag>(simulator));
    }

    static bool eqApplies(unsigned eqIdx)
    {
        if (!enableMicrobes)
            return false;
        
        // True if bacteria equation applies
        return eqIdx == contiBacteriaEqIdx;
    }

    static Scalar eqWeight([[maybe_unused]] unsigned eqIdx)
    {
        assert(eqApplies(eqIdx));

        // Should be OK if bacteria concentration if in kg/m3 or similar.
        return static_cast<Scalar>(1.0);
    }

    template <class LhsEval>
    static void addStorage(Dune::FieldVector<LhsEval, numEq>& storage,
                           const IntensiveQuantities& intQuants)
    {
        if (!enableMicrobes)
            return;

        const LhsEval solidBiofilm =
                Toolbox::template decay<LhsEval>(intQuants.referencePorosity())
                * Toolbox::template decay<LhsEval>(intQuants.bacteriaConcentration());
        storage[contiBacteriaEqIdx] += solidBiofilm*1e-6;
        
        // // Calculate volume water in SM3
        // const auto& fs = intQuants.fluidState();
        // LhsEval surfaceVolumeWater =
        //     Toolbox::template decay<LhsEval>(fs.saturation(waterPhaseIdx))
        //     * Toolbox::template decay<LhsEval>(fs.invB(waterPhaseIdx))
        //     * Toolbox::template decay<LhsEval>(intQuants.porosity());

        // // Avoid singular matrix if no water is present.
        // surfaceVolumeWater = max(surfaceVolumeWater, 1e-10);

        // // Microbes suspended in water phase
        // const LhsEval massMicrobes = 
        //     surfaceVolumeWater * Toolbox::template decay<LhsEval>(intQuants.bacteriaConcentration());
        // storage[contiBacteriaEqIdx] += massMicrobes;
    }

    static void computeFlux(RateVector& flux,
                            const ElementContext& elemCtx,
                            unsigned scvfIdx,
                            unsigned timeIdx)
    {
        if (!enableMicrobes)
            return;
        
        // const auto& extQuants = elemCtx.extensiveQuantities(scvfIdx, timeIdx);

        // const unsigned upIdx = extQuants.upstreamIndex(waterPhaseIdx);
        // const unsigned inIdx = extQuants.interiorIndex();
        // const auto& up = elemCtx.intensiveQuantities(upIdx, timeIdx);


        // if (upIdx == inIdx) {
        //     flux[contiBacteriaEqIdx] = extQuants.volumeFlux(waterPhaseIdx) * up.bacteriaConcentration();
        // }
        // else {
        //     flux[contiBacteriaEqIdx] = extQuants.volumeFlux(waterPhaseIdx) * decay<Scalar>(up.bacteriaConcentration());
        // }
    }

    static void addSource(RateVector& source,
                            const ElementContext& elemCtx,
                            unsigned dofIdx,
                            unsigned timeIdx)
    {
        if (!enableMicrobes)
            return;
        // Get bacteria parameters
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, dofIdx, timeIdx);
        Scalar rhob = biofilmDensity(satnumRegionIdx);
        Scalar mu = maxGrowthRate(satnumRegionIdx);
        Scalar Kn = halfVelocityCoeff(satnumRegionIdx);
        Scalar Y = yieldCoeff(satnumRegionIdx);
        Scalar kd = decayCoeff(satnumRegionIdx);

        // Convert Rsw to concentration to use in source term
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        const auto& fs = intQuants.fluidState();
        const auto& Rsw = fs.Rsw();
        const auto& rhow = fs.density(waterPhaseIdx);
        unsigned pvtRegionIndex = fs.pvtRegionIndex();

        const auto& xG = RswToMassFraction(pvtRegionIndex, Rsw);

        // Get saturation, porosity, bacteria concentration and inverse Bg for convenience
        const Evaluation& poroRef = intQuants.referencePorosity();
        const Evaluation& poro = intQuants.porosity();
        const Evaluation& sw = fs.saturation(waterPhaseIdx);
        const Evaluation& cBact = intQuants.bacteriaConcentration() * poroRef;
        Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, pvtRegionIndex);

        // Calculate bacteria growth rate
        Evaluation kg = mu * (xG * rhow / (xG * rhow + Kn));

        // Compute source terms
        // Microbial growth and decay rate
        source[contiBacteriaEqIdx] += (kg - kd) * cBact * 1e-6;

        // Microbial consumption of dissolved gas is proportional to bacterial growth rate
        unsigned activeGasCompIdx = Indices::canonicalToActiveComponentIndex(gasCompIdx);
        source[activeGasCompIdx] -= cBact * rhob * kg / (Y * rho_gRef);
    }

    static const TabulatedFunction& permporoTable(const ElementContext& elemCtx,
                                                  unsigned scvIdx,
                                                  unsigned timeIdx)
    {
        unsigned satnumRegionIdx = elemCtx.problem().satnumRegionIndex(elemCtx, scvIdx, timeIdx);
        return params_.permporoTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& permporoTable(unsigned satnumRegionIdx)
    {
        return params_.permporoTable_[satnumRegionIdx];
    }

    static const TabulatedFunction& pefactTable(unsigned satnumRegionIdx)
    {
        return params_.pefactTable_[satnumRegionIdx];
    }

    static const Scalar biofilmDensity(unsigned satnumRegionIdx)
    {
        return params_.biofilmDensity_[satnumRegionIdx];
    }

    static const Scalar maxGrowthRate(unsigned satnumRegionIdx)
    {
        return params_.maxGrowthRate_[satnumRegionIdx];
    }

    static const Scalar halfVelocityCoeff(unsigned satnumRegionIdx)
    {
        return params_.halfVelocityCoeff_[satnumRegionIdx];
    }

    static const Scalar yieldCoeff(unsigned satnumRegionIdx)
    {
        return params_.yieldCoeff_[satnumRegionIdx];
    }

    static const Scalar decayCoeff(unsigned satnumRegionIdx)
    {
        return params_.decayCoeff_[satnumRegionIdx];
    }

    static bool hasPefactTables()
    {
        if constexpr (enableMicrobes)
            return !params_.pefactTable_.empty();
        else
            return false;
    }

private:
    static BlackOilMicrobesParams<Scalar> params_;

    static Evaluation RswToMassFraction(unsigned regionIdx, const Evaluation& Rsw) {
        Scalar rho_wRef = FluidSystem::referenceDensity(FluidSystem::waterPhaseIdx, regionIdx);
        Scalar rho_gRef = FluidSystem::referenceDensity(FluidSystem::gasPhaseIdx, regionIdx);

        const Evaluation rho_oG = Rsw * rho_gRef;

        return rho_oG/(rho_wRef + rho_oG);
    }

};  // class BlackOilMicrobesModule

template <class TypeTag, bool enableMicrobesV>
BlackOilMicrobesParams< typename BlackOilMicrobesModule<TypeTag, enableMicrobesV>::Scalar >
BlackOilMicrobesModule<TypeTag, enableMicrobesV>::params_;

/*!
* \ingroup BlackOil
* \class Opm::BlackOilMicrobesIntensiveQuantities
*
* \brief Provides the volumetric quantities required for the equations needed by the
*        (simple) microbes extension of the black-oil model.
*/
template <class TypeTag, bool enableMicrobesV = getPropValue<TypeTag, Properties::EnableMicrobes>()>
class BlackOilMicrobesIntensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::IntensiveQuantities>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using Indices = GetPropType<TypeTag, Properties::Indices>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using MicrobesModule = BlackOilMicrobesModule<TypeTag>;

    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };
    static constexpr int bacteriaConcentrationIdx = Indices::bacteriaConcentrationIdx;

public:
    /*!
    * \brief Update intensive quantities associated with module for (simple) microbial interaction
    */

    void MicrobesPropertiesUpdate_(const ElementContext& elemCtx,
                                   unsigned dofIdx,
                                   unsigned timeIdx)
    {
        const auto linearizationType = elemCtx.linearizationType();
        const PrimaryVariables& priVars = elemCtx.primaryVars(dofIdx, timeIdx);

        // update bacteria concentration from primary variables
        bacteriaConcentration_ = priVars.makeEvaluation(bacteriaConcentrationIdx, timeIdx, linearizationType);
        const Evaluation porosityFactor  = min(1.0 - bacteriaConcentration_, 1.0); //phi/phi_0
        unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, timeIdx);
        const auto totVolume = elemCtx.simulator().model().dofTotalVolume(globalDofIdx);
        const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, timeIdx);
        biofilmMass_ = bacteriaConcentration_ * elemCtx.problem().biofilmDensity(dofIdx) * intQuants.referencePorosity();

        const auto& permporoTable = MicrobesModule::permporoTable(elemCtx, dofIdx, timeIdx);

        permPoro_ = permporoTable.eval(porosityFactor);
        for (unsigned phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            if (!FluidSystem::phaseIsActive(phaseIdx))
                continue;

            asImp_().mobility_[phaseIdx] *= permPoro_;
        }
    }

    const Evaluation& bacteriaConcentration() const
    { return bacteriaConcentration_; }

    const Evaluation& biofilmMass() const
    { return biofilmMass_; }

    const Evaluation& permPoro() const
    { return permPoro_; }

protected:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

    Evaluation bacteriaConcentration_;
    Evaluation biofilmMass_;
    Evaluation permPoro_;

};  // class BlackOilMicrobesIntensiveQuantities

template <class TypeTag>
class BlackOilMicrobesIntensiveQuantities<TypeTag, false>
{
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

public:
    void MicrobesPropertiesUpdate_(const ElementContext& elemCtx,
                                   unsigned dofIdx,
                                   unsigned timeIdx)
    {}
    
    const Evaluation& bacteriaConcentration() const
    { throw std::logic_error("bacteriaConcentration() called but MICROBES is disabled!"); }

    const Evaluation& biofilmMass() const
    { throw std::logic_error("biofilmMass() called but MICROBES is disabled!"); }

     const Evaluation& permPoro() const
    { throw std::logic_error("permPoro() called but MICROBES is disabled"); }

};  // class BlackOilMicrobesIntensiveQuantities<TypeTag, false>

/*!
* \ingroup BlackOil
* \class Opm::BlackOilMicrobesExtensiveQuantities
*
* \brief Extensive quantities for (simple) microbes interaction in black-oil model
*/
template <class TypeTag, bool enableMicrobesV = getPropValue<TypeTag, Properties::EnableMicrobes>()>
class BlackOilMicrobesExtensiveQuantities
{
    using Implementation = GetPropType<TypeTag, Properties::ExtensiveQuantities>;

private:
    Implementation& asImp_()
    { return *static_cast<Implementation*>(this); }

};

template <class TypeTag>
class BlackOilMicrobesExtensiveQuantities<TypeTag, false>{};

}  // namespace Opm

#endif