// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*
  Copyright (C) 2011-2013 by Andreas Lauser

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
 * \copydoc Ewoms::VtkMultiPhaseModule
 */
#ifndef EWOMS_VTK_MULTI_PHASE_MODULE_HH
#define EWOMS_VTK_MULTI_PHASE_MODULE_HH

#include <ewoms/vtk/vtkoutputmodule.hh>

#include <opm/core/utility/PropertySystem.hpp>
#include <ewoms/common/parametersystem.hh>

#include <dune/common/fvector.hh>

#include <cstdio>

namespace Opm {
namespace Properties {
// create new type tag for the VTK multi-phase output
NEW_TYPE_TAG(VtkMultiPhase);

// create the property tags needed for the multi phase module
NEW_PROP_TAG(VtkWritePressures);
NEW_PROP_TAG(VtkWriteDensities);
NEW_PROP_TAG(VtkWriteSaturations);
NEW_PROP_TAG(VtkWriteMobilities);
NEW_PROP_TAG(VtkWriteRelativePermeabilities);
NEW_PROP_TAG(VtkWriteViscosities);
NEW_PROP_TAG(VtkWriteAverageMolarMasses);
NEW_PROP_TAG(VtkWritePorosity);
NEW_PROP_TAG(VtkWriteIntrinsicPermeabilities);
NEW_PROP_TAG(VtkWritePotentialGradients);
NEW_PROP_TAG(VtkWriteFilterVelocities);

// set default values for what quantities to output
SET_BOOL_PROP(VtkMultiPhase, VtkWritePressures, true);
SET_BOOL_PROP(VtkMultiPhase, VtkWriteDensities, true);
SET_BOOL_PROP(VtkMultiPhase, VtkWriteSaturations, true);
SET_BOOL_PROP(VtkMultiPhase, VtkWriteMobilities, false);
SET_BOOL_PROP(VtkMultiPhase, VtkWriteRelativePermeabilities, true);
SET_BOOL_PROP(VtkMultiPhase, VtkWriteViscosities, false);
SET_BOOL_PROP(VtkMultiPhase, VtkWriteAverageMolarMasses, false);
SET_BOOL_PROP(VtkMultiPhase, VtkWritePorosity, true);
SET_BOOL_PROP(VtkMultiPhase, VtkWriteIntrinsicPermeabilities, false);
SET_BOOL_PROP(VtkMultiPhase, VtkWritePotentialGradients, false);
SET_BOOL_PROP(VtkMultiPhase, VtkWriteFilterVelocities, false);
} // namespace Properties
} // namespace Opm

namespace Ewoms {
/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for quantities which make sense for all
 *        models which deal with multiple fluid phases in porous media
 *        that don't use flashy concepts like interfacial area.
 *
 * This module deals with the following quantities:
 * - Pressures of all fluid phases
 * - Densities of all fluid phases
 * - Saturations of all fluid phases
 * - Mobilities of all fluid phases
 * - Relative permeabilities of all fluid phases
 * - Viscosities of all fluid phases
 * - Average molar masses of all fluid phases
 * - Porosity of the medium
 * - Norm of the intrinsic permeability of the medium
 */
template<class TypeTag>
class VtkMultiPhaseModule : public VtkOutputModule<TypeTag>
{
    typedef VtkOutputModule<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, DiscVtkOutputModule) DiscVtkOutputModule;

    typedef Ewoms::VtkMultiWriter<GridView> VtkMultiWriter;

    enum { dim = GridView::dimension };
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases) };

    typedef typename ParentType::ScalarBuffer ScalarBuffer;
    typedef typename ParentType::PhaseBuffer PhaseBuffer;

    typedef Dune::FieldVector<Scalar, dim> VelocityVector;
    typedef Dune::BlockVector<VelocityVector> VelocityField;
    typedef std::array<VelocityField, numPhases> PhaseVectorField;

public:
    VtkMultiPhaseModule(const Problem &problem)
        : ParentType(problem)
    {
    }

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePressures, "Include the phase pressures in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteDensities, "Include the phase densities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteSaturations, "Include the phase saturations in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteMobilities, "Include the phase mobilities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteRelativePermeabilities, "Include the phase relative permeabilities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteViscosities, "Include component phase viscosities in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteAverageMolarMasses, "Include the average phase mass in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePorosity, "Include the porosity in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteIntrinsicPermeabilities, "Include the intrinsic permeability in the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteFilterVelocities, "Include in the filter velocities of the phases the VTK output files");
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWritePotentialGradients, "Include the phase pressure potential gradients in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers(VtkMultiWriter &writer)
    {
        if (pressureOutput_()) this->resizePhaseBuffer_(pressure_);
        if (densityOutput_()) this->resizePhaseBuffer_(density_);
        if (saturationOutput_()) this->resizePhaseBuffer_(saturation_);
        if (mobilityOutput_()) this->resizePhaseBuffer_(mobility_);
        if (relativePermeabilityOutput_()) this->resizePhaseBuffer_(relativePermeability_);
        if (viscosityOutput_()) this->resizePhaseBuffer_(viscosity_);
        if (averageMolarMassOutput_()) this->resizePhaseBuffer_(averageMolarMass_);

        if (porosityOutput_()) this->resizeScalarBuffer_(porosity_);
        if (intrinsicPermeabilityOutput_()) this->resizeScalarBuffer_(intrinsicPermeability_);

        if (velocityOutput_()) {
            Scalar nDof = this->problem_.model().numDof();
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                velocity_[phaseIdx].resize(nDof);
                velocity_[phaseIdx] = 0;
            }
            this->resizePhaseBuffer_(velocityWeight_);
        }

        if (potentialGradientOutput_()) {
            Scalar nDof = this->problem_.model().numDof();
            for (int phaseIdx = 0; phaseIdx < numPhases; ++ phaseIdx) {
                potentialGradient_[phaseIdx].resize(nDof);
                potentialGradient_[phaseIdx] = 0;
            }
            this->resizePhaseBuffer_(potentialWeight_);
        }
    }

    /*!
     * \brief Modify the internal buffers according to the volume
     *        variables seen on an element
     */
    void processElement(const ElementContext &elemCtx)
    {
        for (int i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            int I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto &volVars = elemCtx.volVars(i, /*timeIdx=*/0);
            const auto &fs = volVars.fluidState();

            if (porosityOutput_()) porosity_[I] = volVars.porosity();
            if (intrinsicPermeabilityOutput_()) {
                const auto &K = elemCtx.problem().intrinsicPermeability(elemCtx, i, /*timeIdx=*/0);
                intrinsicPermeability_[I] = K[0][0];
            }

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                if (pressureOutput_()) pressure_[phaseIdx][I] = fs.pressure(phaseIdx);
                if (densityOutput_()) density_[phaseIdx][I] = fs.density(phaseIdx);
                if (saturationOutput_()) saturation_[phaseIdx][I] = fs.saturation(phaseIdx);
                if (mobilityOutput_()) mobility_[phaseIdx][I] = volVars.mobility(phaseIdx);
                if (relativePermeabilityOutput_()) relativePermeability_[phaseIdx][I] = volVars.relativePermeability(phaseIdx);
                if (viscosityOutput_()) viscosity_[phaseIdx][I] = fs.viscosity(phaseIdx);
                if (averageMolarMassOutput_()) averageMolarMass_[phaseIdx][I] = fs.averageMolarMass(phaseIdx);
            }
        }

        if (potentialGradientOutput_()) {
            // calculate velocities if requested by the problem
            for (int faceIdx = 0; faceIdx < elemCtx.numInteriorFaces(/*timeIdx=*/0); ++ faceIdx) {
                const auto &fluxVars = elemCtx.fluxVars(faceIdx, /*timeIdx=*/0);

                int i = fluxVars.interiorIndex();
                int I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    Scalar weight = fluxVars.extrusionFactor();

                    Dune::FieldVector<Scalar, dim> pGrad(fluxVars.potentialGrad(phaseIdx));
                    pGrad *= weight;

                    potentialGradient_[phaseIdx][I] += pGrad;
                    potentialWeight_[phaseIdx][I] += weight;
                } // end for all phases
            } // end for all faces
        }

        if (velocityOutput_()) {
            // calculate velocities if requested by the problem
            for (int faceIdx = 0; faceIdx < elemCtx.numInteriorFaces(/*timeIdx=*/0); ++ faceIdx) {
                const auto &fluxVars = elemCtx.fluxVars(faceIdx, /*timeIdx=*/0);

                int i = fluxVars.interiorIndex();
                int I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);

                int j = fluxVars.exteriorIndex();
                int J = elemCtx.globalSpaceIndex(j, /*timeIdx=*/0);

                for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                    Scalar weight = std::max(1e-16, std::abs(fluxVars.volumeFlux(phaseIdx)));
                    Valgrind::CheckDefined(fluxVars.extrusionFactor());
                    assert(fluxVars.extrusionFactor() > 0);
                    weight *= fluxVars.extrusionFactor();

                    Dune::FieldVector<Scalar, dim> v(fluxVars.filterVelocity(phaseIdx));
                    v *= weight;

                    velocity_[phaseIdx][I] += v;
                    velocity_[phaseIdx][J] += v;

                    velocityWeight_[phaseIdx][I] += weight;
                    velocityWeight_[phaseIdx][J] += weight;
                } // end for all phases
            } // end for all faces
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(VtkMultiWriter &writer)
    {
        if (pressureOutput_()) this->commitPhaseBuffer_(writer, "pressure_%s", pressure_);
        if (densityOutput_()) this->commitPhaseBuffer_(writer, "density_%s", density_);
        if (saturationOutput_()) this->commitPhaseBuffer_(writer, "saturation_%s", saturation_);
        if (mobilityOutput_()) this->commitPhaseBuffer_(writer, "mobility_%s", mobility_);
        if (relativePermeabilityOutput_()) this->commitPhaseBuffer_(writer, "relativePerm_%s", relativePermeability_);
        if (viscosityOutput_()) this->commitPhaseBuffer_(writer, "viscosity_%s", viscosity_);
        if (averageMolarMassOutput_()) this->commitPhaseBuffer_(writer, "averageMolarMass_%s", averageMolarMass_);

        if (porosityOutput_()) this->commitScalarBuffer_(writer, "porosity", porosity_);
        if (intrinsicPermeabilityOutput_()) this->commitScalarBuffer_(writer, "intrinsicPerm", intrinsicPermeability_);

        if (velocityOutput_()) {
            int nDof = this->problem_.model().numDof();

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (int i = 0; i < nDof; ++i)
                    velocity_[phaseIdx][i] /= velocityWeight_[phaseIdx][i];
                // commit the phase velocity
                char name[512];
                snprintf(name, 512, "filterVelocity_%s", FluidSystem::phaseName(phaseIdx));

                DiscVtkOutputModule::attachDofData_(writer, velocity_[phaseIdx], name, dim);
            }
        }

        if (potentialGradientOutput_()) {
            int nDof = this->problem_.model().numDof();

            for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
                // first, divide the velocity field by the
                // respective finite volume's surface area
                for (int i = 0; i < nDof; ++i)
                    potentialGradient_[phaseIdx][i] /= potentialWeight_[phaseIdx][i];
                // commit the phase velocity
                char name[512];
                snprintf(name, 512, "gradP_%s", FluidSystem::phaseName(phaseIdx));

                DiscVtkOutputModule::attachDofData_(writer, potentialGradient_[phaseIdx], name, dim);
            }
        }
    }

private:
    static bool pressureOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWritePressures); }

    static bool densityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteDensities); }

    static bool saturationOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteSaturations); }

    static bool mobilityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteMobilities); }

    static bool relativePermeabilityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteRelativePermeabilities); }

    static bool viscosityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteViscosities); }

    static bool averageMolarMassOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteAverageMolarMasses); }

    static bool porosityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWritePorosity); }

    static bool intrinsicPermeabilityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteIntrinsicPermeabilities); }

    static bool velocityOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWriteFilterVelocities); }

    static bool potentialGradientOutput_()
    { return EWOMS_GET_PARAM(TypeTag, bool, VtkWritePotentialGradients); }

    PhaseBuffer pressure_;
    PhaseBuffer density_;
    PhaseBuffer saturation_;
    PhaseBuffer mobility_;
    PhaseBuffer relativePermeability_;
    PhaseBuffer viscosity_;
    PhaseBuffer averageMolarMass_;

    ScalarBuffer porosity_;
    ScalarBuffer intrinsicPermeability_;

    PhaseVectorField velocity_;
    PhaseBuffer velocityWeight_;

    PhaseVectorField potentialGradient_;
    PhaseBuffer potentialWeight_;
};

} // namespace Ewoms

#endif