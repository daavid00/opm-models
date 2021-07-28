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
 * \copydoc Opm::VtkWaModule
 */
#ifndef EWOMS_VTK_WA_HH
#define EWOMS_VTK_WA_HH

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>

#include <opm/material/common/MathToolbox.hpp>
#include <opm/material/common/Valgrind.hpp>

#include <dune/common/fvector.hh>

#include <cstdio>

namespace Opm::Properties {

namespace TTag {

// create new type tag for the VTK wa output
struct VtkWa {};

} // namespace TTag

// create the property tags needed for the wa
template<class TypeTag, class MyTypeTag>
struct VtkWriteWa { using type = UndefinedProperty; };

// set default values for what quantities to output
template<class TypeTag>
struct VtkWriteWa<TypeTag, TTag::VtkWa> { static constexpr bool value = true; };

} // namespace Opm::Properties

namespace Opm {

/*!
 * \ingroup Vtk
 *
 * \brief VTK output module for wettability alteration
 *
 * This module deals with the following quantities:
 * - Wettability of the medium
 * - Norm of the intrinsic permeability of the medium
 */
template<class TypeTag>
class VtkWaModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using DiscBaseOutputModule = GetPropType<TypeTag, Properties::DiscBaseOutputModule>;

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { dimWorld = GridView::dimensionworld };
    enum { numPhases = getPropValue<TypeTag, Properties::NumPhases>() };

    using ScalarBuffer = typename ParentType::ScalarBuffer;
    using VectorBuffer = typename ParentType::VectorBuffer;
    using TensorBuffer = typename ParentType::TensorBuffer;
    using PhaseBuffer = typename ParentType::PhaseBuffer;

    using DimVector = Dune::FieldVector<Scalar, dimWorld>;

    using PhaseVectorBuffer = std::array<VectorBuffer, numPhases>;

public:
    VtkWaModule(const Simulator& simulator)
        : ParentType(simulator)
    {}

    /*!
     * \brief Register all run-time parameters for the multi-phase VTK output module.
     */
    static void registerParameters()
    {
        EWOMS_REGISTER_PARAM(TypeTag, bool, VtkWriteWa,
                             "Include the wettability alteration in the VTK output files");
    }

    /*!
     * \brief Allocate memory for the scalar fields we would like to
     *        write to the VTK file.
     */
    void allocBuffers()
    {
        if (waOutput_()) this->resizeScalarBuffer_(wa_);
    }

    /*!
     * \brief Modify the internal buffers according to the intensive quantities seen on
     *        an element
     */
    void processElement(const ElementContext& elemCtx)
    {
        if (!EWOMS_GET_PARAM(TypeTag, bool, EnableVtkOutput))
            return;

        const auto& problem = elemCtx.problem();
        for (unsigned i = 0; i < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++i) {
            unsigned I = elemCtx.globalSpaceIndex(i, /*timeIdx=*/0);
            const auto& intQuants = elemCtx.intensiveQuantities(i, /*timeIdx=*/0);
            const auto& fs = intQuants.fluidState();
            if (waOutput_()) wa_[I] = Opm::getValue(intQuants.wa());
        }
    }

    /*!
     * \brief Add all buffers to the VTK output writer.
     */
    void commitBuffers(BaseOutputWriter& baseWriter)
    {
        VtkMultiWriter *vtkWriter = dynamic_cast<VtkMultiWriter*>(&baseWriter);
        if (!vtkWriter)
            return;
        if (waOutput_())
            this->commitScalarBuffer_(baseWriter, "wa", wa_);
    }

private:
    static bool waOutput_()
    {
        static bool val = EWOMS_GET_PARAM(TypeTag, bool, VtkWriteWa);
        return val;
    }

    ScalarBuffer wa_;
};

} // namespace Opm

#endif
