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
* \copydoc Opm::VtkBlackOilMicrobesModule
*/
#ifndef EWOMS_VTK_BLACK_OIL_MICROBES_MODULE_HH
#define EWOMS_VTK_BLACK_OIL_MICROBES_MODULE_HH

#include <opm/material/densead/Math.hpp>

#include "vtkmultiwriter.hh"
#include "baseoutputmodule.hh"

#include <opm/models/utils/propertysystem.hh>
#include <opm/models/utils/parametersystem.hh>
#include <opm/models/blackoil/blackoilproperties.hh>

#include <dune/common/fvector.hh>

#include <cstdio>

namespace Opm::Properties {

namespace TTag {

// create new type tag for the VTK multi-phase output
struct VtkBlackOilMicrobes {};

}  // namespace TTag

// create the property tags needed for the microbes output module
template<class TypeTag, class MyTypeTag>
struct VtkWriteBacteriaConcentration { using type = UndefinedProperty; };

// set default values for what quantities to output
template<class TypeTag>
struct VtkWriteBacteriaConcentration<TypeTag, TTag::VtkBlackOilMicrobes> { static constexpr bool value = true; };

}  // namespace Opm::properties

namespace Opm {
/*!
* \ingroup Vtk
*
* \brief VTK output module for the black-oil microbes model
*/
template <class TypeTag>
class VtkBlackOilMicrobesModule : public BaseOutputModule<TypeTag>
{
    using ParentType = BaseOutputModule<TypeTag>;

    using Simulator = GetPropType<TypeTag, Properties::Simulator>;
    using GridView = GetPropType<TypeTag, Properties::GridView>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Evaluation = GetPropType<TypeTag, Properties::Evaluation>;
    using ElementContext = GetPropType<TypeTag, Properties::ElementContext>;

    static const int vtkFormat = getPropValue<TypeTag, Properties::VtkOutputFormat>();
    using VtkMultiWriter = ::Opm::VtkMultiWriter<GridView, vtkFormat>;

    enum { enableMicrobes = getPropValue<TypeTag, Properties::EnableMicrobes>() };

    using ScalarBuffer = typename ParentType::ScalarBuffer;

public:
    VtkBlackOilMicrobesModule(const Simulator& simulator)
        : ParentType(simulator)
    { }

    /*!
    * \brief Register all run-time parameters for the multi-phase VTK output
    * module.
    */
    static void registerParameters()
    {
        if (!enableMicrobes)
            return;

        Parameters::registerParam<TypeTag, Properties::VtkWriteBacteriaConcentration>
                             ("Include the concentration of microbes in the water phase "
                             "in the VTK output files");
    }

    /*!
    * \brief Allocate memory for the scalar fields we would like to
    *        write to the VTK file.
    */
    void allocBuffers()
    {
        if (!Parameters::get<TypeTag, Parameters::EnableVtkOutput>())
            return;

        if (!enableMicrobes)
            return;

        if (bacteriaConcentrationOutput_())
            this->resizeScalarBuffer_(bacteriaConcentration_);
    }

    /*!
    * \brief Modify the internal buffers according to the intensive quantities relevant for
    *        an element
    */
    void processElement(const ElementContext& elemCtx)
    {
        if (!Parameters::get<TypeTag, Parameters::EnableVtkOutput>())
            return;

        if (!enableMicrobes)
            return;

        for (unsigned dofIdx = 0; dofIdx < elemCtx.numPrimaryDof(/*timeIdx=*/0); ++dofIdx) {
            const auto& intQuants = elemCtx.intensiveQuantities(dofIdx, /*timeIdx=*/0);
            unsigned globalDofIdx = elemCtx.globalSpaceIndex(dofIdx, /*timeIdx=*/0);

            if (bacteriaConcentrationOutput_())
                bacteriaConcentration_[globalDofIdx] =
                    scalarValue(intQuants.bacteriaConcentration());
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

        if (!enableMicrobes)
            return;

        if (bacteriaConcentrationOutput_())
            this->commitScalarBuffer_(baseWriter, "bacteria concentration", bacteriaConcentration_);

    }

private:
    static bool bacteriaConcentrationOutput_()
    {
        static bool val = Parameters::get<TypeTag, Properties::VtkWriteBacteriaConcentration>();
        return val;
    }

    ScalarBuffer bacteriaConcentration_;

};
}  // namespace Opm

#endif