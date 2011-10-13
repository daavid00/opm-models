// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
 *   Copyright (C) 2008-2009 by Klaus Mosthaf                                *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the non-isothermal two-phase box model.
 *
 */
#ifndef DUMUX_2PNI_LOCAL_RESIDUAL_HH
#define DUMUX_2PNI_LOCAL_RESIDUAL_HH

#include <dumux/boxmodels/2p/2plocalresidual.hh>

#include <dumux/boxmodels/2pni/2pnivolumevariables.hh>
#include <dumux/boxmodels/2pni/2pnifluxvariables.hh>

#include <dumux/boxmodels/2pni/2pniproperties.hh>
namespace Dumux
{
/*!
 * \ingroup TwoPNIModel
 * \ingroup BoxLocalResidual
 * \brief Element-wise calculation of the Jacobian matrix for problems
 *        using the non-isothermal two-phase box model.
 */
template<class TypeTag>
class TwoPNILocalResidual : public TwoPLocalResidual<TypeTag>
{
    typedef TwoPLocalResidual<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, VolumeVariables) VolumeVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluxVariables) FluxVariables;

    typedef typename GET_PROP_TYPE(TypeTag, TwoPIndices) Indices;
    enum {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        energyEqIdx = Indices::energyEqIdx,
        temperatureIdx = Indices::temperatureIdx,
    };

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    enum { dimWorld = GridView::dimensionworld };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, ElementContext) ElementContext;
    typedef Dune::FieldVector<Scalar, dimWorld> Vector;

public:
    /*!
     * \brief Evaluate the amount all conservation quantities
     *        (e.g. phase mass) within a sub-control volume.
     *
     * The result should be averaged over the volume (e.g. phase mass
     * inside a sub control volume divided by the volume)
     *
     *  \param result The storage of the conservation quantitiy (mass or energy) within the sub-control volume
     *  \param scvIdx The SCV (sub-control-volume) index
     *  \param usePrevSol Evaluate function with solution of current or previous time step
     */
    void computeStorage(PrimaryVariables &storage,
                        const ElementContext &elemCtx,
                        int scvIdx,
                        int historyIdx) const
    {
        // compute the storage term for phase mass
        ParentType::computeStorage(storage, elemCtx, scvIdx, historyIdx);

        const VolumeVariables &volVars =
            elemCtx.volVars(scvIdx, historyIdx);

        storage[energyEqIdx] = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx)
        {
            // add the internal energy of the phase
            storage[energyEqIdx] +=
                volVars.porosity() * (
                    volVars.fluidState().density(phaseIdx)
                    * volVars.fluidState().internalEnergy(phaseIdx)
                    * volVars.fluidState().saturation(phaseIdx));
        };

        // handle the heat capacity of the solid
        storage[energyEqIdx] +=
            volVars.fluidState().temperature(/*phaseIdx=*/0)
            * volVars.heatCapacitySolid()
            * (1 - volVars.porosity());
    }

    /*!
     * \brief Evaluates the advective mass flux and the heat flux
     * over a face of a subcontrol volume and writes the result in
     * the flux vector.
     *
     * \param flux The advective flux over the SCV (sub-control-volume) face for each component
     * \param fluxData The flux variables at the current SCV face
     *
     * This method is called by compute flux (base class)
     */
    void computeAdvectiveFlux(PrimaryVariables &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx) const
    {
        // advective mass flux
        ParentType::computeAdvectiveFlux(flux, elemCtx, scvfIdx);

        const auto &fluxVars = elemCtx.fluxVars(scvfIdx);
        const auto &evalPointFluxVars = elemCtx.evalPointFluxVars(scvfIdx);
        
        // advective heat flux in all phases
        flux[energyEqIdx] = 0;
        for (int phaseIdx = 0; phaseIdx < numPhases; ++phaseIdx) {
            // vertex data of the upstream and the downstream vertices
            const VolumeVariables &up = elemCtx.volVars(evalPointFluxVars.upstreamIdx(phaseIdx));
            const VolumeVariables &dn = elemCtx.volVars(evalPointFluxVars.downstreamIdx(phaseIdx));

            flux[energyEqIdx] +=
                fluxVars.filterVelocityNormal(phaseIdx)
                * (fluxVars.upstreamWeight(phaseIdx)
                   * up.fluidState().enthalpy(phaseIdx)
                   * up.fluidState().density(phaseIdx)
                   +
                   fluxVars.downstreamWeight(phaseIdx)
                   * dn.fluidState().enthalpy(phaseIdx)
                   * dn.fluidState().density(phaseIdx));
        }
    }

    /*!
     * \brief Adds the diffusive heat flux to the flux vector over
     *        the face of a sub-control volume.
     *
     * \param flux The diffusive flux over the SCV (sub-control-volume) face for each conservation quantity (mass, energy)
     * \param fluxData The flux variables at the current SCV face
     *
     * This method is called by compute flux (base class)
     */
    void computeDiffusiveFlux(PrimaryVariables &flux,
                              const ElementContext &elemCtx,
                              int scvfIdx) const
    {
        // diffusive mass flux
        ParentType::computeDiffusiveFlux(flux, elemCtx, scvfIdx);

        const auto &fluxVars = elemCtx.fluxVars(scvfIdx);

        // diffusive heat flux
        flux[energyEqIdx] +=
            -fluxVars.temperatureGradNormal()
            * fluxVars.heatConductivity();
    }
};

}

#endif
