// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2012 by Andreas Lauser                                    *
 *   Copyright (C) 2012 by Bernd Flemisch                                    *
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
 * \brief A Ncp specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_NCP_NEWTON_CONTROLLER_HH
#define DUMUX_NCP_NEWTON_CONTROLLER_HH

#include "ncpproperties.hh"

#include <dumux/boxmodels/common/boxnewtoncontroller.hh>

#include <algorithm>

namespace Dumux {

template <class TypeTag>
class NcpNewtonChop
{
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Model) Model;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;

    enum { numPhases =  GET_PROP_VALUE(TypeTag, NumPhases) };
    enum { numComponents =  GET_PROP_VALUE(TypeTag, NumComponents) };
    enum { fugacity0Idx = Indices::fugacity0Idx };
    enum { saturation0Idx = Indices::saturation0Idx };
    enum { pressure0Idx = Indices::pressure0Idx };

public:
    static void chop(SolutionVector &uCurrentIter,
                     const SolutionVector &uLastIter, 
                     const Model &model)
    {
        for (size_t i = 0; i < uLastIter.size(); ++i) {
            for (unsigned phaseIdx = 0; phaseIdx < numPhases - 1; ++phaseIdx)
                saturationChop_(uCurrentIter[i][saturation0Idx + phaseIdx],
                                uLastIter[i][saturation0Idx + phaseIdx],
                                model);
            pressureChop_(uCurrentIter[i][pressure0Idx], 
                          uLastIter[i][pressure0Idx],
                          model);
            
            // fugacities
            for (int compIdx = 0; compIdx < numComponents; ++compIdx) {
                Scalar &val = uCurrentIter[i][fugacity0Idx + compIdx];
                Scalar oldVal = uLastIter[i][fugacity0Idx + compIdx];

                // allow the mole fraction of the component to change
                // at most 70% (assuming composition independent
                // fugacity coefficients)
                Scalar minPhi = model.minActivityCoeff(i, compIdx);
                Scalar maxDelta = 0.7 * minPhi;

                clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);

                // do not allow mole fractions lager than 101% or
                // smaller than -1%
                val = std::max(-0.01*minPhi, val);
                val = std::min(1.01*minPhi, val);

            }
        }
    }

private:
    static void clampValue_(Scalar &val, Scalar minVal, Scalar maxVal)
    {
        val = std::max(minVal, std::min(val, maxVal));
    }

    static void pressureChop_(Scalar &val, Scalar oldVal, const Model &model)
    {
        // limit pressure updates to 20% per iteration
        clampValue_(val, oldVal*0.8, oldVal*1.2);
    }

    static void saturationChop_(Scalar &val, Scalar oldVal, const Model &model)
    {
        // limit saturation updates to 20% per iteration
        const Scalar maxDelta = 0.20;
        clampValue_(val, oldVal - maxDelta, oldVal + maxDelta);
    }

};

/*!
 * \ingroup Newton
 * \brief A Ncp specific controller for the newton solver.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class NcpNewtonController : public BoxNewtonController<TypeTag>
{
    typedef BoxNewtonController<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GlobalEqVector) GlobalEqVector;
    typedef typename GET_PROP_TYPE(TypeTag, SolutionVector) SolutionVector;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    enum {
        numEq = GET_PROP_VALUE(TypeTag, NumEq),
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag, NumComponents)
    };

    typedef NcpNewtonChop<TypeTag> NewtonChop;

public:
    NcpNewtonController(Problem &problem)
        : ParentType(problem)
    {
        choppedIterations_ = GET_PARAM_FROM_GROUP(TypeTag, int, Newton, ChoppedIterations);
        Dune::FMatrixPrecision<Scalar>::set_singular_limit(1e-35);
    }

    void newtonUpdate(SolutionVector &uCurrentIter,
                      const SolutionVector &uLastIter,
                      const GlobalEqVector &deltaU)
    {
        // make sure not to swallow non-finite values at this point
        if (!std::isfinite(deltaU.two_norm2()))
            DUNE_THROW(NumericalProblem, "Non-finite update!\n");

        // compute the vertex and element colors for partial
        // reassembly
        if (this->enablePartialReassemble_) {
            const Scalar minReasmTol = 1e-2*this->tolerance_;
            const Scalar maxReasmTol = 1e1*this->tolerance_;
            Scalar reassembleTol = std::max(minReasmTol, std::min(maxReasmTol, this->error_/1e4));

            this->model_().jacobianAssembler().updateDiscrepancy(uLastIter, deltaU);
            this->model_().jacobianAssembler().computeColors(reassembleTol);
        }

        if (this->useLineSearch_)
            this->lineSearchUpdate_(uCurrentIter, uLastIter, deltaU);
        else {
            for (size_t i = 0; i < uLastIter.size(); ++i) {
                for (int j = 0; j < numEq; ++j) {
                    uCurrentIter[i][j] = uLastIter[i][j] - deltaU[i][j];
                }
            }
           
            if (this->numSteps_ < choppedIterations_) {
                // put crash barriers along the update path at the
                // beginning...
                NewtonChop::chop(uCurrentIter, uLastIter, this->model_());
            }
        }
    }

protected:
    int choppedIterations_;
};
}

#endif
