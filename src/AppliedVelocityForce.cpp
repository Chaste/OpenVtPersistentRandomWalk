/*

Copyright (c) 2005-2024, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#include "AppliedVelocityForce.hpp"
#include "AbstractCentreBasedCellPopulation.hpp"

template<unsigned DIM>
AppliedVelocityForce<DIM>::AppliedVelocityForce()
    : AbstractForce<DIM>(),
	  mMovementParameter(0.01)
{
}

template<unsigned DIM>
AppliedVelocityForce<DIM>::~AppliedVelocityForce()
{
}

template<unsigned DIM>
void AppliedVelocityForce<DIM>::SetMovementParameter(double movementParameter)
{
    assert(movementParameter > 0.0);
    mMovementParameter = movementParameter;
}

template<unsigned DIM>
double AppliedVelocityForce<DIM>::GetMovementParameter()
{
    return mMovementParameter;
}

template<unsigned DIM>
void AppliedVelocityForce<DIM>::AddForceContribution(AbstractCellPopulation<DIM>& rCellPopulation)
{
    // Throw an exception message if not using a subclass of AbstractCentreBasedCellPopulation
    if (dynamic_cast<AbstractCentreBasedCellPopulation<DIM,DIM>*>(&rCellPopulation) == nullptr)
    {
        EXCEPTION("AppliedVelocityForce is to be used with subclasses of AbstractCentreBasedCellPopulation only");
    }

    double dt = SimulationTime::Instance()->GetTimeStep();

    // Iterate over the nodes
    for (typename AbstractMesh<DIM, DIM>::NodeIterator node_iter = rCellPopulation.rGetMesh().GetNodeIteratorBegin();
         node_iter != rCellPopulation.rGetMesh().GetNodeIteratorEnd();
         ++node_iter)
    {
        CellPtr p_cell = rCellPopulation.GetCellUsingLocationIndex(node_iter->GetIndex());

        c_vector<double, DIM> applied_velocity;
        
        switch (DIM)
        {
            case 1:
                applied_velocity[0] = p_cell->GetCellData()->GetItem("applied_velocity_x");
                break;
            case 2:
                applied_velocity[0] = p_cell->GetCellData()->GetItem("applied_velocity_x");
                applied_velocity[1] = p_cell->GetCellData()->GetItem("applied_velocity_y");
                break;
            case 3:
                applied_velocity[0] = p_cell->GetCellData()->GetItem("applied_velocity_x");
                applied_velocity[1] = p_cell->GetCellData()->GetItem("applied_velocity_y");
                applied_velocity[2] = p_cell->GetCellData()->GetItem("applied_velocity_z");
                break;
            default:
                NEVER_REACHED;
        }
        
        c_vector<double, DIM> force_contribution;
        for (unsigned i=0; i<DIM; i++)
        {
            force_contribution[i] = applied_velocity[i]/dt;
        }
        node_iter->AddAppliedForceContribution(force_contribution);
    }
}

template<unsigned DIM>
void AppliedVelocityForce<DIM>::OutputForceParameters(out_stream& rParamsFile)
{
    *rParamsFile << "\t\t\t<MovementParameter>" << mMovementParameter << "</MovementParameter> \n";

    // Call direct parent class
    AbstractForce<DIM>::OutputForceParameters(rParamsFile);
}

/////////////////////////////////////////////////////////////////////////////
// Explicit instantiation
/////////////////////////////////////////////////////////////////////////////

template class AppliedVelocityForce<1>;
template class AppliedVelocityForce<2>;
template class AppliedVelocityForce<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(AppliedVelocityForce)
