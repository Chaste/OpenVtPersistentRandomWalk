/*

Copyright (c) 2005-2025, University of Oxford.
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

#include "PersistentRandomWalkModifier.hpp"
#include "RandomNumberGenerator.hpp"

template<unsigned DIM>
PersistentRandomWalkModifier<DIM>::PersistentRandomWalkModifier()
    : AbstractCellBasedSimulationModifier<DIM>(),
      mVelocityChangeRate(10.0),
      mMotilityParmeter(0.1),
      mBiasParmeter(0.5),
      mBiasDirection(unit_vector<double>(DIM, 0))
{
}

template<unsigned DIM>
PersistentRandomWalkModifier<DIM>::~PersistentRandomWalkModifier()
{
}

template<unsigned DIM>
void PersistentRandomWalkModifier<DIM>::UpdateAtEndOfTimeStep(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    UpdateCellData(rCellPopulation);
}

template<unsigned DIM>
void PersistentRandomWalkModifier<DIM>::SetupSolve(AbstractCellPopulation<DIM,DIM>& rCellPopulation, std::string outputDirectory)
{
    c_vector<double, DIM> applied_velocity = mMotilityParmeter*mBiasDirection;

    // Next iterate over the population and update the velocity if needed
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {
        switch (DIM)
        {
            case 1:
                cell_iter->GetCellData()->SetItem("applied_velocity_x",applied_velocity[0]);
                break;
            case 2:
                cell_iter->GetCellData()->SetItem("applied_velocity_x",applied_velocity[0]);
                cell_iter->GetCellData()->SetItem("applied_velocity_y",applied_velocity[1]);
                break;
            case 3:
                cell_iter->GetCellData()->SetItem("applied_velocity_x",applied_velocity[0]);
                cell_iter->GetCellData()->SetItem("applied_velocity_y",applied_velocity[1]);
                cell_iter->GetCellData()->SetItem("applied_velocity_z",applied_velocity[2]);
                break;
            default:
                NEVER_REACHED;
        }
    }
}

template<unsigned DIM>
void PersistentRandomWalkModifier<DIM>::UpdateCellData(AbstractCellPopulation<DIM,DIM>& rCellPopulation)
{
    // Make sure the cell population is updated
    rCellPopulation.Update();

    double dt = SimulationTime::Instance()->GetTimeStep();

    // Next iterate over the population and update the velocity if needed
    for (typename AbstractCellPopulation<DIM>::Iterator cell_iter = rCellPopulation.Begin();
         cell_iter != rCellPopulation.End();
         ++cell_iter)
    {

        double probability_update = dt*mVelocityChangeRate;

        // Sample ot see if updating velocity 
        double u = RandomNumberGenerator::Instance()->ranf();
        
        if (u < probability_update) // Update the velocity for this cell
        {
            // Pick a random direction 
            
            c_vector<double, DIM> random_vector;

            switch (DIM)
            {
                case 1:
                {
                    double random_direction = -1.0 + 2.0*(RandomNumberGenerator::Instance()->ranf() < 0.5);

                    random_vector(0) = random_direction;
                    break;
                }
                case 2:
                {
                    double random_angle = 2.0*M_PI*RandomNumberGenerator::Instance()->ranf();

                    random_vector(0) = cos(random_angle);
                    random_vector(1) = sin(random_angle);
                    break;
                }
                case 3:
                {
                    /*
                    * Note that to pick a random point on the surface of a sphere, it is incorrect
                    * to select spherical coordinates from uniform distributions on [0, 2*pi) and
                    * [0, pi) respectively, since points picked in this way will be 'bunched' near
                    * the poles. See #2230.
                    */
                    double u = RandomNumberGenerator::Instance()->ranf();
                    double v = RandomNumberGenerator::Instance()->ranf();

                    double random_azimuth_angle = 2*M_PI*u;
                    double random_zenith_angle = std::acos(2*v - 1);

                    random_vector(0) = cos(random_azimuth_angle)*sin(random_zenith_angle);
                    random_vector(1) = sin(random_azimuth_angle)*sin(random_zenith_angle);
                    random_vector(2) = cos(random_zenith_angle);
                    break;
                }
                default:
                    // This can't happen
                    NEVER_REACHED;
            }


            c_vector<double, DIM> applied_velocity = (1-mBiasParmeter)*random_vector + mBiasParmeter*mBiasDirection;

            applied_velocity = applied_velocity / norm_2(applied_velocity);
            applied_velocity *= mMotilityParmeter;
            // Now we have the new velocity, we need to store it in CellData.
        
            switch (DIM)
            {
                case 1:
                    cell_iter->GetCellData()->SetItem("applied_velocity_x",applied_velocity[0]);
                    break;
                case 2:
                    cell_iter->GetCellData()->SetItem("applied_velocity_x",applied_velocity[0]);
                    cell_iter->GetCellData()->SetItem("applied_velocity_y",applied_velocity[1]);
                    break;
                case 3:
                    cell_iter->GetCellData()->SetItem("applied_velocity_x",applied_velocity[0]);
                    cell_iter->GetCellData()->SetItem("applied_velocity_y",applied_velocity[1]);
                    cell_iter->GetCellData()->SetItem("applied_velocity_z",applied_velocity[2]);
                    break;
                default:
                    NEVER_REACHED;
            }
        }
    }
}

template<unsigned DIM>
void PersistentRandomWalkModifier<DIM>::OutputSimulationModifierParameters(out_stream& rParamsFile)
{
    // No parameters to output, so just call method on direct parent class
    AbstractCellBasedSimulationModifier<DIM>::OutputSimulationModifierParameters(rParamsFile);
}

// Explicit instantiation
template class PersistentRandomWalkModifier<1>;
template class PersistentRandomWalkModifier<2>;
template class PersistentRandomWalkModifier<3>;

// Serialization for Boost >= 1.36
#include "SerializationExportWrapperForCpp.hpp"
EXPORT_TEMPLATE_CLASS_SAME_DIMS(PersistentRandomWalkModifier)
