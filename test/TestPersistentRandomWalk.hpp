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

#ifndef TEST01PERSISTENTRANDOMWALK_HPP_
#define TEST01PERSISTENTRANDOMWALK_HPP_

#include <cxxtest/TestSuite.h>

// Must be included before other cell_based headers
#include "CellBasedSimulationArchiver.hpp"

#include "AbstractCellBasedWithTimingsTestSuite.hpp"

#include "RandomMotionForce.hpp"
#include "AppliedVelocityForce.hpp"
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellVolumesWriter.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"
#include "PersistentRandomWalkModifier.hpp"

#include "NodeVelocityWriter.hpp"


#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"



class Test01PersistentRandomWalk: public AbstractCellBasedWithTimingsTestSuite
{
private:

public:
 
    // Edit this to take in XML file to specify parameters or use string generated for JSON.


    /**
     * Simple node based random walk of N cells 
     */
    void TestNodeBasedRandomWalk()
    {
        EXIT_IF_PARALLEL; // Cant access cells with index loop on multiple processors.

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-num_cells"));
        unsigned num_cells = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-num_cells");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-sampling_rate"));
        unsigned sampling_rate = CommandLineArguments::Instance()->GetUnsignedCorrespondingToOption("-sampling_rate");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-end_time"));
        double end_time = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-end_time");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-dt"));
        double dt = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-dt");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-model_type"));
        std::string model_type = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-model_type");

        TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-output_name"));
        std::string output_name = CommandLineArguments::Instance()->GetStringCorrespondingToOption("-output_name");

        std::string base_name = "OpenVT/Test01PersistentRandomWalk/";
        std::string output_dir = base_name + output_name;
    
        PRINT_2_VARIABLES("Number of cells", num_cells);
        PRINT_2_VARIABLES("End time", end_time);
        PRINT_2_VARIABLES("Time step", dt);
        PRINT_2_VARIABLES("Model type", model_type);
        PRINT_2_VARIABLES("Output name", output_name);

        bool is_model_001 = (model_type.compare("MODEL001")==0);
        bool is_model_004 = (model_type.compare("MODEL004")==0);
        PRINT_2_VARIABLES("Is this model 001?", is_model_001);
        PRINT_2_VARIABLES("Is this model 004?", is_model_004);

        // Create a simple mesh
        std::vector<Node<2>*> nodes (num_cells);

        for ( unsigned j = 0; j < num_cells; j++ )
        {
                nodes[j] = new Node<2>(j, false, 0.0,0.0);
        }

        // Convert this to a NodesOnlyMesh
        NodesOnlyMesh<2> mesh;
        mesh.ConstructNodesWithoutMesh(nodes, 1.5);

        // Create cells
        std::vector<CellPtr> cells;
        MAKE_PTR(DifferentiatedCellProliferativeType, p_differentiated_type);
        CellsGenerator<FixedG1GenerationalCellCycleModel, 2> cells_generator;
        cells_generator.GenerateBasicRandom(cells, mesh.GetNumNodes(), p_differentiated_type);

        // Create a node-based cell population
        NodeBasedCellPopulation<2> node_based_cell_population(mesh, cells);

        // Set up cell-based simulation
        OffLatticeSimulation<2> simulator(node_based_cell_population);
        simulator.SetOutputDirectory(output_dir);
        simulator.SetDt(dt);
        simulator.SetSamplingTimestepMultiple(sampling_rate);
        simulator.SetEndTime(end_time);
        simulator.SetUpdateCellPopulationRule(false);

        if(model_type.compare("MODEL001")==0) 
        {
            TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-movement_parameter"));
            double movement_parameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-movement_parameter");

            // Create a Random force law (i.e diffusion) and pass it to the simulation
            MAKE_PTR(RandomMotionForce<2>, p_diffusion_force);
            p_diffusion_force->SetMovementParameter(movement_parameter);
            simulator.AddForce(p_diffusion_force);

            simulator.Solve();

            // Reset for next simulation
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
        else if(model_type.compare("MODEL004")==0) 
        {
            TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-velocity_change_rate"));
            double velocity_change_rate = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-velocity_change_rate");

            TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-motility_parameter"));
            double motility_parameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-motility_parameter");

            TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-bias_parameter"));
            double bias_parameter = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-bias_parameter");

            TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-bias_direction_x"));
            double bias_direction_x = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-bias_direction_x");

            TS_ASSERT(CommandLineArguments::Instance()->OptionExists("-bias_direction_y"));
            double bias_direction_y = CommandLineArguments::Instance()->GetDoubleCorrespondingToOption("-bias_direction_y");

            c_vector<double, 2> bias_direction;
            bias_direction[0] = bias_direction_x;
            bias_direction[1] = bias_direction_y;
            MAKE_PTR(AppliedVelocityForce<2>, p_applied_velocity_force);
            simulator.AddForce(p_applied_velocity_force);
            
            MAKE_PTR(PersistentRandomWalkModifier<2>, p_modifier);
            // TODO set parameters on modifier
            p_modifier->SetVelocityChangeRate(velocity_change_rate);
            p_modifier->SetMotilityParameter(motility_parameter);
            p_modifier->SetBiasParameter(bias_parameter);
            p_modifier->SetBiasDirection(bias_direction);
            simulator.AddSimulationModifier(p_modifier);

            simulator.Solve();

            // Reset for next simulation
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
        else
        {
            NEVER_REACHED;
        }

       
    }
        
        
};

#endif /* TEST01PERSISTENTRANDOMWALK_HPP_ */
