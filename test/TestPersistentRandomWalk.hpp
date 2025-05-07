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
#include "AbstractCellBasedWithTimingsTestSuite.hpp"
#include "CellVolumesWriter.hpp"
#include "CellsGenerator.hpp"
#include "FixedG1GenerationalCellCycleModel.hpp"
#include "ForwardEulerNumericalMethod.hpp"
#include "NodeBasedCellPopulation.hpp"
#include "OffLatticeSimulation.hpp"
#include "SmartPointers.hpp"
#include "DifferentiatedCellProliferativeType.hpp"

#include "NodeVelocityWriter.hpp"


#include "PetscSetupAndFinalize.hpp"
#include "Debug.hpp"



class Test01PersistentRandomWalk: public AbstractCellBasedWithTimingsTestSuite
{
private:

public:

    /**
     * Simple node based random walk of N cells 
     */
    void TestNodeBasedRandomWalk()
    {
        EXIT_IF_PARALLEL; // Cant access cells with index loop on multiple processors.

        std::string base_name = "OpenVT/Test01PersistentRandomWalk/";

        unsigned num_cells = 10;
        double end_time = 1;
        double dt = 1.0/100.0;

        std::string model_types[2] = {"Model001","Model003"};

        for (unsigned model_type_index = 0; model_type_index != 1; model_type_index++)
        {
            std::string model_type = model_types[model_type_index];

            std::string output_dir = base_name + model_type;
        
            PRINT_VARIABLE(output_dir);

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
            simulator.SetEndTime(end_time);
            simulator.SetUpdateCellPopulationRule(false);

            // Create a Random force law (i.e diffusion) and pass it to the simulation
            MAKE_PTR(RandomMotionForce<2>, p_diffusion_force);
            simulator.AddForce(p_diffusion_force);

            simulator.Solve();

            // Reset for next simulation
            SimulationTime::Instance()->Destroy();
            SimulationTime::Instance()->SetStartTime(0.0);
        }
    }
        
        
};

#endif /* TEST01PERSISTENTRANDOMWALK_HPP_ */
