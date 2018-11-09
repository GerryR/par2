/**
* @file maingpu.cu
* @brief Entry point for PAR2
*
* @author Calogero B. Rizzo
*
* @copyright This file is part of the PAR2 software.
*            Copyright (C) 2018 Calogero B. Rizzo
*
* @license This program is free software: you can redistribute it and/or modify
*          it under the terms of the GNU General Public License as published by
*          the Free Software Foundation, either version 3 of the License, or
*          (at your option) any later version.
*
*          This program is distributed in the hope that it will be useful,
*          but WITHOUT ANY WARRANTY; without even the implied warranty of
*          MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*          GNU General Public License for more details.
*
*          You should have received a copy of the GNU General Public License
*          along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/
#include <iostream>
#include <iomanip>
#include <sstream>
#include <chrono>
#include <exception>
#include "Geometry/Point.cuh"
#include "Geometry/CartesianGrid.cuh"
#include "Geometry/FaceField.cuh"
#include "Geometry/CellField.cuh"
#include "Particles/PParticles.cuh"
#include "Utilities/Parameters.h"
#include "Utilities/spdlog/spdlog.h"

#include <thrust/host_vector.h>

namespace spd = spdlog;

void run(int argc, char** argv)
{
    auto console = spd::stdout_color_mt("console");
    console->set_pattern("%v");

    console->info("*********************************************************");
    console->info("*-------------------------------------------------------*");
    console->info("*------------------------- PAR2 ------------------------*");
    console->info("*-------------------------------------------------------*");
    console->info("*********************************************************");
    console->info("");

    std::string configurationFile;

    if (argc == 1)
    {
        console->info("Insert path to YAML configuration file:");
        std::cin >> configurationFile;
        console->info("");
    }
    else if (argc == 2)
    {
        configurationFile = std::string(argv[1]);
    }
    else
    {
        std::stringstream errorMessage;
        errorMessage << "Specify path to the YAML configuration file: "
                     << "(> par2 /path/to/configuration.yaml)";
        throw std::runtime_error(errorMessage.str());
    }

    std::string configurationPath =
        configurationFile.substr(0, configurationFile.find_last_of("/\\") + 1);

    auto logfile = spd::basic_logger_mt("logfile",
                configurationPath + std::string("par2.log"), true);

    logfile->info("Logfile created {}", configurationFile);

    // Start time
    auto t_start = std::chrono::system_clock::now();

    // Load configuration file
    par2::Parameters<PAR2_FLOAT> par(configurationFile);

    // Create grid
    auto grid = par2::grid::build<PAR2_FLOAT>(par.nx(), par.ny(), par.nz(),
                                         par.dx(), par.dy(), par.dz());
    logfile->info("Grid size: {} x {} x {}", par.nx(), par.ny(), par.nz());
    logfile->info("Cell size: {} x {} x {}", par.dx(), par.dy(), par.dz());

    // Create and load velocity field
    console->info("Import from MODFLOW...");
    thrust::host_vector<PAR2_FLOAT> datax, datay, dataz;
    par2::facefield::build<PAR2_FLOAT>(grid, datax, datay, dataz);

    if (par.velType() == "modflow")
    {
        par2::facefield::importFromModflow(grid, datax, datay, dataz,
                                          configurationPath + par.velPath(),
                                          par.rho());
        logfile->info("FTL file: {}", configurationPath + par.velPath());
        logfile->info("Porosity: {}", par.rho());
    }
    else
    {
        throw par2::yaml_invalid_argument("Velocity type must be 'modflow'.");
    }

    // Choose interpolation method
    bool cTrilinear = false;
    if (par.interp() == "trilinear")
    {
        cTrilinear = true;
    }
    else if (par.interp() == "finite difference")
    {
        cTrilinear = false;
    }
    else
    {
        throw par2::yaml_invalid_argument("Interpolation must be either 'trilinear' or 'finite difference'.");
    }
    logfile->info("Trilinear: {}", cTrilinear);

    // Create and initialize particles
    console->info("Device initialization...");
    long int seed = par.seed();
    logfile->info("Seed: {}", seed);
    logfile->info("Molecular diffusion: {}", par.Dm());
    logfile->info("Longitudinal dispersivity: {}", par.alphaL());
    logfile->info("Transverse dispersivity: {}", par.alphaT());
    logfile->info("Particles: {}", par.nParticles());
    par2::PParticles<PAR2_FLOAT> particles(grid, datax, datay, dataz,
                                     par.Dm(), par.alphaL(), par.alphaT(),
                                     par.nParticles(), seed, cTrilinear);

    logfile->info("Box P1: ({}, {}, {})", par.p1x(), par.p1y(), par.p1z());
    logfile->info("Box P2: ({}, {}, {})", par.p2x(), par.p2y(), par.p2z());
    particles.initializeBox(par.p1x(), par.p1y(), par.p1z(),
                            par.p2x(), par.p2y(), par.p2z());

    // Initialize CSV output files
    std::ofstream csvOutStream;
    if (par.csvOutput())
    {
        csvOutStream.open(configurationPath + par.csvPath());

        // Precision of the output
        csvOutStream << std::setprecision(15) << std::fixed;

        if (csvOutStream.is_open())
        {
            console->info("Prepare output ({})...", configurationPath + par.csvPath());

            // Write labels in the first row
            // First two columns contain the current step and time
            csvOutStream << "step, time";
            for (auto i = 0; i < par.csvNumberOfItems(); i++)
            {
                csvOutStream << ", " << par.csvItemLabel(i);
            }
            csvOutStream << std::endl;

            logfile->info("{} ready", par.csvPath());
        }
        else
        {
            throw std::runtime_error(std::string("Could not open file ") +
                        configurationPath + par.csvPath());
        }
    }

    // Keep track of the next step where we want to output a full particles
    // snapshot. If stepSnapshot=-1, no snapshots will be created.
    int stepSnapshot = -1;
    int stepId       = -1;

    if (par.snapshotOutput())
    {
        if (par.snapshotUseSkip())
        {
            stepSnapshot = 0;
        }
        else if (par.snapshotSize() > 0)
        {
            stepId = 0;
            stepSnapshot = par.snapshotStep(stepId);
        }
    }

    // Ready time
    auto t_ready = std::chrono::system_clock::now();

    // Simulation variables
    int steps = par.steps();
    PAR2_FLOAT dt = par.dt();
    PAR2_FLOAT completed = 0.0;

    console->info("");
    console->info("Start simulation...");
    // Start time loop
    // Using steps+1 to take into account the last step.
    for (auto step = 0; step < steps+1; step++)
    {
        // CSV Output
        if (par.csvOutput())
        {
            if (step%par.csvSkip() == 0)
            {
                logfile->info("Writing CSV file (STEP {})", step);

                // Write step and time
                csvOutStream << step << ", " << step*dt;
                for (auto i = 0; i < par.csvNumberOfItems(); i++)
                {
                    csvOutStream << ", ";
                    if (!par.csvItemType(i).compare("after-x"))
                    {
                        csvOutStream << particles.concentrationAfterX(par.csvItemX(i));
                    }
                    else if (!par.csvItemType(i).compare("box"))
                    {
                        csvOutStream << particles.concentrationBox(
                                            par.csvItemP1X(i),
                                            par.csvItemP1Y(i),
                                            par.csvItemP1Z(i),
                                            par.csvItemP2X(i),
                                            par.csvItemP2Y(i),
                                            par.csvItemP2Z(i));
                    }
                    else
                    {
                        throw par2::yaml_invalid_argument("CSV Item type must be either 'after-x' or 'box'.");
                    }
                }
                csvOutStream << std::endl;
            }
        }

        // Snapshot Output
        if (step == stepSnapshot)
        {
            // Export particle positions
            auto snapshotOutputPath = configurationPath + par.snapshotPath(step);

            logfile->info("Writing Snapshot file {} (STEP {})", snapshotOutputPath, step);
            console->info("Writing Snapshot file {}", snapshotOutputPath);

            particles.exportCSV(snapshotOutputPath);

            // Update next snapshot
            if (par.snapshotUseSkip())
            {
                stepSnapshot += par.snapshotSkip();
            }
            else
            {
                stepId++;
                if (stepId == par.snapshotSize())
                {
                    stepSnapshot = -1;
                }
                else
                {
                    stepSnapshot = par.snapshotStep(stepId);
                }
            }
        }

        // Move particles
        particles.move(dt);

        if (step >= completed)
        {
            console->info("{:3.0f}% completed (STEP {})", completed/steps*100.0, step);
            completed += steps/50.0;
        }
    }

    console->info("Simulation completed");
    console->info("");

    // Close CSV file
    if (par.csvOutput())
    {
        csvOutStream.close();
    }

    // End time
    auto t_end = std::chrono::system_clock::now();

    // Print information about computation time
    std::chrono::duration<PAR2_FLOAT> preproc_seconds    = t_ready - t_start;
    std::chrono::duration<PAR2_FLOAT> simulation_seconds = t_end - t_ready;
    std::chrono::duration<PAR2_FLOAT> elapsed_seconds    = t_end - t_start;

    logfile->info("TIME PREPROCESSING: {}s", preproc_seconds.count());
    logfile->info("TIME SIMULATION: {}s", simulation_seconds.count());
    logfile->info("TIME ELAPSED: {}s", elapsed_seconds.count());

    console->info("TIME PREPROCESSING: {}s", preproc_seconds.count());
    console->info("TIME SIMULATION: {}s", simulation_seconds.count());
    console->info("TIME ELAPSED: {}s", elapsed_seconds.count());
}

int main(int argc, char** argv)
{
    int EXIT = EXIT_SUCCESS;

    try
    {
        run(argc, argv);
    }
    catch (const std::exception& e)
    {
        auto console = spdlog::get("console");
        if (console)
        {
            console->critical("Execution terminated with an error:");
            console->critical(e.what());
        }
        auto logfile = spdlog::get("logfile");
        if (logfile)
        {
            logfile->critical("Execution terminated with an error:");
            logfile->critical(e.what());
        }
        EXIT = EXIT_FAILURE;
    }
    // Release and close all loggers
    spd::drop_all();

    return EXIT;
}
