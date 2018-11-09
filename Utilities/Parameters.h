/**
* @file Parameters.h
* @brief Parse the YAML configuration file using the yaml-cpp library
*        and create the simulation parameters.
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

#ifndef PAR2_PARAMETERS_H
#define PAR2_PARAMETERS_H

#include <yaml-cpp/yaml.h>
#include <sstream>
#include <string>
#include <ctime>
#include <iostream>
#include <exception>

namespace par2
{
    /**
    * @struct yaml_invalid_argument
    * @brief Custom exception for YAML invalid argument.
    */
    struct yaml_invalid_argument : public std::runtime_error
    {
        using std::runtime_error::runtime_error;
    };

    /**
    * @class Parameters
    * @brief Interface between the YAML input file and PAR2.
    * @tparam T Float number precision
    */
    template<typename T>
    class Parameters
    {
    public:

        /**
        * @brief Constructor.
        * @param fileName Path to the YAML configuration file.
        */
        Parameters(const std::string& fileName)
        {
            try
            {
                config = YAML::LoadFile(fileName);
            }
            catch(const std::exception& e)
            {
                std::stringstream errorMessage;
                errorMessage << "Error in the YAML configuration file "
                             << fileName << " (" << e.what() << ")";
                throw std::runtime_error(errorMessage.str());
            }
        }

        // GRID PARAMETERS

        /**
        * @brief Number of cells in x direction.
        */
        int nx() const
        {
            return config["grid"]["dimension"][0].as<int>();
        }

        /**
        * @brief Number of cells in y direction.
        */
        int ny() const
        {
            return config["grid"]["dimension"][1].as<int>();
        }

        /**
        * @brief Number of cells in z direction.
        */
        int nz() const
        {
            return config["grid"]["dimension"][2].as<int>();
        }

        /**
        * @brief Size of cell in x direction.
        */
        T dx() const
        {
            return config["grid"]["cell size"][0].as<T>();
        }

        /**
        * @brief Size of cell in y direction.
        */
        T dy() const
        {
            return config["grid"]["cell size"][1].as<T>();
        }

        /**
        * @brief Size of cell in z direction.
        */
        T dz() const
        {
            return config["grid"]["cell size"][2].as<T>();
        }

        // PHYSICS PARAMETERS

        /**
        * @brief Effective porosity.
        */
        T rho() const
        {
            return config["physics"]["porosity"].as<T>();
        }

        /**
        * @brief Effective molecular diffusion.
        */
        T Dm() const
        {
            return config["physics"]["molecular diffusion"].as<T>();
        }

        /**
        * @brief Longitudinal dispersivity.
        */
        T alphaL() const
        {
            return config["physics"]["longitudinal dispersivity"].as<T>();
        }

        /**
        * @brief Transverse dispersivity.
        */
        T alphaT() const
        {
            return config["physics"]["transverse dispersivity"].as<T>();
        }

        /**
        * @brief Type of velocity to import.
        */
        std::string velType() const
        {
            return config["physics"]["velocity"]["type"].as<std::string>();
        }

        /**
        * @brief Path to the file containing the velocity.
        */
        std::string velPath() const
        {
            return config["physics"]["velocity"]["file"].as<std::string>();
        }

        // SIMULATION PARAMETERS

        /**
        * @brief Number of particles.
        */
        int nParticles() const
        {
            return config["simulation"]["particles"]["N"].as<int>();
        }

        /**
        * @brief First point of the initial volume containing the particles (x).
        */
        T p1x() const
        {
            return config["simulation"]["particles"]["start"]["p1"][0].as<T>();
        }

        /**
        * @brief First point of the initial volume containing the particles (y).
        */
        T p1y() const
        {
            return config["simulation"]["particles"]["start"]["p1"][1].as<T>();
        }

        /**
        * @brief First point of the initial volume containing the particles (z).
        */
        T p1z() const
        {
            return config["simulation"]["particles"]["start"]["p1"][2].as<T>();
        }

        /**
        * @brief Second point of the initial volume containing the particles (x).
        */
        T p2x() const
        {
            return config["simulation"]["particles"]["start"]["p2"][0].as<T>();
        }

        /**
        * @brief Second point of the initial volume containing the particles (y).
        */
        T p2y() const
        {
            return config["simulation"]["particles"]["start"]["p2"][1].as<T>();
        }

        /**
        * @brief Second point of the initial volume containing the particles (z).
        */
        T p2z() const
        {
            return config["simulation"]["particles"]["start"]["p2"][2].as<T>();
        }

        /**
        * @brief Interpolation method.
        */
        std::string interp() const
        {
            if (config["simulation"]["interpolation"])
            {
                return config["simulation"]["interpolation"].as<std::string>();
            }
            return std::string("trilinear");
        }

        /**
        * @brief Time step.
        */
        T dt() const
        {
            return config["simulation"]["dt"].as<T>();
        }

        /**
        * @brief Number of steps.
        */
        int steps() const
        {
            return config["simulation"]["steps"].as<int>();
        }

        /**
        * @brief Seed for pseudo-random number generator.
        */
        long int seed() const
        {
            if (config["simulation"]["seed"])
            {
                return config["simulation"]["seed"].as<long int>();
            }
            return std::time(NULL);
        }

        // OUTPUT PARAMETERS

        /**
        * @brief Check if output to CSV is required.
        */
        bool csvOutput() const
        {
            if (config["output"]["csv"])
                return true;
            return false;
        }

        /**
        * @brief CSV file path.
        */
        std::string csvPath() const
        {
            return config["output"]["csv"]["file"].as<std::string>();
        }

        /**
        * @brief Output every 'skip' steps.
        */
        int csvSkip() const
        {
            return config["output"]["csv"]["skip"].as<int>();
        }

        /**
        * @brief Number of items to export.
        */
        int csvNumberOfItems() const
        {
            return config["output"]["csv"]["items"].size();
        }

        /**
        * @brief Label of the output item.
        * @param i ID of the item
        */
        std::string csvItemLabel(int i) const
        {
            assert(i < csvNumberOfItems());
            return config["output"]["csv"]["items"][i]["label"].as<std::string>();
        }

        /**
        * @brief Label of the output item.
        * @param i ID of the item
        */
        std::string csvItemType(int i) const
        {
            assert(i < csvNumberOfItems());
            return config["output"]["csv"]["items"][i]["type"].as<std::string>();
        }

        /**
        * @brief Value x of the output item. Valid only if type=='after-x'.
        * @param i ID of the item
        */
        T csvItemX(int i) const
        {
            assert(!csvItemType(i).compare("after-x"));
            return config["output"]["csv"]["items"][i]["x"].as<T>();
        }

        /**
        * @brief Value p1x of the output item. Valid only if type=='box'.
        * @param i ID of the item
        */
        T csvItemP1X(int i) const
        {
            assert(!csvItemType(i).compare("box"));
            return config["output"]["csv"]["items"][i]["p1"][0].as<T>();
        }

        /**
        * @brief Value p1y of the output item. Valid only if type=='box'.
        * @param i ID of the item
        */
        T csvItemP1Y(int i) const
        {
            assert(!csvItemType(i).compare("box"));
            return config["output"]["csv"]["items"][i]["p1"][1].as<T>();
        }

        /**
        * @brief Value p1z of the output item. Valid only if type=='box'.
        * @param i ID of the item
        */
        T csvItemP1Z(int i) const
        {
            assert(!csvItemType(i).compare("box"));
            return config["output"]["csv"]["items"][i]["p1"][2].as<T>();
        }

        /**
        * @brief Value p2x of the output item. Valid only if type=='box'.
        * @param i ID of the item
        */
        T csvItemP2X(int i) const
        {
            assert(!csvItemType(i).compare("box"));
            return config["output"]["csv"]["items"][i]["p2"][0].as<T>();
        }

        /**
        * @brief Value p2y of the output item. Valid only if type=='box'.
        * @param i ID of the item
        */
        T csvItemP2Y(int i) const
        {
            assert(!csvItemType(i).compare("box"));
            return config["output"]["csv"]["items"][i]["p2"][1].as<T>();
        }

        /**
        * @brief Value p2z of the output item. Valid only if type=='box'.
        * @param i ID of the item
        */
        T csvItemP2Z(int i) const
        {
            assert(!csvItemType(i).compare("box"));
            return config["output"]["csv"]["items"][i]["p2"][2].as<T>();
        }

        /**
        * @brief Check if snapshot output is required.
        */
        bool snapshotOutput() const
        {
            if (config["output"]["snapshot"])
                return true;
            return false;
        }

        /**
        * @brief Check if we are using skip or the explicit vector for snapshots.
        */
        bool snapshotUseSkip() const
        {
            if (config["output"]["snapshot"]["skip"])
                return true;
            return false;
        }

        /**
        * @brief Output snapshot every 'skip' steps.
        */
        int snapshotSkip() const
        {
            return config["output"]["snapshot"]["skip"].as<int>();
        }

        /**
        * @brief Number of snapshots to output.
        */
        int snapshotSize() const
        {
            return config["output"]["snapshot"]["steps"].size();
        }

        /**
        * @brief Step where a snapshot is taken.
        * @param i Position in the steps list
        */
        int snapshotStep(int i) const
        {
            assert(i < snapshotSize());
            return config["output"]["snapshot"]["steps"][i].as<int>();
        }

        /**
        * @brief Snapshot file path for the given step.
        * @param step Current step
        */
        std::string snapshotPath(int step) const
        {
            auto path = config["output"]["snapshot"]["file"].as<std::string>();

            // replace the * with the current step
            path.replace(path.find("*"), 1, std::to_string(step));

            return path;
        }


    private:

        YAML::Node config;

    };
}

#endif //PAR2_PARAMETERS_H
