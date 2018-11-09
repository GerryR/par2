/**
* @file PParticles.cuh
* @brief Header file for PParticles class.
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

#ifndef PAR2_PPARTICLES_CUH
#define PAR2_PPARTICLES_CUH

#include <thrust/device_vector.h>
#include <thrust/host_vector.h>
#include "../Geometry/CartesianGrid.cuh"
#include "../Geometry/FaceField.cuh"
#include "MoveParticle.cuh"
#include <curand_kernel.h>

namespace par2
{
    /**
    * @class PParticles
    * @brief Class managing a cloud of particles. It provides functions
    *        to initialize and move the particles on device.
    * @tparam T Float number precision
    */
    template<typename T>
    class PParticles
    {
    public:

        /**
        * @brief Constructor.
        * @param _grid Grid where the velocity is defined
        * @param _datax Velocity vector (x-direction)
        * @param _datay Velocity vector (y-direction)
        * @param _dataz Velocity vector (z-direction)
        * @param _molecularDiffusion Effective molecular diffusion
        * @param _alphaL Longitudinal dispersivity
        * @param _alphaT Transverse dispersivity
        * @param _nParticles Total number of particles
        * @param _seed Seed for pseudo-random number generator
        * @param _useTrilinearCorrection True if trilinear correction is used
        */
        PParticles(const grid::Grid<T>& _grid,
                   const thrust::host_vector<T> &_datax,
                   const thrust::host_vector<T> &_datay,
                   const thrust::host_vector<T> &_dataz,
                   T _molecularDiffusion,
                   T _alphaL,
                   T _alphaT,
                   unsigned int _nParticles,
                   long int _seed = time(NULL),
                   bool _useTrilinearCorrection = true);

        /**
        * @brief Number of particles in the system.
        * @return Number of particles
        */
        unsigned int size() const;

        /**
        * @brief Initialize particles inside a box defined by
        *        two points p1 and p2. The particles are uniformly
        *        distributed inside the box.
        * @param p1x x-component of p1
        * @param p1y y-component of p1
        * @param p1z z-component of p1
        * @param p2x x-component of p2
        * @param p2y y-component of p2
        * @param p2z z-component of p2
        */
        void initializeBox(T p1x, T p1y, T p1z,
                              T p2x, T p2y, T p2z);

        /**
        * @brief Execute one step of the particle tracking method using
        *        a time step dt.
        * @param dt Time step
        */
        void move(T dt);

        /**
        * @brief Export all the particles in a csv file.
        *        WARNING: this function is computationally expensive
        *        and should be avoided if computation time is critical.
        * @param fileName Path to the csv file
        */
        void exportCSV(const std::string &fileName) const;

        /**
        * @brief Compute the percentage of particles inside a box
        *        defined by two points p1 and p2.
        * @param p1x x-component of p1
        * @param p1y y-component of p1
        * @param p1z z-component of p1
        * @param p2x x-component of p2
        * @param p2y y-component of p2
        * @param p2z z-component of p2
        * @return Percentage of particles inside the box.
        */
        T concentrationBox(T p1x, T p1y, T p1z,
                           T p2x, T p2y, T p2z) const;

        /**
        * @brief Compute the percentage of particles that crossed
        *        a plane orthogonal to the x-axis (i.e., px > xplane).
        * @param xplane Location of the plane
        * @return Percentage of particles after the plane.
        */
        T concentrationAfterX(T xplane) const;

    private:

        // Thrust functor for one step of particle tracking
        MoveParticle<T> moveParticle;
        // Grid and physical variables
        grid::Grid<T> grid;
        T molecularDiffusion;
        T alphaL;
        T alphaT;
        unsigned int nParticles;
        // Velocity field stored on device (facefield)
        thrust::device_vector<T> datax;
        thrust::device_vector<T> datay;
        thrust::device_vector<T> dataz;
        // Velocity field stored on device (cornerfield or cellfield)
        thrust::device_vector<T> cdatax;
        thrust::device_vector<T> cdatay;
        thrust::device_vector<T> cdataz;
        // Particle positions stored on device
        thrust::device_vector<T> cx;
        thrust::device_vector<T> cy;
        thrust::device_vector<T> cz;
        // Vector of curand states stored on device
        thrust::device_vector<curandState_t> states;
        // If true, use trilinear interpolation
        bool useTrilinearCorrection;
        // Max number of particles simulated in each kernel
        const int maxParticles = 65536;
    };

}

#include "PParticles.cu"

#endif //PAR2_PPARTICLES_CUH
