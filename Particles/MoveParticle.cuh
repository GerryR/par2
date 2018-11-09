/**
* @file PParticles.cu
* @brief MoveParticle functor.
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

#ifndef PAR2_MOVEPARTICLE_CUH
#define PAR2_MOVEPARTICLE_CUH

#include <thrust/tuple.h>
#include <thrust/random/linear_congruential_engine.h>
#include <thrust/random/normal_distribution.h>
#include <curand_kernel.h>
#include "../Geometry/CartesianGrid.cuh"
#include "../Geometry/FaceField.cuh"
#include "../Geometry/CornerField.cuh"
#include "../Geometry/Vector.cuh"

namespace par2
{
    /**
    * @struct MoveParticle
    * @brief Thrust functor for one step of the particle tracking method.
    * @tparam T Float number precision
    */
    template<typename T>
    struct MoveParticle
    {
        // Raw pointer to velocity vectors (facefield)
        T* datax;
        T* datay;
        T* dataz;
        // Raw pointer to velocity vectors (cornerfield or cellfield)
        T* cdatax;
        T* cdatay;
        T* cdataz;
        // Grid and physical variables
        grid::Grid<T> grid;
        T dt;
        T molecularDiffusion;
        T alphaL;
        T alphaT;
        // Raw pointer to curand state vector
        curandState_t* states;
        // If true, use trilinear interpolation
        bool useTrilinearCorrection;

        /**
        * @brief Constructor.
        * @param _grid Grid where the velocity is defined
        */
        MoveParticle(const grid::Grid<T>& _grid) : grid(_grid), dt(0.0) {};

        /**
        * @brief Initialize the functor.
        * @param _datax Velocity vector (x-direction)
        * @param _datay Velocity vector (y-direction)
        * @param _dataz Velocity vector (z-direction)
        * @param _molecularDiffusion Effective molecular diffusion
        * @param _alphaL Longitudinal dispersivity
        * @param _alphaT Transverse dispersivity
        * @param _states Curand states
        * @param _useTrilinearCorrection True if trilinear correction is used
        */
        void initialize(thrust::device_vector<T> &_datax,
                        thrust::device_vector<T> &_datay,
                        thrust::device_vector<T> &_dataz,
                        thrust::device_vector<T> &_cdatax,
                        thrust::device_vector<T> &_cdatay,
                        thrust::device_vector<T> &_cdataz,
                        T _molecularDiffusion,
                        T _alphaL,
                        T _alphaT,
                        thrust::device_vector<curandState_t> &_states,
                        bool _useTrilinearCorrection)
        {
            datax = thrust::raw_pointer_cast(_datax.data());
            datay = thrust::raw_pointer_cast(_datay.data());
            dataz = thrust::raw_pointer_cast(_dataz.data());

            cdatax = thrust::raw_pointer_cast(_cdatax.data());
            cdatay = thrust::raw_pointer_cast(_cdatay.data());
            cdataz = thrust::raw_pointer_cast(_cdataz.data());

            molecularDiffusion = _molecularDiffusion;
            alphaL = _alphaL;
            alphaT = _alphaT;
            states = thrust::raw_pointer_cast(_states.data());

            useTrilinearCorrection = _useTrilinearCorrection;
        }

        /**
        * @brief Set time step.
        * @param _dt Time step
        */
        void setTimeStep(T _dt)
        {
            dt = _dt;
        }

        using Position = thrust::tuple<T, T, T, unsigned int>;

        /**
        * @brief Execute one step of the particle tracking method
        *        on one particle.
        * @param p Initial position of the particle
        * @return Final position of the particle
        */
        __device__
        Position operator()(Position p) const
        {
            int idx, idy, idz;
            grid::idPoint(grid,
                          thrust::get<0>(p),
                          thrust::get<1>(p),
                          thrust::get<2>(p),
                          &idx, &idy, &idz);

            int id = grid::mergeId(grid, idx, idy, idz);

            bool idValid = grid::validId<T>(grid, idx, idy, idz);

            // Velocity term (linear interpolation)
            T vlx, vly, vlz;
            facefield::in<T>
                        (datax, datay, dataz, grid, idx, idy, idz,
                         idValid,
                         thrust::get<0>(p),
                         thrust::get<1>(p),
                         thrust::get<2>(p),
                         &vlx, &vly, &vlz);

            T vcx, vcy, vcz;
            if (useTrilinearCorrection)
            {
                // Velocity correction div(D) (trilinear interpolation)
                cornerfield::velocityCorrection<T>
                            (cdatax, cdatay, cdataz, grid, idx, idy, idz,
                             idValid,
                             thrust::get<0>(p),
                             thrust::get<1>(p),
                             thrust::get<2>(p),
                             molecularDiffusion, alphaL, alphaT,
                             &vcx, &vcy, &vcz);
            }
            else
            {
                // Velocity correction div(D) (block-centered finite difference)
                vcx = idValid ? cdatax[id] : 0;
                vcy = idValid ? cdatay[id] : 0;
                vcz = idValid ? cdataz[id] : 0;
            }

            // Displacement Matrix (trilinear interpolation)
            T B00, B11, B22, B01, B02, B12;
            cornerfield::displacementMatrix<T>
                        (cdatax, cdatay, cdataz, grid, idx, idy, idz,
                         idValid,
                         thrust::get<0>(p),
                         thrust::get<1>(p),
                         thrust::get<2>(p),
                         molecularDiffusion, alphaL, alphaT, dt,
                         &B00, &B11, &B22, &B01, &B02, &B12);

            // Random Displacement
            T xi0 = curand_normal_double(&states[thrust::get<3>(p)]);
            T xi1 = curand_normal_double(&states[thrust::get<3>(p)]);
            T xi2 = curand_normal_double(&states[thrust::get<3>(p)]);

            // Update positions
            T dpx, dpy, dpz;

            dpx = (idValid ? ((vlx + vcx)*dt + (B00*xi0 + B01*xi1 + B02*xi2)) : 0);
            dpy = (idValid ? ((vly + vcy)*dt + (B01*xi0 + B11*xi1 + B12*xi2)) : 0);
            if (grid.nz != 1)
                dpz = (idValid ? ((vlz + vcz)*dt + (B02*xi0 + B12*xi1 + B22*xi2)) : 0);

            // Closed condition: particles cannot exit the domain
            thrust::get<0>(p) += (grid::validX(grid, thrust::get<0>(p) + dpx) ? dpx : 0);
            thrust::get<1>(p) += (grid::validY(grid, thrust::get<1>(p) + dpy) ? dpy : 0);
            if (grid.nz != 1)
                thrust::get<2>(p) += (grid::validZ(grid, thrust::get<2>(p) + dpz) ? dpz : 0);

            return p;
        }
    };
}

#endif //PAR2_MOVEPARTICLE_CUH
