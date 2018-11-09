/**
* @file PParticles.cu
* @brief Implementation file for PParticles class.
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

#include "../Geometry/FaceField.cuh"
#include "../Geometry/CornerField.cuh"

#include <thrust/execution_policy.h>
#include <thrust/count.h>
#include <thrust/fill.h>
#include <thrust/sort.h>
#include <fstream>
#include <algorithm>

namespace par2
{
    template<typename T>
    struct InitCURAND
    {
        unsigned long long seed;
        curandState_t *states;
        InitCURAND(unsigned long long _seed, thrust::device_vector<curandState_t> &_states)
        {
            seed = _seed;
            states = thrust::raw_pointer_cast(_states.data());
        }

        __device__
        void operator()(unsigned int i)
        {
            curand_init(seed, i, 0, &states[i]);
        }
    };

    template<typename T>
    struct InitVolume
    {
        curandState_t* states;
        T p1x, p1y, p1z;
        T p2x, p2y, p2z;
        InitVolume(thrust::device_vector<curandState_t> &_states,
                    T _p1x, T _p1y, T _p1z,
                    T _p2x, T _p2y, T _p2z)
        {
            states = thrust::raw_pointer_cast(_states.data());
            p1x = _p1x;
            p1y = _p1y;
            p1z = _p1z;
            p2x = _p2x;
            p2y = _p2y;
            p2z = _p2z;
        }

        using Position = thrust::tuple<T, T, T>;

        __device__
        Position operator()(unsigned int i) const
        {
            Position p;

            thrust::get<0>(p) = p1x + (p2x-p1x)*curand_uniform(&states[i]);
            thrust::get<1>(p) = p1y + (p2y-p1y)*curand_uniform(&states[i]);
            thrust::get<2>(p) = p1z + (p2z-p1z)*curand_uniform(&states[i]);

            return p;
        }
    };

    template<typename T>
    PParticles<T>::PParticles(const grid::Grid<T> &_grid,
                           const thrust::host_vector<T> &_datax,
                           const thrust::host_vector<T> &_datay,
                           const thrust::host_vector<T> &_dataz,
                           T _molecularDiffusion,
                           T _alphaL,
                           T _alphaT,
                           unsigned int _nParticles,
                           long int _seed,
                           bool _useTrilinearCorrection)
            : nParticles(_nParticles), molecularDiffusion(_molecularDiffusion),
              alphaL(_alphaL), alphaT(_alphaT), grid(_grid), moveParticle(_grid),
              useTrilinearCorrection(_useTrilinearCorrection)
    {
        cx.resize(nParticles);
        cy.resize(nParticles);
        cz.resize(nParticles);

        datax = _datax;
        datay = _datay;
        dataz = _dataz;

        thrust::host_vector<T> _cdatax, _cdatay, _cdataz;
        if (useTrilinearCorrection)
        {
            par2::cornerfield::build(grid, _cdatax);
            par2::cornerfield::build(grid, _cdatay);
            par2::cornerfield::build(grid, _cdataz);

            par2::cornerfield::computeCornerVelocities(grid, _datax, _datay, _dataz,
                                                    _cdatax, _cdatay, _cdataz);
        }
        else
        {
            par2::cellfield::build(grid, _cdatax);
            par2::cellfield::build(grid, _cdatay);
            par2::cellfield::build(grid, _cdataz);

            par2::cellfield::computeDriftCorrection(grid, _datax, _datay, _dataz,
                                                    _cdatax, _cdatay, _cdataz,
                                                    molecularDiffusion, alphaL, alphaT);
        }
        cdatax = _cdatax;
        cdatay = _cdatay;
        cdataz = _cdataz;

        states.resize(maxParticles);
        thrust::counting_iterator<unsigned int> count(0);
        thrust::for_each(count, count+maxParticles, InitCURAND<T>(_seed, states));
        
        moveParticle.initialize(datax,
                                datay,
                                dataz,
                                cdatax,
                                cdatay,
                                cdataz,
                                molecularDiffusion,
                                alphaL,
                                alphaT,
                                states,
                                useTrilinearCorrection);

        cudaDeviceSynchronize();

    }

    template<typename T>
    unsigned int PParticles<T>::size() const
    {
        return nParticles;
    }

    template<typename T>
    void PParticles<T>::initializeBox(T p1x, T p1y, T p1z,
                                      T p2x, T p2y, T p2z)
    {
        thrust::counting_iterator<unsigned int> count(0);
        auto pBeg = thrust::make_zip_iterator(
            thrust::make_tuple(cx.begin(), cy.begin(), cz.begin()));

        auto functor = InitVolume<T>(states, p1x, p1y, p1z, p2x, p2y, p2z);

        for (auto i = 0; i*maxParticles < nParticles; i++)
        {
            unsigned int kernelSize = maxParticles;
            if (kernelSize > nParticles - i*maxParticles)
            {
                kernelSize = nParticles - i*maxParticles;
            }
            thrust::transform(count,
                              count + kernelSize,
                              pBeg + i*maxParticles,
                              functor);
        }
    }

    template<typename T>
    void PParticles<T>::move(T dt)
    {
        thrust::counting_iterator<unsigned int> count(0);
        moveParticle.setTimeStep(dt);

        for (auto i = 0; i*maxParticles < nParticles; i++)
        {
            unsigned int kernelSize = maxParticles;
            if (kernelSize > nParticles - i*maxParticles)
            {
                kernelSize = nParticles - i*maxParticles;
            }

            auto pBeg = thrust::make_zip_iterator(
                thrust::make_tuple(cx.begin() + i*maxParticles,
                                   cy.begin() + i*maxParticles,
                                   cz.begin() + i*maxParticles,
                                   count));
            //auto pEnd = thrust::make_zip_iterator(
            //    thrust::make_tuple(cx.end(),   cy.end()  , cz.end()  , count+kernelSize));

            thrust::transform(pBeg, pBeg + kernelSize, pBeg, moveParticle);
        }
        cudaDeviceSynchronize();

    }

    template<typename T>
    void PParticles<T>::exportCSV(const std::string &fileName) const
    {
        // Copy to host memory
        thrust::host_vector<T> hx = cx;
        thrust::host_vector<T> hy = cy;
        thrust::host_vector<T> hz = cz;

        std::ofstream outStream;
        outStream.open(fileName);
        if (outStream.is_open())
        {
            outStream << "id,x coord,y coord,z coord" << std::endl;
            for (unsigned int i = 0; i < nParticles; i++)
            {
                outStream << i << "," << hx[i] << "," << hy[i] << "," << hz[i]
                          << std::endl;
            }
        }
        else
        {
            throw std::runtime_error(std::string("Could not open file ") + fileName);
        }
        outStream.close();
    }

    template<typename T>
    struct isInside
    {
        T plane;

        T p1x, p1y, p1z;
        T p2x, p2y, p2z;
        isInside(T _p1x, T _p1y, T _p1z,
                 T _p2x, T _p2y, T _p2z)
        {
            p1x = _p1x;
            p1y = _p1y;
            p1z = _p1z;
            p2x = _p2x;
            p2y = _p2y;
            p2z = _p2z;
        }

        using Position = thrust::tuple<T, T, T>;

        __device__
        bool operator()(Position p) const
        {
            return (p1x <= thrust::get<0>(p) && thrust::get<0>(p) <= p2x) &&
                   (p1y <= thrust::get<1>(p) && thrust::get<1>(p) <= p2y) &&
                   (p1z <= thrust::get<2>(p) && thrust::get<2>(p) <= p2z);
        }
    };

    template<typename T>
    T PParticles<T>::concentrationBox(T p1x, T p1y, T p1z,
                                      T p2x, T p2y, T p2z) const
    {
        auto pBeg = thrust::make_zip_iterator(
            thrust::make_tuple(cx.begin(), cy.begin(), cz.begin()));
        auto pEnd = thrust::make_zip_iterator(
            thrust::make_tuple(cx.end(),   cy.end()  , cz.end()  ));

        return thrust::count_if(pBeg, pEnd,
                    isInside<T>(p1x, p1y, p1z, p2x, p2y, p2z))/T(nParticles);
    }

    template<typename T>
    struct isAfter
    {
        T plane;

        isAfter(T _plane) : plane(_plane) {};

        __device__
        bool operator()(T x)
        {
            return x > plane;
        }
    };

    template<typename T>
    T PParticles<T>::concentrationAfterX(T xplane) const
    {
        return thrust::count_if(cx.begin(), cx.end(),
                    isAfter<T>(xplane))/T(nParticles);
    }

}
