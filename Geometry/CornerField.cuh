/**
* @file CornerField.cuh
* @brief Header file for cornerfield.
*        A cornerfield is a field that is defined at the each corner
*        of every cells of the grid.
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

#ifndef PAR2_CORNERFIELD_CUH
#define PAR2_CORNERFIELD_CUH

#include <fstream>
#include <iostream>

#include "Vector.cuh"
#include "CartesianGrid.cuh"
#include "FaceField.cuh"

namespace par2
{

    namespace cornerfield
    {
        /**
        * @brief Build the vector containing a cornerfield.
        * @param g Grid where the field is defined
        * @param data Vector that contains the cornerfield values
        * @param v Init value for data
        * @tparam T Float number precision
        * @tparam Vector Container for data vectors
        */
        template<typename T, class Vector>
        void build(const grid::Grid<T>& g, Vector& data, T v = 0)
        {
            int size = (g.nx+1)*(g.ny+1)*(g.nz+1);
            data.resize(size, v);
        }

        /**
        * @brief Get the position on the data vector for the corner
        *        defined by the IDs in each direction.
        * @param g Grid where the field is defined
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        * @return The position on the data vector
        */
        template<typename T>
        __host__ __device__
        int mergeId(const grid::Grid<T>& g, int idx, int idy, int idz)
        {
            return idz*(g.ny+1)*(g.nx+1) + idy*(g.nx+1) + idx;
        }

        /**
        * @brief Identification for each corner. Each cell has
        *        8 corners.
        */
        enum Direction
        {
            C000 = 0, // CXYZ
            C100,
            C010,
            C110,
            C001,
            C101,
            C011,
            C111
        };

        /**
        * @brief Get the value of the cornerfield at given corner.
        *        The location is found from the IDs of the cell and the
        *        direction of the corner.
        * @param data Vector that contains the cornerfield values
        * @param g Grid where the field is defined
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        * @tparam dir Direction of the interface
        * @return The value of the cornerfield at the corner
        */
        template<typename T, int dir>
        __host__ __device__
        T get(const T* data, const grid::Grid<T>& g,
              int idx, int idy, int idz)
        {
            int id;
            switch (dir)
            {
                case Direction::C000:
                    id = cornerfield::mergeId(g, idx, idy, idz);
                    return data[id];
                case Direction::C100:
                    id = cornerfield::mergeId(g, idx+1, idy, idz);
                    return data[id];
                case Direction::C010:
                    id = cornerfield::mergeId(g, idx, idy+1, idz);
                    return data[id];
                case Direction::C110:
                    id = cornerfield::mergeId(g, idx+1, idy+1, idz);
                    return data[id];
                case Direction::C001:
                    id = cornerfield::mergeId(g, idx, idy, idz+1);
                    return data[id];
                case Direction::C101:
                    id = cornerfield::mergeId(g, idx+1, idy, idz+1);
                    return data[id];
                case Direction::C011:
                    id = cornerfield::mergeId(g, idx, idy+1, idz+1);
                    return data[id];
                case Direction::C111:
                    id = cornerfield::mergeId(g, idx+1, idy+1, idz+1);
                    return data[id];

            }
            return 0;
        }

        /**
        * @brief Get the value of the field at given point using trilinear
        *        interpolation.
        * @param cdatax Vector that contains the x-values of the field
        * @param cdatay Vector that contains the y-values of the field
        * @param cdataz Vector that contains the z-values of the field
        * @param g Grid where the field is defined
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @param idValid True if the position is inside the grid
        * @param px Position x-coordinate
        * @param py Position y-coordinate
        * @param pz Position z-coordinate
        * @param vx Result x-coordinate
        * @param vy Result y-coordinate
        * @param vz Result z-coordinate
        * @tparam T Float number precision
        */
        template<typename T>
        __host__ __device__
        void in(const T* cdatax, const T* cdatay, const T* cdataz, const grid::Grid<T>& g,
                int idx, int idy, int idz, bool idValid, T px, T py, T pz, T* vx, T* vy, T* vz)
        {
            using Direction = cornerfield::Direction;

            T cx, cy, cz;
            grid::centerOfCell<T>(g, idx, idy, idz, &cx, &cy, &cz);

            T x = (px - cx)/g.dx + 0.5;
            T y = (py - cy)/g.dy + 0.5;
            T z = (pz - cz)/g.dz + 0.5;

            *vx = idValid ? interpolation::trilinear(x, y, z,
                                get<T, Direction::C000>(cdatax, g, idx, idy, idz),
                                get<T, Direction::C100>(cdatax, g, idx, idy, idz),
                                get<T, Direction::C010>(cdatax, g, idx, idy, idz),
                                get<T, Direction::C110>(cdatax, g, idx, idy, idz),
                                get<T, Direction::C001>(cdatax, g, idx, idy, idz),
                                get<T, Direction::C101>(cdatax, g, idx, idy, idz),
                                get<T, Direction::C011>(cdatax, g, idx, idy, idz),
                                get<T, Direction::C111>(cdatax, g, idx, idy, idz)) : 0;

            *vy = idValid ? interpolation::trilinear(x, y, z,
                                get<T, Direction::C000>(cdatay, g, idx, idy, idz),
                                get<T, Direction::C100>(cdatay, g, idx, idy, idz),
                                get<T, Direction::C010>(cdatay, g, idx, idy, idz),
                                get<T, Direction::C110>(cdatay, g, idx, idy, idz),
                                get<T, Direction::C001>(cdatay, g, idx, idy, idz),
                                get<T, Direction::C101>(cdatay, g, idx, idy, idz),
                                get<T, Direction::C011>(cdatay, g, idx, idy, idz),
                                get<T, Direction::C111>(cdatay, g, idx, idy, idz)) : 0;

            *vz = idValid ? interpolation::trilinear(x, y, z,
                                get<T, Direction::C000>(cdataz, g, idx, idy, idz),
                                get<T, Direction::C100>(cdataz, g, idx, idy, idz),
                                get<T, Direction::C010>(cdataz, g, idx, idy, idz),
                                get<T, Direction::C110>(cdataz, g, idx, idy, idz),
                                get<T, Direction::C001>(cdataz, g, idx, idy, idz),
                                get<T, Direction::C101>(cdataz, g, idx, idy, idz),
                                get<T, Direction::C011>(cdataz, g, idx, idy, idz),
                                get<T, Direction::C111>(cdataz, g, idx, idy, idz)) : 0;
        }

        /**
        * @brief Compute velocity correction for the particle tracking.
        * @param cdatax Vector that contains the x-values of the field
        * @param cdatay Vector that contains the y-values of the field
        * @param cdataz Vector that contains the z-values of the field
        * @param g Grid where the field is defined
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @param idValid True if the position is inside the grid
        * @param px Position x-coordinate
        * @param py Position y-coordinate
        * @param pz Position z-coordinate
        * @param Dm Effective molecular diffusion
        * @param alphaL Longitudinal dispersivity
        * @param alphaT Transverse dispersivity
        * @param vx Result x-coordinate
        * @param vy Result y-coordinate
        * @param vz Result z-coordinate
        * @tparam T Float number precision
        */
        template<typename T>
        __host__ __device__
        void velocityCorrection(const T* cdatax, const T* cdatay, const T* cdataz, const grid::Grid<T>& g,
                int idx, int idy, int idz,
                bool idValid,
                T px, T py, T pz,
                T Dm, T alphaL, T alphaT,
                T* vx, T* vy, T* vz)
        {
            using Direction = cornerfield::Direction;

            T cx, cy, cz;
            grid::centerOfCell<T>(g, idx, idy, idz, &cx, &cy, &cz);

            T x = (px - cx)/g.dx + 0.5;
            T y = (py - cy)/g.dy + 0.5;
            T z = (pz - cz)/g.dz + 0.5;

            T vx000 = idValid ? get<T, Direction::C000>(cdatax, g, idx, idy, idz) : 1;
            T vx100 = idValid ? get<T, Direction::C100>(cdatax, g, idx, idy, idz) : 1;
            T vx010 = idValid ? get<T, Direction::C010>(cdatax, g, idx, idy, idz) : 1;
            T vx110 = idValid ? get<T, Direction::C110>(cdatax, g, idx, idy, idz) : 1;
            T vx001 = idValid ? get<T, Direction::C001>(cdatax, g, idx, idy, idz) : 1;
            T vx101 = idValid ? get<T, Direction::C101>(cdatax, g, idx, idy, idz) : 1;
            T vx011 = idValid ? get<T, Direction::C011>(cdatax, g, idx, idy, idz) : 1;
            T vx111 = idValid ? get<T, Direction::C111>(cdatax, g, idx, idy, idz) : 1;

            T vy000 = idValid ? get<T, Direction::C000>(cdatay, g, idx, idy, idz) : 1;
            T vy100 = idValid ? get<T, Direction::C100>(cdatay, g, idx, idy, idz) : 1;
            T vy010 = idValid ? get<T, Direction::C010>(cdatay, g, idx, idy, idz) : 1;
            T vy110 = idValid ? get<T, Direction::C110>(cdatay, g, idx, idy, idz) : 1;
            T vy001 = idValid ? get<T, Direction::C001>(cdatay, g, idx, idy, idz) : 1;
            T vy101 = idValid ? get<T, Direction::C101>(cdatay, g, idx, idy, idz) : 1;
            T vy011 = idValid ? get<T, Direction::C011>(cdatay, g, idx, idy, idz) : 1;
            T vy111 = idValid ? get<T, Direction::C111>(cdatay, g, idx, idy, idz) : 1;

            T vz000 = idValid ? get<T, Direction::C000>(cdataz, g, idx, idy, idz) : 1;
            T vz100 = idValid ? get<T, Direction::C100>(cdataz, g, idx, idy, idz) : 1;
            T vz010 = idValid ? get<T, Direction::C010>(cdataz, g, idx, idy, idz) : 1;
            T vz110 = idValid ? get<T, Direction::C110>(cdataz, g, idx, idy, idz) : 1;
            T vz001 = idValid ? get<T, Direction::C001>(cdataz, g, idx, idy, idz) : 1;
            T vz101 = idValid ? get<T, Direction::C101>(cdataz, g, idx, idy, idz) : 1;
            T vz011 = idValid ? get<T, Direction::C011>(cdataz, g, idx, idy, idz) : 1;
            T vz111 = idValid ? get<T, Direction::C111>(cdataz, g, idx, idy, idz) : 1;

            // If velocity is zero in one corner, add eps to vx component to avoid nan
            const T toll = 0.01*Dm/alphaL;

            vx000 = (vx000 < toll && vy000 < toll && vz000 < toll) ? toll : vx000;
            vx100 = (vx100 < toll && vy100 < toll && vz100 < toll) ? toll : vx100;
            vx010 = (vx010 < toll && vy010 < toll && vz010 < toll) ? toll : vx010;
            vx110 = (vx110 < toll && vy110 < toll && vz110 < toll) ? toll : vx110;
            vx001 = (vx001 < toll && vy001 < toll && vz001 < toll) ? toll : vx001;
            vx101 = (vx101 < toll && vy101 < toll && vz101 < toll) ? toll : vx101;
            vx011 = (vx011 < toll && vy011 < toll && vz011 < toll) ? toll : vx011;
            vx111 = (vx111 < toll && vy111 < toll && vz111 < toll) ? toll : vx111;

            T vnorm000 = par2::vector::norm(vx000, vy000, vz000);
            T vnorm100 = par2::vector::norm(vx100, vy100, vz100);
            T vnorm010 = par2::vector::norm(vx010, vy010, vz010);
            T vnorm110 = par2::vector::norm(vx110, vy110, vz110);
            T vnorm001 = par2::vector::norm(vx001, vy001, vz001);
            T vnorm101 = par2::vector::norm(vx101, vy101, vz101);
            T vnorm011 = par2::vector::norm(vx011, vy011, vz011);
            T vnorm111 = par2::vector::norm(vx111, vy111, vz111);

            T dDxxx = interpolation::trilinearDevX(x, y, z, g.dx,
                        (alphaT*vnorm000 + Dm) + (alphaL - alphaT)*vx000*vx000/vnorm000,
                        (alphaT*vnorm100 + Dm) + (alphaL - alphaT)*vx100*vx100/vnorm100,
                        (alphaT*vnorm010 + Dm) + (alphaL - alphaT)*vx010*vx010/vnorm010,
                        (alphaT*vnorm110 + Dm) + (alphaL - alphaT)*vx110*vx110/vnorm110,
                        (alphaT*vnorm001 + Dm) + (alphaL - alphaT)*vx001*vx001/vnorm001,
                        (alphaT*vnorm101 + Dm) + (alphaL - alphaT)*vx101*vx101/vnorm101,
                        (alphaT*vnorm011 + Dm) + (alphaL - alphaT)*vx011*vx011/vnorm011,
                        (alphaT*vnorm111 + Dm) + (alphaL - alphaT)*vx111*vx111/vnorm111);

            T dDyyy = interpolation::trilinearDevY(x, y, z, g.dy,
                        (alphaT*vnorm000 + Dm) + (alphaL - alphaT)*vy000*vy000/vnorm000,
                        (alphaT*vnorm100 + Dm) + (alphaL - alphaT)*vy100*vy100/vnorm100,
                        (alphaT*vnorm010 + Dm) + (alphaL - alphaT)*vy010*vy010/vnorm010,
                        (alphaT*vnorm110 + Dm) + (alphaL - alphaT)*vy110*vy110/vnorm110,
                        (alphaT*vnorm001 + Dm) + (alphaL - alphaT)*vy001*vy001/vnorm001,
                        (alphaT*vnorm101 + Dm) + (alphaL - alphaT)*vy101*vy101/vnorm101,
                        (alphaT*vnorm011 + Dm) + (alphaL - alphaT)*vy011*vy011/vnorm011,
                        (alphaT*vnorm111 + Dm) + (alphaL - alphaT)*vy111*vy111/vnorm111);

            T dDzzz = interpolation::trilinearDevZ(x, y, z, g.dz,
                        (alphaT*vnorm000 + Dm) + (alphaL - alphaT)*vz000*vz000/vnorm000,
                        (alphaT*vnorm100 + Dm) + (alphaL - alphaT)*vz100*vz100/vnorm100,
                        (alphaT*vnorm010 + Dm) + (alphaL - alphaT)*vz010*vz010/vnorm010,
                        (alphaT*vnorm110 + Dm) + (alphaL - alphaT)*vz110*vz110/vnorm110,
                        (alphaT*vnorm001 + Dm) + (alphaL - alphaT)*vz001*vz001/vnorm001,
                        (alphaT*vnorm101 + Dm) + (alphaL - alphaT)*vz101*vz101/vnorm101,
                        (alphaT*vnorm011 + Dm) + (alphaL - alphaT)*vz011*vz011/vnorm011,
                        (alphaT*vnorm111 + Dm) + (alphaL - alphaT)*vz111*vz111/vnorm111);

            T dDxyx = interpolation::trilinearDevX(x, y, z, g.dx,
                        (alphaL - alphaT)*vx000*vy000/vnorm000,
                        (alphaL - alphaT)*vx100*vy100/vnorm100,
                        (alphaL - alphaT)*vx010*vy010/vnorm010,
                        (alphaL - alphaT)*vx110*vy110/vnorm110,
                        (alphaL - alphaT)*vx001*vy001/vnorm001,
                        (alphaL - alphaT)*vx101*vy101/vnorm101,
                        (alphaL - alphaT)*vx011*vy011/vnorm011,
                        (alphaL - alphaT)*vx111*vy111/vnorm111);

            T dDxyy = interpolation::trilinearDevY(x, y, z, g.dy,
                        (alphaL - alphaT)*vx000*vy000/vnorm000,
                        (alphaL - alphaT)*vx100*vy100/vnorm100,
                        (alphaL - alphaT)*vx010*vy010/vnorm010,
                        (alphaL - alphaT)*vx110*vy110/vnorm110,
                        (alphaL - alphaT)*vx001*vy001/vnorm001,
                        (alphaL - alphaT)*vx101*vy101/vnorm101,
                        (alphaL - alphaT)*vx011*vy011/vnorm011,
                        (alphaL - alphaT)*vx111*vy111/vnorm111);

            T dDxzx = interpolation::trilinearDevX(x, y, z, g.dx,
                        (alphaL - alphaT)*vx000*vz000/vnorm000,
                        (alphaL - alphaT)*vx100*vz100/vnorm100,
                        (alphaL - alphaT)*vx010*vz010/vnorm010,
                        (alphaL - alphaT)*vx110*vz110/vnorm110,
                        (alphaL - alphaT)*vx001*vz001/vnorm001,
                        (alphaL - alphaT)*vx101*vz101/vnorm101,
                        (alphaL - alphaT)*vx011*vz011/vnorm011,
                        (alphaL - alphaT)*vx111*vz111/vnorm111);

            T dDxzz = interpolation::trilinearDevZ(x, y, z, g.dz,
                        (alphaL - alphaT)*vx000*vz000/vnorm000,
                        (alphaL - alphaT)*vx100*vz100/vnorm100,
                        (alphaL - alphaT)*vx010*vz010/vnorm010,
                        (alphaL - alphaT)*vx110*vz110/vnorm110,
                        (alphaL - alphaT)*vx001*vz001/vnorm001,
                        (alphaL - alphaT)*vx101*vz101/vnorm101,
                        (alphaL - alphaT)*vx011*vz011/vnorm011,
                        (alphaL - alphaT)*vx111*vz111/vnorm111);

            T dDyzy = interpolation::trilinearDevY(x, y, z, g.dy,
                        (alphaL - alphaT)*vy000*vz000/vnorm000,
                        (alphaL - alphaT)*vy100*vz100/vnorm100,
                        (alphaL - alphaT)*vy010*vz010/vnorm010,
                        (alphaL - alphaT)*vy110*vz110/vnorm110,
                        (alphaL - alphaT)*vy001*vz001/vnorm001,
                        (alphaL - alphaT)*vy101*vz101/vnorm101,
                        (alphaL - alphaT)*vy011*vz011/vnorm011,
                        (alphaL - alphaT)*vy111*vz111/vnorm111);

            T dDyzz = interpolation::trilinearDevZ(x, y, z, g.dz,
                        (alphaL - alphaT)*vy000*vz000/vnorm000,
                        (alphaL - alphaT)*vy100*vz100/vnorm100,
                        (alphaL - alphaT)*vy010*vz010/vnorm010,
                        (alphaL - alphaT)*vy110*vz110/vnorm110,
                        (alphaL - alphaT)*vy001*vz001/vnorm001,
                        (alphaL - alphaT)*vy101*vz101/vnorm101,
                        (alphaL - alphaT)*vy011*vz011/vnorm011,
                        (alphaL - alphaT)*vy111*vz111/vnorm111);

            *vx = dDxxx + dDxyy + dDxzz;
            *vy = dDxyx + dDyyy + dDyzz;
            *vz = dDxzx + dDyzy + dDzzz;

        }

        /**
        * @brief Compute displacement matrix for the particle tracking.
        * @param cdatax Vector that contains the x-values of the field
        * @param cdatay Vector that contains the y-values of the field
        * @param cdataz Vector that contains the z-values of the field
        * @param g Grid where the field is defined
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @param idValid True if the position is inside the grid
        * @param px Position x-coordinate
        * @param py Position y-coordinate
        * @param pz Position z-coordinate
        * @param Dm Effective molecular diffusion
        * @param alphaL Longitudinal dispersivity
        * @param alphaT Transverse dispersivity
        * @param dt Time step
        * @param B00 Matrix component 00
        * @param B00 Matrix component 11
        * @param B00 Matrix component 22
        * @param B00 Matrix component 01
        * @param B00 Matrix component 02
        * @param B00 Matrix component 12
        * @tparam T Float number precision
        */
        template<typename T>
        __host__ __device__
        void displacementMatrix(const T* cdatax, const T* cdatay, const T* cdataz, const grid::Grid<T>& g,
                int idx, int idy, int idz,
                bool idValid,
                T px, T py, T pz,
                T Dm, T alphaL, T alphaT, T dt,
                T* B00, T* B11, T* B22, T* B01, T* B02, T* B12)
        {
            T vx, vy, vz;
            cornerfield::in<T>(cdatax, cdatay, cdataz, g, idx, idy, idz, idValid,
                               px, py, pz,
                               &vx, &vy, &vz);

            const T toll = 0.01*Dm/alphaL;

            vx = (vx < toll) ? toll : vx;

            T vnorm2 = vector::norm2<T>(vx, vy, vz);
            T vnorm  = sqrt(vnorm2);

            // Compute dispersion Matrix
            T alpha = alphaT*vnorm + Dm;
            T beta  = (alphaL - alphaT)/vnorm; // Check threshold needed

            // Eigenvectors
            // vx0 == vx, etc.

            T vx1 = -vy;
            T vy1 = vx;
            T vz1 = 0.0;

            T vx2 = -vz*vx;
            T vy2 = -vz*vy;
            T vz2 = vx*vx + vy*vy;

            // vnorm2_0 == vnorm2
            T vnorm2_1 = vector::norm2<T>(vx1, vy1, vz1);
            T vnorm2_2 = vector::norm2<T>(vx2, vy2, vz2);

            T gamma0 = sqrt(alpha + beta*vnorm2)/vnorm2;
            T gamma1 = sqrt(alpha)/vnorm2_1;
            T gamma2 = sqrt(alpha)/vnorm2_2;

            // Pre-coefficient
            T coeff = sqrt(2.0*dt);

            // Dispersion Matrix
            *B00 = coeff*(gamma0*vx *vx +
                          gamma1*vx1*vx1 +
                          gamma2*vx2*vx2);

            *B11 = coeff*(gamma0*vy *vy +
                          gamma1*vy1*vy1 +
                          gamma2*vy2*vy2);

            *B22 = coeff*(gamma0*vz *vz +
                          gamma1*vz1*vz1 +
                          gamma2*vz2*vz2);

            *B01 = coeff*(gamma0*vx *vy +
                          gamma1*vx1*vy1 +
                          gamma2*vx2*vy2);

            // B10 == B01

            *B02 = coeff*(gamma0*vx *vz +
                          gamma1*vx1*vz1 +
                          gamma2*vx2*vz2);

            // B20 == B02

            *B12 = coeff*(gamma0*vy *vz +
                          gamma1*vy1*vz1 +
                          gamma2*vy2*vz2);

            // B21 == B12

        }

        /**
        * @brief Compute cornerfield from facefield.
        * @param g Grid where the field is defined
        * @param datax Vector that contains the values on the interfaces
        *        orthogonal to the x-axis
        * @param datay Vector that contains the values on the interfaces
        *        orthogonal to the y-axis
        * @param dataz Vector that contains the values on the interfaces
        *        orthogonal to the z-axis
        * @param cdatax Vector that contains the x-values of the cornerfield
        * @param cdatay Vector that contains the y-values of the cornerfield
        * @param cdataz Vector that contains the z-values of the cornerfield
        * @tparam T Float number precision
        * @tparam Vector Container for data vectors
        */
        template<typename T, class Vector>
        void computeCornerVelocities(
                          const grid::Grid<T>& g,
                          const Vector& datax,
                          const Vector& datay,
                          const Vector& dataz,
                          Vector& cdatax,
                          Vector& cdatay,
                          Vector& cdataz)
        {
            using Direction = grid::Direction;

            auto dataxPtr = datax.data();
            auto datayPtr = datay.data();
            auto datazPtr = dataz.data();

            for (auto idz = 0; idz < g.nz+1; idz++)
            {
                for (auto idy = 0; idy < g.ny+1; idy++)
                {
                    for (auto idx = 0; idx < g.nx+1; idx++)
                    {
                        int id = cornerfield::mergeId(g, idx, idy, idz);

                        T vx, vy, vz;
                        int fx, fy, fz;
                        int tx, ty, tz;

                        // Velocity x-direction
                        vx = 0;
                        fx = 0;
                        tx = 0;

                        if (idx == g.nx)
                        {
                            tx = 1;
                        }

                        for (auto ty = 0; ty <= 1; ty++)
                        {
                            for (auto tz = 0; tz <= 1; tz++)
                            {
                                if (grid::validId(g, idx-tx, idy-ty, idz-tz))
                                {
                                    if (tx == 0)
                                    {
                                        vx += facefield::get<T, Direction::XM>(
                                            dataxPtr, datayPtr, datazPtr, g, idx-tx, idy-ty, idz-tz);
                                    }
                                    else
                                    {
                                        vx += facefield::get<T, Direction::XP>(
                                            dataxPtr, datayPtr, datazPtr, g, idx-tx, idy-ty, idz-tz);
                                    }
                                    fx += 1;
                                }
                            }
                        }

                        cdatax[id] = vx/fx;

                        // Velocity y-direction
                        vy = 0;
                        fy = 0;
                        ty = 0;

                        if (idy == g.ny)
                        {
                            ty = 1;
                        }

                        for (auto tx = 0; tx <= 1; tx++)
                        {
                            for (auto tz = 0; tz <= 1; tz++)
                            {
                                if (grid::validId(g, idx-tx, idy-ty, idz-tz))
                                {
                                    if (ty == 0)
                                    {
                                        vy += facefield::get<T, Direction::YM>(
                                            dataxPtr, datayPtr, datazPtr, g, idx-tx, idy-ty, idz-tz);
                                    }
                                    else
                                    {
                                        vy += facefield::get<T, Direction::YP>(
                                            dataxPtr, datayPtr, datazPtr, g, idx-tx, idy-ty, idz-tz);
                                    }
                                    fy += 1;
                                }
                            }
                        }
                        cdatay[id] = vy/fy;

                        // Velocity z-direction
                        vz = 0;
                        fz = 0;
                        tz = 0;

                        if (idz == g.nz)
                        {
                            tz = 1;
                        }

                        for (auto tx = 0; tx <= 1; tx++)
                        {
                            for (auto ty = 0; ty <= 1; ty++)
                            {
                                if (grid::validId(g, idx-tx, idy-ty, idz-tz))
                                {
                                    if (tz == 0)
                                    {
                                        vz += facefield::get<T, Direction::ZM>(
                                            dataxPtr, datayPtr, datazPtr, g, idx-tx, idy-ty, idz-tz);
                                    }
                                    else
                                    {
                                        vz += facefield::get<T, Direction::ZP>(
                                            dataxPtr, datayPtr, datazPtr, g, idx-tx, idy-ty, idz-tz);
                                    }
                                    fz += 1;
                                }
                            }
                        }
                        cdataz[id] = vz/fz;
                    }
                }
            }
        }
    }

}

#endif //PAR2_CORNERFIELD_CUH
