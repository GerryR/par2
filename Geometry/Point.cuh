/**
* @file Point.cuh
* @brief Header file for point.
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

#ifndef PAR2_POINT_CUH
#define PAR2_POINT_CUH

namespace par2
{

    namespace point
    {
        /**
        * @brief Compute the distance between two points.
        * @param px0 x-component of p0
        * @param py0 y-component of p0
        * @param pz0 z-component of p0
        * @param px1 x-component of p1
        * @param py1 y-component of p1
        * @param pz1 z-component of p1
        * @return Distance between p0 and p1
        */
        template<typename T>
        __host__ __device__
        T distance(T px0, T py0, T pz0, T px1, T py1, T pz1)
        {
            return sqrt((px0 - px1) * (px0 - px1) + (py0 - py1) * (py0 - py1) + (pz0 - pz1) * (pz0 - pz1));
        };

        /**
        * @brief Compute the square of the distance between two points.
        * @param px0 x-component of p0
        * @param py0 y-component of p0
        * @param pz0 z-component of p0
        * @param px1 x-component of p1
        * @param py1 y-component of p1
        * @param pz1 z-component of p1
        * @return Square of the distance between p0 and p1
        */
        template<typename T>
        __host__ __device__
        T distance2(T px0, T py0, T pz0, T px1, T py1, T pz1)
        {
            return (px0 - px1) * (px0 - px1) + (py0 - py1) * (py0 - py1) + (pz0 - pz1) * (pz0 - pz1);
        };

        /**
        * @brief Summation of the coordinates of two points.
        * @param px0 x-component of p0
        * @param py0 y-component of p0
        * @param pz0 z-component of p0
        * @param px1 x-component of p1
        * @param py1 y-component of p1
        * @param pz1 z-component of p1
        * @param rx x-component of p0 + p1
        * @param ry y-component of p0 + p1
        * @param rz z-component of p0 + p1
        */
        template<typename T>
        __host__ __device__
        void plus(T px0, T py0, T pz0, T px1, T py1, T pz1, T* rx, T* ry, T* rz)
        {
            *rx = px0 + px1;
            *ry = py0 + py1;
            *rz = pz0 + pz1;
        };

        /**
        * @brief Subtraction the coordinates of two points.
        * @param px0 x-component of p0
        * @param py0 y-component of p0
        * @param pz0 z-component of p0
        * @param px1 x-component of p1
        * @param py1 y-component of p1
        * @param pz1 z-component of p1
        * @param rx x-component of p0 - p1
        * @param ry y-component of p0 - p1
        * @param rz z-component of p0 - p1
        */
        template<typename T>
        __host__ __device__
        void minus(T px0, T py0, T pz0, T px1, T py1, T pz1, T* rx, T* ry, T* rz)
        {
            *rx = px0 - px1;
            *ry = py0 - py1;
            *rz = pz0 - pz1;
        };

        /**
        * @brief Multiplication between a point and a scalar.
        * @param px0 x-component of p0
        * @param py0 y-component of p0
        * @param pz0 z-component of p0
        * @param val scalar
        * @param rx x-component of val*p0
        * @param ry y-component of val*p0
        * @param rz z-component of val*p0
        */
        template<typename T>
        __host__ __device__
        void multiply(T px0, T py0, T pz0, T val, T* rx, T* ry, T* rz)
        {
            *rx = px0*val;
            *ry = py0*val;
            *rz = pz0*val;
        };

        /**
        * @brief Division between a point and a scalar.
        * @param px0 x-component of p0
        * @param py0 y-component of p0
        * @param pz0 z-component of p0
        * @param val scalar
        * @param rx x-component of p0/val
        * @param ry y-component of p0/val
        * @param rz z-component of p0/val
        */
        template<typename T>
        __host__ __device__
        void divide(T px0, T py0, T pz0, T val, T* rx, T* ry, T* rz)
        {
            *rx = px0/val;
            *ry = py0/val;
            *rz = pz0/val;
        };

    };

}

#endif //PAR2_POINT_CUH
