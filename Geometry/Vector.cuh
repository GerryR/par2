/**
* @file Vector.cuh
* @brief Header file for vector.
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

#ifndef PAR2_VECTOR_CUH
#define PAR2_VECTOR_CUH

namespace par2
{

    namespace vector
    {
        /**
        * @brief Compute the square of the Euclidean norm.
        * @param vx x-component
        * @param vy y-component
        * @param vz z-component
        * @return Square of the Euclidean norm
        */
        template<typename T>
        __host__ __device__
        T norm2(T vx, T vy, T vz)
        {
            return vx*vx + vy*vy + vz*vz;
        }

        /**
        * @brief Compute the Euclidean norm.
        * @param vx x-component
        * @param vy y-component
        * @param vz z-component
        * @return Euclidean norm
        */
        template<typename T>
        __host__ __device__
        T norm(T vx, T vy, T vz)
        {
            return sqrt(norm2(vx, vy, vz));
        }

        /**
        * @brief Dot product between two vectors.
        * @param vx0 x-component of v0
        * @param vy0 y-component of v0
        * @param vz0 z-component of v0
        * @param vx1 x-component of v1
        * @param vy1 y-component of v1
        * @param vz1 z-component of v1
        * @return v0 dot v1
        */
        template<typename T>
        __host__ __device__
        T dot(T vx0, T vy0, T vz0, T vx1, T vy1, T vz1)
        {
            return vx0 * vx1 + vy0 * vy1 + vz0 * vz1;
        }

    }

}

#endif //PAR2_VECTOR_CUH
