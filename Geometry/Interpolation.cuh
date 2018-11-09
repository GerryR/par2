/**
* @file Point.cuh
* @brief Header file for interpolation.
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

#ifndef PAR2_INTERPOLATION_CUH
#define PAR2_INTERPOLATION_CUH

namespace par2
{

    namespace interpolation
    {
        /**
        * @brief Normalized linear interpolation.
        * @param x Value between 0 and 1
        * @param v0 Value in 0
        * @param v1 Value in 1
        * @return Linear interpolation
        */
        template<typename T>
        __host__ __device__
        T linear(T x, T v0, T v1)
        {
            return v0*(1-x) + v1*x;
        }

        /**
        * @brief Normalized trilinear interpolation.
        * @param x Value between 0 and 1
        * @param y Value between 0 and 1
        * @param z Value between 0 and 1
        * @param v000 Value in 000
        * @param v100 Value in 100
        * @param v010 Value in 010
        * @param v110 Value in 110
        * @param v001 Value in 001
        * @param v101 Value in 101
        * @param v011 Value in 011
        * @param v111 Value in 111
        * @return Trilinear interpolation
        */
        template<typename T>
        __host__ __device__
        T trilinear(T x, T y, T z,
                    T v000, T v100, T v010, T v110,
                    T v001, T v101, T v011, T v111)
        {
            T v00 = linear(x, v000, v100);
            T v01 = linear(x, v001, v101);
            T v10 = linear(x, v010, v110);
            T v11 = linear(x, v011, v111);

            T v0 = linear(y, v00, v10);
            T v1 = linear(y, v01, v11);

            return linear(z, v0, v1);
        }

        /**
        * @brief Normalized trilinear interpolation of x-derivative.
        * @param x Value between 0 and 1
        * @param y Value between 0 and 1
        * @param z Value between 0 and 1
        * @param dx Block size in x-direction
        * @param v000 Value in 000
        * @param v100 Value in 100
        * @param v010 Value in 010
        * @param v110 Value in 110
        * @param v001 Value in 001
        * @param v101 Value in 101
        * @param v011 Value in 011
        * @param v111 Value in 111
        * @return Trilinear interpolation of x-derivative
        */
        template<typename T>
        __host__ __device__
        T trilinearDevX(T x, T y, T z, T dx,
                        T v000, T v100, T v010, T v110,
                        T v001, T v101, T v011, T v111)
        {
            // vYZ
            T v00 = (v100 - v000)/dx;
            T v01 = (v101 - v001)/dx;
            T v10 = (v110 - v010)/dx;
            T v11 = (v111 - v011)/dx;

            // vZ
            T v0 = linear(y, v00, v10);
            T v1 = linear(y, v01, v11);

            return linear(z, v0, v1);
        }

        /**
        * @brief Normalized trilinear interpolation of y-derivative.
        * @param x Value between 0 and 1
        * @param y Value between 0 and 1
        * @param z Value between 0 and 1
        * @param dy Block size in y-direction
        * @param v000 Value in 000
        * @param v100 Value in 100
        * @param v010 Value in 010
        * @param v110 Value in 110
        * @param v001 Value in 001
        * @param v101 Value in 101
        * @param v011 Value in 011
        * @param v111 Value in 111
        * @return Trilinear interpolation of y-derivative
        */
        template<typename T>
        __host__ __device__
        T trilinearDevY(T x, T y, T z, T dy,
                        T v000, T v100, T v010, T v110,
                        T v001, T v101, T v011, T v111)
        {
            // vXZ
            T v00 = (v010 - v000)/dy;
            T v01 = (v011 - v001)/dy;
            T v10 = (v110 - v100)/dy;
            T v11 = (v111 - v101)/dy;

            // vZ
            T v0 = linear(x, v00, v10);
            T v1 = linear(x, v01, v11);

            return linear(z, v0, v1);
        }

        /**
        * @brief Normalized trilinear interpolation of z-derivative.
        * @param x Value between 0 and 1
        * @param y Value between 0 and 1
        * @param z Value between 0 and 1
        * @param dz Block size in z-direction
        * @param v000 Value in 000
        * @param v100 Value in 100
        * @param v010 Value in 010
        * @param v110 Value in 110
        * @param v001 Value in 001
        * @param v101 Value in 101
        * @param v011 Value in 011
        * @param v111 Value in 111
        * @return Trilinear interpolation of z-derivative
        */
        template<typename T>
        __host__ __device__
        T trilinearDevZ(T x, T y, T z, T dz,
                        T v000, T v100, T v010, T v110,
                        T v001, T v101, T v011, T v111)
        {
            // vXY
            T v00 = (v001 - v000)/dz;
            T v01 = (v011 - v010)/dz;
            T v10 = (v101 - v100)/dz;
            T v11 = (v111 - v110)/dz;

            // vY
            T v0 = linear(x, v00, v10);
            T v1 = linear(x, v01, v11);

            return linear(y, v0, v1);
        }

    }
}


#endif //PAR2_INTERPOLATION_CUH
