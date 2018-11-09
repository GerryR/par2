/**
* @file FaceField.cuh
* @brief Header file for facefield.
*        A facefield is a variable that is defined at the center of
*        of each cell interface of the grid.
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

#ifndef PAR2_FACEFIELD_CUH
#define PAR2_FACEFIELD_CUH

#include <fstream>
#include <iostream>
#include <algorithm>

#include "CartesianGrid.cuh"
#include "Interpolation.cuh"

namespace par2
{

    namespace facefield
    {
        /**
        * @brief Build the vectors containing a facefield.
        * @param g Grid where the field is defined
        * @param datax Vector that contains the values on the interfaces
        *        orthogonal to the x-axis
        * @param datay Vector that contains the values on the interfaces
        *        orthogonal to the y-axis
        * @param dataz Vector that contains the values on the interfaces
        *        orthogonal to the z-axis
        * @param vx Init value for datax
        * @param vy Init value for datay
        * @param vz Init value for dataz
        * @tparam T Float number precision
        * @tparam Vector Container for data vectors
        */
        template<typename T, class Vector>
        void build(const grid::Grid<T>& g, Vector& datax, Vector& datay, Vector& dataz,
                               T vx = 0, T vy = 0, T vz = 0)
        {
            int size = (g.nx+1)*(g.ny+1)*(g.nz+1);
            datax.resize(size, vx);
            datay.resize(size, vy);
            dataz.resize(size, vz);
        }

        /**
        * @brief Get the position on the data vectors for the interface
        *        defined by the IDs in each direction.
        * @param g Grid where the field is defined
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        * @return The position on the data vectors
        */
        template<typename T>
        __host__ __device__
        int mergeId(const grid::Grid<T>& g, int idx, int idy, int idz)
        {
            return idz*(g.ny+1)*(g.nx+1) + idy*(g.nx+1) + idx;
        }

        /**
        * @brief Get the value of the facefield at given interface.
        *        The interface is found from the IDs of the cell and the
        *        direction of the interface.
        * @param datax Vector that contains the values on the interfaces
        *        orthogonal to the x-axis
        * @param datay Vector that contains the values on the interfaces
        *        orthogonal to the y-axis
        * @param dataz Vector that contains the values on the interfaces
        *        orthogonal to the z-axis
        * @param g Grid where the field is defined
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        * @tparam dir Direction of the interface
        * @return The value of the facefield at the interface
        */
        template<typename T, int dir>
        __host__ __device__
        T get(const T* datax, const T* datay, const T* dataz, const grid::Grid<T>& g,
              int idx, int idy, int idz)
        {
            int id;
            switch (dir)
            {
                case grid::XP:
                    id = facefield::mergeId(g, idx+1, idy, idz);
                    //id = (idz)*(g.nx+1)*(g.ny+1) + (idy)*(g.nx+1) + (idx+1);
                    return datax[id];
                case grid::XM:
                    id = facefield::mergeId(g, idx  , idy, idz);
                    //id = (idz)*(g.nx+1)*(g.ny+1) + (idy)*(g.nx+1) + (idx  );
                    return datax[id];
                case grid::YP:
                    id = facefield::mergeId(g, idx, idy+1, idz);
                    //id = (idz)*(g.nx+1)*(g.ny+1) + (idy+1)*(g.nx+1) + (idx);
                    return datay[id];
                case grid::YM:
                    id = facefield::mergeId(g, idx, idy  , idz);
                    //id = (idz)*(g.nx+1)*(g.ny+1) + (idy  )*(g.nx+1) + (idx);
                    return datay[id];
                case grid::ZP:
                    id = facefield::mergeId(g, idx, idy, idz+1);
                    //id = (idz+1)*(g.nx+1)*(g.ny+1) + (idy)*(g.nx+1) + (idx);
                    return dataz[id];
                case grid::ZM:
                    id = facefield::mergeId(g, idx, idy, idz  );
                    //id = (idz  )*(g.nx+1)*(g.ny+1) + (idy)*(g.nx+1) + (idx);
                    return dataz[id];
            }
            return 0;
        }

        /**
        * @brief Get the value of the field at given point using linear
        *        interpolation.
        * @param datax Vector that contains the values on the interfaces
        *        orthogonal to the x-axis
        * @param datay Vector that contains the values on the interfaces
        *        orthogonal to the y-axis
        * @param dataz Vector that contains the values on the interfaces
        *        orthogonal to the z-axis
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
        void in(const T* datax, const T* datay, const T* dataz, const grid::Grid<T>& g,
                int idx, int idy, int idz, bool idValid, T px, T py, T pz, T* vx, T* vy, T* vz)
        {
            using Direction = grid::Direction;

            T cx, cy, cz;
            grid::centerOfCell<T>(g, idx, idy, idz, &cx, &cy, &cz);

            T Dx = px - cx;
            T Dy = py - cy;
            T Dz = pz - cz;

            T vxp = idValid ? get<T, Direction::XP>(datax, datay, dataz, g, idx, idy, idz) : 1;
            T vxm = idValid ? get<T, Direction::XM>(datax, datay, dataz, g, idx, idy, idz) : 1;
            T vyp = idValid ? get<T, Direction::YP>(datax, datay, dataz, g, idx, idy, idz) : 1;
            T vym = idValid ? get<T, Direction::YM>(datax, datay, dataz, g, idx, idy, idz) : 1;
            T vzp = idValid ? get<T, Direction::ZP>(datax, datay, dataz, g, idx, idy, idz) : 1;
            T vzm = idValid ? get<T, Direction::ZM>(datax, datay, dataz, g, idx, idy, idz) : 1;

            *vx = interpolation::linear<T>(Dx/g.dx + 0.5, vxm, vxp);
            *vy = interpolation::linear<T>(Dy/g.dy + 0.5, vym, vyp);
            *vz = interpolation::linear<T>(Dz/g.dz + 0.5, vzm, vzp);
        }

        /**
        * @brief Import velocity field from modflow.
        * @param g Grid where the field is defined
        * @param datax Vector that contains the values on the interfaces
        *        orthogonal to the x-axis
        * @param datay Vector that contains the values on the interfaces
        *        orthogonal to the y-axis
        * @param dataz Vector that contains the values on the interfaces
        *        orthogonal to the z-axis
        * @param fileName Path to the modflow file
        * @param rho Porosity
        * @tparam T Float number precision
        * @tparam Vector Container for data vectors
        */
        template<typename T, class Vector>
        void importFromModflow(const grid::Grid<T>& g,
                               Vector& datax, Vector& datay, Vector& dataz,
                               const std::string fileName,
                               const T rho)
        {
            std::ifstream inStream;
            inStream.open(fileName, std::ifstream::in);
            if (inStream.is_open())
            {
                std::string line;

                bool is2D = g.nz == 1;
                T area, val;
                int id;

                std::string prefix;

                prefix = " 'QXX";

                while (line.compare(0, prefix.size(), prefix) != 0)
                {
                    std::getline(inStream, line);
                }

                area = g.dy*g.dz;

                for (auto z = 0; z < g.nz; z++)
                {
                    for (auto y = 0; y < g.ny; y++)
                    {
                        for (auto x = 0; x < g.nx; x++)
                        {
                            id = z*(g.nx+1)*(g.ny+1) + y*(g.nx+1) + (x+1);
                            inStream >> val;
                            datax[id] = val/area/rho;   // velocity x
                        }
                    }
                }

                prefix = " 'QYY";

                while (line.compare(0, prefix.size(), prefix) != 0)
                {
                    std::getline(inStream, line);
                }

                area = g.dx*g.dz;

                for (auto z = 0; z < g.nz; z++)
                {
                    for (auto y = 0; y < g.ny; y++)
                    {
                        for (auto x = 0; x < g.nx; x++)
                        {
                            id = z*(g.nx+1)*(g.ny+1) + (y+1)*(g.nx+1) + x;
                            inStream >> val;
                            datay[id] = val/area/rho;   // velocity y
                        }
                    }
                }

                if (!is2D)
                {
                    prefix = " 'QZZ";

                    while (line.compare(0, prefix.size(), prefix) != 0)
                    {
                        std::getline(inStream, line);
                    }

                    area = g.dx*g.dy;

                    for (auto z = 0; z < g.nz; z++)
                    {
                        for (auto y = 0; y < g.ny; y++)
                        {
                            for (auto x = 0; x < g.nx; x++)
                            {
                                id = (z+1)*(g.nx+1)*(g.ny+1) + y*(g.nx+1) + x;
                                inStream >> val;
                                dataz[id] = val/area/rho;   // velocity z
                            }
                        }
                    }
                }
            }
            else
            {
                throw std::runtime_error(std::string("Could not open file ") + fileName);
            }
            inStream.close();
        }

        /**
        * @brief Export velocity field to VTK (unusued).
        * @param g Grid where the field is defined
        * @param datax Vector that contains the values on the interfaces
        *        orthogonal to the x-axis
        * @param datay Vector that contains the values on the interfaces
        *        orthogonal to the y-axis
        * @param dataz Vector that contains the values on the interfaces
        *        orthogonal to the z-axis
        * @param fileName Path to the vtk file
        * @tparam T Float number precision
        * @tparam Vector Container for data vectors
        */
        template<typename T, class Vector>
        void exportVTK(const grid::Grid<T>& g, Vector& datax, Vector& datay, Vector& dataz, const std::string fileName)
        {
            std::ofstream outStream;
            outStream.open(fileName);
            if (outStream.is_open())
            {
                outStream << "# vtk DataFile Version 2.0" << std::endl;
                outStream << "Velocity Field" << std::endl;
                outStream << "ASCII" << std::endl;
                outStream << "DATASET STRUCTURED_POINTS" << std::endl;
                outStream << "DIMENSIONS " << g.nx << " " << g.ny << " " << g.nz << std::endl;
                outStream << "ORIGIN " << g.dx << " " << g.dy << " " << g.dz << std::endl;
                outStream << "SPACING " << g.dx << " " << g.dy << " " << g.dz << std::endl;
                outStream << "POINT_DATA " << g.nx * g.ny * g.nz << std::endl;

                outStream << std::endl;
                outStream << "SCALARS velocityX double" << std::endl;
                outStream << "LOOKUP_TABLE default" << std::endl;
                for (auto idz = 0; idz < g.nz; idz++)
                {
                    for (auto idy = 0; idy < g.ny; idy++)
                    {
                        for (auto idx = 0; idx < g.nx; idx++)
                        {
                            T px, py, pz;
                            grid::centerOfCell(g, idx, idy, idz, &px, &py, &pz);

                            // Velocity term (linear interpolation)
                            T vlx, vly, vlz;
                            facefield::in<T>
                                        (datax.data(), datay.data(), dataz.data(),
                                         g, idx, idy, idz,
                                         true,
                                         px, py, pz,
                                         &vlx, &vly, &vlz);

                            outStream << vlx << std::endl;
                        }
                    }
                }

                outStream << std::endl;
                outStream << "SCALARS velocityY double" << std::endl;
                outStream << "LOOKUP_TABLE default" << std::endl;
                for (auto idz = 0; idz < g.nz; idz++)
                {
                    for (auto idy = 0; idy < g.ny; idy++)
                    {
                        for (auto idx = 0; idx < g.nx; idx++)
                        {
                            T px, py, pz;
                            grid::centerOfCell(g, idx, idy, idz, &px, &py, &pz);

                            // Velocity term (linear interpolation)
                            T vlx, vly, vlz;
                            facefield::in<T>
                                        (datax.data(), datay.data(), dataz.data(),
                                         g, idx, idy, idz,
                                         true,
                                         px, py, pz,
                                         &vlx, &vly, &vlz);

                            outStream << vly << std::endl;
                        }
                    }
                }

                outStream << std::endl;
                outStream << "SCALARS velocityZ double" << std::endl;
                outStream << "LOOKUP_TABLE default" << std::endl;
                for (auto idz = 0; idz < g.nz; idz++)
                {
                    for (auto idy = 0; idy < g.ny; idy++)
                    {
                        for (auto idx = 0; idx < g.nx; idx++)
                        {
                            T px, py, pz;
                            grid::centerOfCell(g, idx, idy, idz, &px, &py, &pz);

                            // Velocity term (linear interpolation)
                            T vlx, vly, vlz;
                            facefield::in<T>
                                        (datax.data(), datay.data(), dataz.data(),
                                         g, idx, idy, idz,
                                         true,
                                         px, py, pz,
                                         &vlx, &vly, &vlz);

                            outStream << vlz << std::endl;
                        }
                    }
                }
            }
            else
            {
                throw std::runtime_error(std::string("Could not open file ") + fileName);
            }
            outStream.close();
        }

    };
};


#endif //PAR2_FACEFIELD_CUH
