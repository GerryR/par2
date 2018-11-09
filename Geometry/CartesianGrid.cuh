/**
* @file CellVariable.cuh
* @brief Header file for grid.
*        Data stuctures and functions of a Cartesian grid
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

#ifndef PAR2_CARTESIANGRID_CUH
#define PAR2_CARTESIANGRID_CUH

#include <stddef.h>
#include "Point.cuh"

namespace par2
{

    namespace grid
    {
        /**
        * @struct Grid
        * @brief Struct containing the Cartesian grid parameters.
        * @tparam T Float number precision
        */
        template <typename T>
        struct Grid
        {
            int nx, ny, nz;
            T dx, dy, dz;
            T px, py, pz;

            Grid() {};

            __host__ __device__
            Grid(const Grid& g)
                : nx(g.nx), ny(g.ny), nz(g.nz),
                  dx(g.dx), dy(g.dy), dz(g.dz),
                  px(g.px), py(g.py), pz(g.pz) {};
        };

        /**
        * @brief Build the grid.
        * @param _nx Number of cells along x
        * @param _ny Number of cells along y
        * @param _nz Number of cells along z
        * @param _dx Size of cell along x
        * @param _dy Size of cell along y
        * @param _dz Size of cell along z
        * @param _px Origin along x
        * @param _py Origin along y
        * @param _pz Origin along z
        * @tparam T Float number precision
        */
        template <typename T>
        Grid<T> build(int _nx, int _ny, int _nz, T _dx, T _dy, T _dz, T _px=0, T _py=0, T _pz=0)
        {
            Grid<T> g;
            g.nx = _nx; g.ny = _ny; g.nz = _nz;
            g.dx = _dx; g.dy = _dy; g.dz = _dz;
            g.px = _px; g.py = _py; g.pz = _pz;
            return g;
        }

        /**
        * @brief Check if x is inside the grid.
        * @param g Grid parameters
        * @param x Coordinate
        * @tparam T Float number precision
        * @return True if x is inside the grid
        */
        template<typename T>
        __host__ __device__
        bool validX(const Grid<T>& g, const T x)
        {
            return g.px < x && x < g.px + g.dx*g.nx;
        }

        /**
        * @brief Check if y is inside the grid.
        * @param g Grid parameters
        * @param y Coordinate
        * @tparam T Float number precision
        * @return True if y is inside the grid
        */
        template<typename T>
        __host__ __device__
        bool validY(const Grid<T>& g, const T y)
        {
            return g.py < y && y < g.py + g.dy*g.ny;
        }

        /**
        * @brief Check if z is inside the grid.
        * @param g Grid parameters
        * @param z Coordinate
        * @tparam T Float number precision
        * @return True if z is inside the grid
        */
        template<typename T>
        __host__ __device__
        bool validZ(const Grid<T>& g, const T z)
        {
            return g.pz < z && z < g.pz + g.dz*g.nz;
        }

        /**
        * @brief Get the total number of cells.
        * @param g Grid parameters
        * @tparam T Float number precision
        * @return Number of cells
        */
        template<typename T>
        __host__ __device__
        int numberOfCells(const Grid<T>& g)
        {
            return g.nx*g.ny*g.nz;
        }

        /**
        * @brief Get the total number of faces.
        * @param g Grid parameters
        * @tparam T Float number precision
        * @return Number of faces
        */
        template<typename T>
        __host__ __device__
        int numberOfFaces(const Grid<T>& g)
        {
            return g.nx*g.ny*g.nz + g.nx*g.ny + g.ny*g.nz + g.nx*g.nz;
        }

        /**
        * @brief Get the ID of the cell containing at given position
        * @param g Grid parameters
        * @param px Position x-coordinate
        * @param py Position y-coordinate
        * @param pz Position z-coordinate
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        */
        template<typename T>
        __host__ __device__
        void idPoint(const Grid<T>& g, T px, T py, T pz, int* idx, int* idy, int* idz)
        {
            *idx = floor((px-g.px)/g.dx);
            *idy = floor((py-g.py)/g.dy);
            *idz = floor((pz-g.pz)/g.dz);
        }

        /**
        * @brief Check if a given ID is valid (3 Ids)
        * @param g Grid parameters
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        * @return True if the ID is valid
        */
        template<typename T>
        __host__ __device__
        bool validId(const Grid<T>& g, int idx, int idy, int idz)
        {
            return idx >= 0 && idx < g.nx && idy >= 0 && idy < g.ny && idz >= 0 && idz < g.nz;
        }

        /**
        * @brief Check if a given ID is valid (unique Id)
        * @param g Grid parameters
        * @param id Unique ID
        * @tparam T Float number precision
        * @return True if the ID is valid
        */
        template<typename T>
        __host__ __device__
        bool validId(const Grid<T>& g, int id)
        {
            return id >= 0 && id < numberOfCells(g);
        }

        /**
        * @brief Given the unique ID of a cell, find the corresponding
        *        IDs along each principal direction.
        * @param g Grid parameters
        * @param id Unique ID
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        */
        template<typename T>
        __host__ __device__
        void splitId(const Grid<T>& g, int id, int* idx, int* idy, int* idz)
        {
            *idz = id / (g.nx*g.ny);
            id = id % (g.nx*g.ny);
            *idy = id / g.nx;
            id = id % g.nx;
            *idx = id;
        }

        /**
        * @brief Given the IDs along each principal direction,
        *        find the corresponding unique ID of a cell without
        *        performing validation.
        * @param g Grid parameters
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        * @return Unique ID
        */
        template<typename T>
        __host__ __device__
        int mergeId(const Grid<T>& g, int idx, int idy, int idz)
        {
            return idz*g.ny*g.nx + idy*g.nx + idx;
        }

        /**
        * @brief Given the IDs along each principal direction,
        *        find the corresponding unique ID of a cell
        *        performing validation. If the IDs are not valid,
        *        return the total number of cells.
        * @param g Grid parameters
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        * @return Unique ID
        */
        template<typename T>
        __host__ __device__
        int uniqueId(const Grid<T>& g, int idx, int idy, int idz)
        {
            return validId(g, idx, idy, idz) ? mergeId(g, idx, idy, idz) : numberOfCells(g);
        }

        /**
        * @brief Compute the volume of a single cell.
        * @param g Grid parameters
        * @tparam T Float number precision
        * @return Volume of a cell
        */
        template<typename T>
        __host__ __device__
        T volumeCell(const Grid<T>& g)
        {
            return g.dx*g.dy*g.dz;
        }

        /**
        * @brief Compute the volume of the whole grid.
        * @param g Grid parameters
        * @tparam T Float number precision
        * @return Volume of the grid
        */
        template<typename T>
        __host__ __device__
        T volumeGrid(const Grid<T>& g)
        {
            return volumeCell(g) * numberOfCells(g);
        }

        /**
        * @brief Compute the center point of a cell.
        * @param g Grid parameters
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @param px Center x-coordinate
        * @param py Center y-coordinate
        * @param pz Center z-coordinate
        * @tparam T Float number precision
        * @return Volume of the grid
        */
        template<typename T>
        __host__ __device__
        void centerOfCell(const Grid<T>& g, int idx, int idy, int idz, T* px, T* py, T* pz)
        {
            *px = g.px + idx*g.dx + 0.5*g.dx;
            *py = g.py + idy*g.dy + 0.5*g.dy;
            *pz = g.pz + idz*g.dz + 0.5*g.dz;
        }

        /**
        * @brief Direction needed to find the neighbors cells.
        */
        enum Direction
        {
            XP = 0, // XP: cell increasing x
            XM,     // XM: cell decreasing x
            YP,
            YM,
            ZP,
            ZM
        };

        /**
        * @brief Find the ID of a neighbor cell.
        * @param g Grid parameters
        * @param idx ID along x
        * @param idy ID along y
        * @param idz ID along z
        * @tparam T Float number precision
        * @tparam dir Direction of the neighbor
        * @return Volume of the grid
        */
        template<typename T, int dir>
        __host__ __device__
        int idNeighbor(const Grid<T>& g, int idx, int idy, int idz)
        {
            switch (dir)
            {
                case XP:
                    return (idx != g.nx-1) ? mergeId(g, idx+1, idy, idz) : numberOfCells(g);
                case XM:
                    return (idx != 0    ) ? mergeId(g, idx-1, idy, idz) : numberOfCells(g);
                case YP:
                    return (idy != g.ny-1) ? mergeId(g, idx, idy+1, idz) : numberOfCells(g);
                case YM:
                    return (idy != 0    ) ? mergeId(g, idx, idy-1, idz) : numberOfCells(g);
                case ZP:
                    return (idz != g.nz-1) ? mergeId(g, idx, idy, idz+1) : numberOfCells(g);
                case ZM:
                    return (idz != 0    ) ? mergeId(g, idx, idy, idz-1) : numberOfCells(g);
            }
            return numberOfCells(g);
        }
    };

}


#endif //PAR2_CARTESIANGRID_CUH
