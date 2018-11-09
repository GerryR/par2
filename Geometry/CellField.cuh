/**
* @file CellField.cuh
* @brief Header file for cellfield.
*        A cellfield is a field that is defined at the center of
*        every cells of the grid.
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

#ifndef PAR2_CELLFIELD_CUH
#define PAR2_CELLFIELD_CUH

#include <fstream>
#include <iostream>
#include <algorithm>

#include "Vector.cuh"
#include "CartesianGrid.cuh"
#include "FaceField.cuh"

namespace par2
{

    namespace cellfield
    {
        /**
        * @brief Build the vector containing a cellfield.
        * @param g Grid where the field is defined
        * @param data Vector that contains the cellfield values
        * @param v Init value for data
        * @tparam T Float number precision
        * @tparam Vector Container for data vectors
        */
        template<typename T, class Vector>
        void build(const grid::Grid<T>& g, Vector& data, T v = 0)
        {
            int size = g.nx*g.ny*g.nz;
            data.resize(size, v);
        }

        /**
        * @brief Compute drift correction in each cell.
        * @param g Grid where the field is defined
        * @param datax Vector that contains the values on the interfaces
        *        orthogonal to the x-axis
        * @param datay Vector that contains the values on the interfaces
        *        orthogonal to the y-axis
        * @param dataz Vector that contains the values on the interfaces
        *        orthogonal to the z-axis
        * @param driftx Vector that contains the x-values of the correction
        * @param drifty Vector that contains the y-values of the correction
        * @param driftz Vector that contains the z-values of the correction
        * @param Dm Effective molecular diffusion
        * @param alphaL Longitudinal dispersivity
        * @param alphaT Transverse dispersivity
        * @tparam T Float number precision
        * @tparam Vector Container for data vectors
        */
        template<typename T, class Vector>
        void computeDriftCorrection(const grid::Grid<T>& g,
                          const Vector& datax,
                          const Vector& datay,
                          const Vector& dataz,
                          Vector& driftx,
                          Vector& drifty,
                          Vector& driftz,
                          T Dm,
                          T alphaL,
                          T alphaT)
        {
            Vector D11, D22, D33;
            Vector D12, D13, D23;
            build(g, D11);
            build(g, D22);
            build(g, D33);
            build(g, D12);
            build(g, D13);
            build(g, D23);

            // Compute tensor D in all the cells
            for (auto idz = 0; idz < g.nz; idz++)
            {
                for (auto idy = 0; idy < g.ny; idy++)
                {
                    for (auto idx = 0; idx < g.nx; idx++)
                    {
                        T cx, cy, cz;
                        grid::centerOfCell<T>(g, idx, idy, idz, &cx, &cy, &cz);

                        T vx, vy, vz;
                        par2::facefield::in<T>(datax.data(), datay.data(), dataz.data(),
                                            g, idx, idy, idz, true, cx, cy, cz, &vx, &vy, &vz);

                        int id = grid::mergeId(g, idx, idy, idz);
                        T vnorm = par2::vector::norm(vx, vy, vz);

                        D11[id] = (alphaT*vnorm + Dm) + (alphaL - alphaT)*vx*vx/vnorm;
                        D22[id] = (alphaT*vnorm + Dm) + (alphaL - alphaT)*vy*vy/vnorm;
                        D33[id] = (alphaT*vnorm + Dm) + (alphaL - alphaT)*vz*vz/vnorm;
                        D12[id] = (alphaL - alphaT)*vx*vy/vnorm;
                        D13[id] = (alphaL - alphaT)*vx*vz/vnorm;
                        D23[id] = (alphaL - alphaT)*vy*vz/vnorm;
                    }
                }
            }



            // Compute drift correction term
            for (auto idz = 0; idz < g.nz; idz++)
            {
                for (auto idy = 0; idy < g.ny; idy++)
                {
                    for (auto idx = 0; idx < g.nx; idx++)
                    {
                        int id1, id2;

                        // Derivatives in x-direction
                        T dD11x, dD12x, dD13x;

                        T ddx;
                        int idx1, idx2;
                        if (idx == 0)
                        {
                            ddx = g.dx;
                            idx1 = idx;
                            idx2 = idx+1;
                        }
                        else if (idx == g.nx-1)
                        {
                            ddx = g.dx;
                            idx1 = idx-1;
                            idx2 = idx;
                        }
                        else
                        {
                            ddx = 2*g.dx;
                            idx1 = idx-1;
                            idx2 = idx+1;
                        }

                        id1 = grid::mergeId(g, idx1, idy, idz);
                        id2 = grid::mergeId(g, idx2, idy, idz);

                        dD11x = (D11[id2] - D11[id1])/ddx;
                        dD12x = (D12[id2] - D12[id1])/ddx;
                        dD13x = (D13[id2] - D13[id1])/ddx;

                        // Derivatives in y-direction
                        T dD12y, dD22y, dD23y;

                        T ddy;
                        int idy1, idy2;
                        if (idy == 0)
                        {
                            ddy = g.dy;
                            idy1 = idy;
                            idy2 = idy+1;
                        }
                        else if (idy == g.ny-1)
                        {
                            ddy = g.dy;
                            idy1 = idy-1;
                            idy2 = idy;
                        }
                        else
                        {
                            ddy = 2*g.dy;
                            idy1 = idy-1;
                            idy2 = idy+1;
                        }

                        id1 = grid::mergeId(g, idx, idy1, idz);
                        id2 = grid::mergeId(g, idx, idy2, idz);

                        dD12y = (D12[id2] - D12[id1])/ddy;
                        dD22y = (D22[id2] - D22[id1])/ddy;
                        dD23y = (D23[id2] - D23[id1])/ddy;

                        // Derivatives in y-direction
                        T dD13z, dD23z, dD33z;

                        // TODO fix 2D case
                        T ddz;
                        int idz1, idz2;
                        if (idz == 0)
                        {
                            ddz = g.dz;
                            idz1 = idz;
                            idz2 = idz+1;
                        }
                        else if (idz == g.nz-1)
                        {
                            ddz = g.dz;
                            idz1 = idz-1;
                            idz2 = idz;
                        }
                        else
                        {
                            ddz = 2*g.dz;
                            idz1 = idz-1;
                            idz2 = idz+1;
                        }

                        id1 = grid::mergeId(g, idx, idy, idz1);
                        id2 = grid::mergeId(g, idx, idy, idz2);

                        dD13z = (D13[id2] - D13[id1])/ddz;
                        dD23z = (D23[id2] - D23[id1])/ddz;
                        dD33z = (D33[id2] - D33[id1])/ddz;

                        // Compute drift coefficient in the cell
                        int id = grid::mergeId(g, idx, idy, idz);

                        driftx[id] = dD11x + dD12y + dD13z;
                        drifty[id] = dD12x + dD22y + dD23z;
                        driftz[id] = dD13x + dD23y + dD33z;
                    }
                }
            }

        }

    };
};

#endif //PAR2_CELLFIELD_CUH
