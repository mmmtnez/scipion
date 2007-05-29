/***************************************************************************
 *
 * Authors:     Roberto Marabini
 *
 * Unidad de  Bioinformatica of Centro Nacional de Biotecnologia , CSIC
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 *
 *  All comments concerning this program package may be sent to the
 *  e-mail address 'xmipp@cnb.uam.es'
 ***************************************************************************/

#ifndef _REFINEMENT_HH
#define _REFINEMENT_HH

#include <data/image.h>
#include <data/volume.h>
#include <data/funcs.h>
#include <data/matrix2d.h>

#include "projection.h"

/**@name Shift refinement */
//@{

/**Correlates two projections and finds the maximun of the correlation matrix.
*/
void calculate_and_find_correlation_max_proj(Projection const &proj1,
        Projection const &proj2,
        Projection & proj_tmp,
        double &shift_X, double &shift_Y,
        double const max_step,
        int ref_trans_after, int act_proj);

/**Correlates two matrices  and finds the maximun of the correlation matrix.
   This center may not be at  an integer position.
   The routine works as follows:
   \begin{enumerate}
   \item Search for the maximun with pixel acuraccy inside the window
   \item Calculate the gravity centre of the corelation
         in a neighborhood such as maximum/sqrt(2) > value
   \item Look for the gravity centre in this neighborhood
   \end{enumerate}
   {\bf Note:} The neighborhood is circular
*/
template <class T>
void calculate_and_find_correlation_max_mat(Matrix2D<T> const &mat1,
        Matrix2D<T> const &mat2,
        Matrix2D<T> & mat_temp,
        double &shift_X, double &shift_Y,
        double const max_step)
{

    //WRAP is OK for crystal but may be not for single particles

    //calculate correlation matrix
    correlation_matrix(mat1, mat2, mat_temp);
    mat_temp.setXmippOrigin();

    //search for the maximun inside window "window"
    int max_step_int = (int)ROUND(max_step);
    int imax, jmax;

    Matrix2D<int> window(2*max_step_int + 1, 2*max_step_int + 1);
    window.setXmippOrigin();

    T temp_max = -100000000;
    SPEED_UP_temps;//speedup temporal variables

    //******* search for the maximun with pixel acuraccy **************/
    FOR_ALL_ELEMENTS_IN_MATRIX2D(window)
    {
        if (temp_max < mat_temp(i, j))
        {
            temp_max = mat_temp(i, j);
            imax = i;
            jmax = j;
        }
    }
    //       #define DEBUG_calculate_and_find_correlation_max_mat
#ifdef  DEBUG_calculate_and_find_correlation_max_mat
    cout << "\n imax, jmax: " << imax << " " << jmax << endl;
    cout.flush();
#endif

    /******* Calculate the gravity centre of the corelation        ****/
    /******* in a neighborhood such as  maximum/sqrt(2) > value   ****/
    int n_max;
    int stop_loop;
    int i, j, i_actual, j_actual;
    int radius;
    stop_loop = false;
    n_max = 0;
    while (!stop_loop)
    {
        n_max ++;
        radius = n_max * n_max;
        for (i = -n_max; i <= n_max; i++)
            for (j = -n_max; j <= n_max; j++)
            {
                i_actual = i + imax;
                j_actual = j + jmax;
                if ((i*i + j*j) > radius)
                    continue;
                if (!mat_temp.outside(i_actual, j_actual))
                {
                    if (temp_max / 1.414 > mat_temp(i_actual, j_actual))
                        stop_loop = true;
                }
                else
                {
                    stop_loop = true;
                    cout << "\nWarning(calculate_and_find_correlation_max_mat)"
                    << "\n some points neede to determine the maxima"
                    << "\n are not available" << endl;
                }
            }
    }

    //       #define DEBUG_calculate_and_find_correlation_max_mat
#ifdef DEBUG_calculate_and_find_correlation_max_mat
    cout << "\n n_max, temp_max: " << n_max << " " << temp_max << endl;
    cout.flush();
#endif

    /*** We have the neighborhood => looking for the gravity centre ***/

    double jj_max ;
    double ii_max ;
    double sum_corr ;
    //n_max and radius has been calculated above
    jj_max = ii_max = sum_corr = 0.;
    for (i = -n_max; i <= n_max; i++)
        for (j = -n_max; j <= n_max; j++)
        {
            i_actual = i + imax;
            j_actual = j + jmax;
            if ((i*i + j*j) > radius)
                continue;
            if (!mat_temp.outside(i_actual, j_actual))
            {
                ii_max += i_actual * mat_temp(i_actual, j_actual);
                jj_max += j_actual * mat_temp(i_actual, j_actual);
                sum_corr += mat_temp(i_actual, j_actual);
            }
        }
    shift_X = jj_max / sum_corr;      /***This is the gravity centre ***/
    shift_Y = ii_max / sum_corr;

    //       #define DEBUG_calculate_and_find_correlation_max_mat
#ifdef DEBUG_calculate_and_find_correlation_max_mat
    cout << "\n  shift_XX  "   << shift_X  << endl;
    cout << "\n  shift_Y  "   << shift_Y  << endl;
    cout << "\n  jj_max   "   << jj_max   << endl;
    cout << "\n  ii_max   "   << ii_max   << endl;
    cout << "\n  sum_corr "   << sum_corr << endl;
    cout.flush();
#endif

}
#undef DEBUG_calculate_and_find_correlation_max_mat
//@}

#endif
