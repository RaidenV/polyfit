/*-----------------------------------------------------------------------------

Name        polyfit.h

Purpose     Single-header polynomial fitting of maximum order 5.

Notes       Based on work performed by Nate Domin:

            https://github.com/natedomin/polyfit

History     12 Jan 15  ND   Created
            22 Apr 20  AFB  Modified for C++
-----------------------------------------------------------------------------*/

#ifndef POLYFIT_H
#define POLYFIT_H

#include <vector>

/*-----------------------------------------------------------------------------

Name        polyfit

Purpose     Given a series of x and y values, find a polynomial of N <= 5 order
            which sufficiently approximates the points.

Inputs      x               dependent variables
            y               independent variables
            order           desired order of polynomial
            coefs           output coefficients
            origin          True if the origin should be forced as a point

Returns     int             0 if successful, anything else if not

History     12 Jan 15  ND   Created
            22 Apr 20  AFB  Modified for C++, added forced origin
-----------------------------------------------------------------------------*/
int polyfit( std::vector<double> x,
             std::vector<double> y,
             unsigned int               order,
             std::vector<double>&       coefs,
             bool                       origin = false )
{

    if ( x.size() != y.size() )
        return 1;

    // force through origin
    if ( origin )
    {
        double prevVal = x[0];
        std::vector<double>::iterator xit = x.begin();
        std::vector<double>::iterator yit = y.begin();
        int pos = 0;

        for ( xit = x.begin(); xit != x.end(); ++xit )
        {
            if ( *xit == 0 )
            {
                *yit = 0;
                break;
            }
            // if the current iterator value is greater than zero, and the
            // previous value was less than zero, we should insert zero
            if ( *xit > 0 && prevVal < 0 )
            {
                x.insert( xit, 0 );
                y.insert( yit, 0 );
                break;
            }
            prevVal = *xit;
            pos++;
            yit++;
        }
    }

    // Declarations...
    // ----------------------------------
    enum {maxOrder = 5};

    double B[maxOrder+1] = {0.0f};
    double P[((maxOrder+1) * 2)+1] = {0.0f};
    double A[(maxOrder + 1)*2*(maxOrder + 1)] = {0.0f};

    double sx, sy, powx;

    unsigned int ii, jj, kk;

    // Verify initial conditions....
    // ----------------------------------

    unsigned int countOfElements = x.size();

    // This method requires that the countOfElements >
    // (order+1)
    if (countOfElements <= order)
        return -1;

    // This method has imposed an arbitrary bound of
    // order <= maxOrder.  Increase maxOrder if necessary.
    if (order > maxOrder)
        return -1;

    // Begin Code...
    // ----------------------------------

    // Identify the column vector
    for (ii = 0; ii < countOfElements; ii++)
    {
        sx    = x[ii];
        sy    = y[ii];
        powx = 1;

        for (jj = 0; jj < (order + 1); jj++)
        {
            B[jj] = B[jj] + (sy * powx);
            powx  = powx * sx;
        }
    }

    // Initialize the PowX array
    P[0] = countOfElements;

    // Compute the sum of the Powers of X
    for (ii = 0; ii < countOfElements; ii++)
    {
        sx    = x[ii];
        powx = x[ii];

        for (jj = 1; jj < ((2 * (order + 1)) + 1); jj++)
        {
            P[jj] = P[jj] + powx;
            powx  = powx * sx;
        }
    }

    // Initialize the reduction matrix
    //
    for (ii = 0; ii < (order + 1); ii++)
    {
        for (jj = 0; jj < (order + 1); jj++)
        {
            A[(ii * (2 * (order + 1))) + jj] = P[ii+jj];
        }

        A[(ii*(2 * (order + 1))) + (ii + (order + 1))] = 1;
    }

    // Move the Identity matrix portion of the redux matrix
    // to the left side (find the inverse of the left side
    // of the redux matrix
    for (ii = 0; ii < (order + 1); ii++)
    {
        sx = A[(ii * (2 * (order + 1))) + ii];
        if (sx != 0)
        {
            for (kk = 0; kk < (2 * (order + 1)); kk++)
            {
                A[(ii * (2 * (order + 1))) + kk] =
                    A[(ii * (2 * (order + 1))) + kk] / sx;
            }

            for (jj = 0; jj < (order + 1); jj++)
            {
                if ((jj - ii) != 0)
                {
                    sy = A[(jj * (2 * (order + 1))) + ii];
                    for (kk = 0; kk < (2 * (order + 1)); kk++)
                    {
                        A[(jj * (2 * (order + 1))) + kk] =
                            A[(jj * (2 * (order + 1))) + kk] -
                            sy * A[(ii * (2 * (order + 1))) + kk];
                    }
                }
            }
        }
        else
        {
            // Cannot work with singular matrices
            return -1;
        }
    }

    coefs.resize(order + 1);

    // Calculate and Identify the coefficients
    for (ii = 0; ii < (order + 1); ii++)
    {
        for (jj = 0; jj < (order + 1); jj++)
        {
            sx = 0;
            for (kk = 0; kk < (order + 1); kk++)
            {
                sx = sx + (A[(ii * (2 * (order + 1))) + (kk + (order + 1))] *
                    B[kk]);
            }
            coefs[ii] = sx;
        }
    }

    return 0;
}

#endif // POLYFIT_H

