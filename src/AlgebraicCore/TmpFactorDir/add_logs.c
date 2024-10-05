//   Copyright (c)  1997-2006  John Abbott,  Anna M Bigatti

//   This file is part of the source of CoCoALib, the CoCoA Library.
//
//   CoCoALib is free software: you can redistribute it and/or modify
//   it under the terms of the GNU General Public License as published by
//   the Free Software Foundation, either version 3 of the License, or
//   (at your option) any later version.
//
//   CoCoALib is distributed in the hope that it will be useful,
//   but WITHOUT ANY WARRANTY; without even the implied warranty of
//   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//   GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with CoCoALib.  If not, see <http://www.gnu.org/licenses/>.

#include "add_logs.h"

#include <math.h>

/***************************************************************************/
/* Compute log(1+x) without losing too much precision if x is small.       */
/* Some systems already have this function (perhaps called log1p),         */
/* nonetheless we define it here as it is not present on all systems.      */
/* Implementation is simple rather than especially fast/accurate/...       */

static double log1plus(double x)
{
  if (x > 0.01) return log(1.0+x);
  return x*(1-x*(1/2.0-x*(1/3.0-x*(1/4.0-x/5.0))));
}

/***************************************************************************/
/* Assuming log1=log(x), and log2=log(y), compute log(x+y).                */
/* Replacing log1plus(x) by log(1+x) causes noticeable loss of precision   */
/* only if the larger of log1 and log2 is very small.  Since add_logs is   */
/* not called very often, extreme speed is not essential.                  */

double add_logs(double log1, double log2)
{
  if (log1 < log2 - 40) return log2;
  if (log2 < log1 - 40) return log1;
  if (log1 < log2) return log2 + log1plus(exp(log1-log2));
  return log1 + log1plus(exp(log2-log1));
}
