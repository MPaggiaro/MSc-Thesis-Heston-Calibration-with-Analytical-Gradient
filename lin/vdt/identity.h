/*
 * identity.h
 * For testing purposes
 * 
 *  Created on: Jun 23, 2012
 *      Author: Danilo Piparo, Thomas Hauth, Vincenzo Innocente
 */

/* 
 * VDT is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser Public License for more details.
 * 
 * You should have received a copy of the GNU Lesser Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

#ifndef IDENTITY_H_
#define IDENTITY_H_

namespace vdt{

// Identity double precision ---------------------------------------------------
inline double identity(double x){return x;}
inline double fast_identity(double x){return x;}
inline double identity2D (double x, double /*y*/){return x;}


// Identity float precision ---------------------------------------------------
inline float identityf(float x){return x;}
inline float fast_identityf(float x){return x;}
inline float identity2Df (float x, float /*y*/){return x;}
//------------------------------------------------------------------------------

void identityv(const uint32_t size, double const * __restrict iarray, double* __restrict oarray);
void identityfv(const uint32_t size, float const * __restrict iarray, float* __restrict oarray);
void identity2Dv(const uint32_t size, double const * __restrict iarray, double const * __restrict iarray2, double* __restrict oarray);
void identity2Dfv(const uint32_t size, float const * __restrict iarray, float const * __restrict iarray2, float* __restrict oarray);

} //vdt namespace

#endif /* IDENTITY_H_ */
