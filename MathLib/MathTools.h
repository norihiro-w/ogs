/*
 * \file MathTools.h
 *
 *  Created on: Jan 13, 2010
 *      Author: TF
 */

#ifndef MATHTOOLS_H_
#define MATHTOOLS_H_

#include <cstddef>
#include <vector>
#include <cmath>
#include <limits>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace MathLib {

/**
 * standard inner product in R^3
 * \param v0 array of type T representing the vector
 * \param v1 array of type T representing the vector
 * \param n the size of the array
 * */
template<class T> inline
double scpr(const T* v0, const T* v1, size_t n)
{
	long double res(0.0);
	for (size_t k(0); k<n; k++)
		res += v0[k] * v1[k];
	return (double) res;
}

#ifdef _OPENMP
double scpr(double const * const v, double const * const w, OPENMP_LOOP_TYPE n);
#endif

/**
 * computes the cross (or vector) product of the 3d vectors u and v
 * the result is given in the vector r
 */
void crossProd (const double u[3], const double v[3], double r[3]);

/**
 * normalize the vector
 */
void normalizeVector(const double* u, size_t n, double* r);
void normalizeVector(double* u, size_t n);

/**
 * calcProjPntToLineAndDists computes the orthogonal projection
 * of a point p to the line described by the points a and b,
 * \f$g(\lambda) = a + \lambda (b - a)\f$,
 * the distance between p and the projected point
 * and the distances between the projected point and the end
 * points a, b of the line
 * \param p the (mesh) point
 * \param a first point of line
 * \param b second point of line
 * \param lambda the projected point described by the line equation above
 * \param d0 distance to the line point a
 * \returns the distance between p and the orthogonal projection of p
 */
double calcProjPntToLineAndDists(const double p[3], const double a[3],
		const double b[3], double &lambda, double &d0);


/** squared dist between double arrays p0 and p1 (size of arrays is 3) */
double sqrDist(const double* p0, const double* p1);

/** linear normalisation of val from [min, max] into [0,1] */
float normalize(float min, float max, float val);

/**
 * Let \f$p_0, p_1, p_2 \in \mathbb{R}^3\f$. The function getAngle
 * computes the angle between the edges \f$(p_0,p_1)\f$ and \f$(p_1,p_2)\f$
 * @param p0 start point of edge 0
 * @param p1 end point of edge 0 and start point of edge 1
 * @param p2 end point of edge 1
 * @return the angle between the edges
 */
double getAngle (const double p0[3], const double p1[3], const double p2[3]);

/**
 * simple power function that takes as a second argument an integer instead of a float
 * @param base basis of the expression
 * @param exp exponent of the expression
 * @return base^exp
 */
template <typename T> inline
T fastpow (T base, size_t exp)
{
	T result (base);
	if (exp == 0) result = static_cast<T>(1);
	for (size_t k(1); k<exp; k++) {
		result *= base;
	}
	return result;
}

/**
 * Error function
 * @param x
 * @return
 */
double erf(double x);

/**
 * Complementary error function
 * @param x
 * @return
 */
double erfc(double x);

} // namespace

#endif /* MATHTOOLS_H_ */
