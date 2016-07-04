/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	Header-file for tutmam_alpha_d_iteration.c
*	Do not change this.
*/
#include  "udf.h"
#include  "tutmam_settings.h"

#ifndef _TUTMAM_ALPHA_D_ITERATION_H
#define _TUTMAM_ALPHA_D_ITERATION_H

	/* Substitutes 2x1 vector from another */
	void substVecFromVec(real *out, real *in);

	/* Substitutes 2x1 vector from scalars */
	void substVecFromScal(real *out, real in0, real in1);

	/* Returns 1 if the solution is not converged */
	int isNotConverged(real *AAndBNow, real *AAndB, real convCriterion);
	
	/* Transponates a 2x2 matrix */
	void transponateMatrix(real *out, real *in);
	
	/* Multiplies two 2x2 matrices */
	void multiplyMatrices(real *out, real *in1, real *in2);
	
	/* Multiplies 2x2 matrix with 2x1 vector */
	void multiplyMatrixAndVec(real *out, real *mat, real *vec);
	
	/* Adding lambda to the diagonal of 2x2 matrix */
	void addLambdaToMatrix(real *out, real *mat, real l);
	
	/* Calculating inverse matrix of 2x2 matrix */
	void inverseMatrix(real *out, real *m);
	
	/* Checking if alpha and d are within limits */
	int isAlphaAndDWithinLimits(real *alphaAndD, real *alphaAndDMinInIterations, real *alphaAndDaxInIterations);

	/* Finds alpha and d using the Levenberg-Marquardt iteration algorithm */
	void iterateAlphaAndD(real *alphaAndD, real *AAndB, real *alphaAndDGuess);

#endif
