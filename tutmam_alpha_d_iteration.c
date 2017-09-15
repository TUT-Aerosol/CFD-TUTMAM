/*
*	TUT Modal Aerosol Model for CFD
*	Miska Olin (miska.olin@tut.fi), Tampere University of Technology
* 
*	This file provides iteration algorithm for alpha and d calculation
*	You should not need to change this
*/
#include  "udf.h"
#include  "tutmam_alpha_d_iteration.h"
#include  "tutmam_settings.h"
#include  "tutmam_func.h"

/* Substitutes 2x1 vector from another */
void substVecFromVec(real *out, real *in) {
	out[0] = in[0];
	out[1] = in[1];
	
	return;
}

/* Substitutes 2x1 vector from scalars */
void substVecFromScal(real *out, real in0, real in1) {
	out[0] = in0;
	out[1] = in1;
	
	return;
}

/* Returns 1 if the solution is not converged */
int isNotConverged(real *AAndBNow, real *AAndB, real convCriterion) {
	if (MAX(fabs(AAndBNow[0]/AAndB[0]-1.0),fabs(AAndBNow[1]/AAndB[1]-1.0)) > convCriterion) {
		return 1;
		
	} else {
		return 0;
	}
}

/* Transponates a 2x2 matrix */
void transponateMatrix(real *out, real *in) {
	out[0] = in[0];
	out[1] = in[2];
	out[2] = in[1];
	out[3] = in[3];
	
	return;
}

/* Multiplies two 2x2 matrices */
void multiplyMatrices(real *out, real *in1, real *in2) {
	out[0] = in1[0]*in2[0]+in1[1]*in2[2];
	out[1] = in1[0]*in2[1]+in1[1]*in2[3];
	out[2] = in1[2]*in2[0]+in1[3]*in2[2];
	out[3] = in1[2]*in2[1]+in1[3]*in2[3];
	
	return;
}

/* Multiplies 2x2 matrix with 2x1 vector */
void multiplyMatrixAndVec(real *out, real *mat, real *vec) {
	out[0] = mat[0]*vec[0]+mat[1]*vec[1];
	out[1] = mat[2]*vec[0]+mat[3]*vec[1];
	
	return;
}

/* Adding lambda to the diagonal of 2x2 matrix */
void addLambdaToMatrix(real *out, real *mat, real l) {
	out[0] = mat[0]+l;
	out[1] = mat[1];
	out[2] = mat[2];
	out[3] = mat[3]+l;
	
	return;
}

/* Calculating inverse matrix of 2x2 matrix */
void inverseMatrix(real *out, real *m) {
	real det = m[0]*m[3]-m[1]*m[2]; /* determinant, must not be zero, as it never should be */
	out[0] = m[3]/det;
	out[1] = -m[1]/det;
	out[2] = -m[2]/det;
	out[3] = m[0]/det;
	
	return;
}

/* Checking if alpha and d are within limits */
int isAlphaAndDWithinLimits(real *alphaAndD, real *alphaAndDMinInIterations, real *alphaAndDMaxInIterations) {
	if (alphaAndD[0] < alphaAndDMinInIterations[0]) {
		return 0;
	}
	
	if (alphaAndD[0] > alphaAndDMaxInIterations[0]) {
		return 0;
	}
	
	if (alphaAndD[1] < alphaAndDMinInIterations[1]) {
		return 0;
	}
	
	if (alphaAndD[1] > alphaAndDMaxInIterations[1]) {
		return 0;
	}
	
	/* checking that alpha and d are not inf or nan */
	if (alphaAndD[0] > -1000.0 && alphaAndD[1] > -1000.0) {
		return 1;
	}
	
	Message("NaN or Inf found in Levenberg-Marquardt iteration algorithm!\n");
	return 0;	
}

/* Finds alpha and d using the Levenberg-Marquardt iteration algorithm */
void iterateAlphaAndD(real *alphaAndD, real *AAndB, real *alphaAndDGuess) {
	real initialAlpha[nTutmamModes] = { initialGSD };	/* initial alpha value */
	real J[4];											/* Jacobian matrix */
	real JT[4];											/* Transpose of J */
	real JTJ[4];										/* JT*J */
	real r[2];											/* distance between the real solution and the iterated one */
	real rNew[2];										/* new value for r */
	real costGradient[2];								/* cost gradient */
	real g[4];											/* g */
	real invG[4];										/* inverse of g */
	real AAndBMax[2];									/* maximum A and B, encountered with max alpha and max d */
	real AAndBNow[2];									/* present value of A and B */
	real AAndBNew[2];									/* new value of A and B */
	real alphaAndDMax[2];								/* maximum alpha and d allowed globally */
	real alphaAndDNow[2];								/* present values for alpha and d */
	real alphaAndDNew[2];								/* new values for alpha and d */
	real alphaAndDMinInIterations[2];					/* minimum values for alpha and d allowed during iterations */
	real alphaAndDMaxInIterations[2];					/* maximum values for alpha and d allowed during iterations */
	real costNow;										/* present cost */
	real costNew;										/* new cost */
	real delta[2];										/* increment for alpha and d in one step */
	real lambda = 0.001;								/* lambda */
	real lambdaDown = 10;								/* with a succesful step, lambda is divided by this factor */
	real lambdaUp = 5;									/* with an unsuccesful step, lambda is multiplied by this factor */
	real convCriterion = 0.001;							/* relative convergence criterion */
	int iterMax = 40;									/* maximum number of iterations allowed */
	int i = 1;											/* iteration count */
	
	/* A and B are now defined with these values */
	if (AAndB[0] > AAndB[1] || AAndB[0] < 1.0 || AAndB[1] < 1.0) {
		alphaAndD[0] = 0.0;
		alphaAndD[1] = pow(AAndB[0],1.3483146);
		return;
	}
	
	/* with these values, alpha and d must be at the maximums */
	substVecFromScal(alphaAndDMax,maxAlpha,maxD2/powerLawD1);
	calculateAAndBParameters(AAndBMax,alphaAndDMax);
	
	if (AAndB[1] > AAndBMax[1] || AAndB[0] > AAndBMax[0]) {
		substVecFromVec(alphaAndD,alphaAndDMax);
		return;
	}
		
	/* present alpha, d, A, and B */
	substVecFromVec(alphaAndDNow,alphaAndDGuess);
	calculateAAndBParameters(AAndBNow,alphaAndDNow);
	
	/* limiting alpha and d with these values */
	substVecFromScal(alphaAndDMinInIterations,-10000.0,0.0);
	substVecFromScal(alphaAndDMaxInIterations,10000.0,1000.0);
	
	/* iterations */
	while (isNotConverged(AAndBNow,AAndB,convCriterion) && i < iterMax ) {
		/* computing matrices used in calculations of new guesses for alpha and d */
		calculateJAAndB(J,alphaAndDNow);
		transponateMatrix(JT,J);
		multiplyMatrices(JTJ,JT,J);
		
		calculateAAndBParameters(AAndBNow,alphaAndDNow);
		r[0] = AAndB[0]-AAndBNow[0];
		r[1] = AAndB[1]-AAndBNow[1];
		
		multiplyMatrixAndVec(costGradient,JT,r);
		addLambdaToMatrix(g,JTJ,lambda);
		
		costNow = SQR(r[0]) + SQR(r[1]);
		
		inverseMatrix(invG,g);
		multiplyMatrixAndVec(delta,invG,costGradient);
		
		/* new guesses for alpha and d */
		alphaAndDNew[0] = alphaAndDNow[0] + delta[0];
		alphaAndDNew[1] = alphaAndDNow[1] + delta[1];
		
		/* corresponding A and B */
		calculateAAndBParameters(AAndBNew,alphaAndDNew);
		
		/* corresponding distance */
		rNew[0] = AAndB[0]-AAndBNew[0];
		rNew[1] = AAndB[1]-AAndBNew[1];
		
		/* corresponding cost */
		costNew = SQR(rNew[0]) + SQR(rNew[1]);
		
		/* if lower cost is obtained and alpha and d are within limits, the step is succesful */
		if (costNew < costNow && isAlphaAndDWithinLimits(alphaAndDNew,alphaAndDMinInIterations,alphaAndDMaxInIterations)) {
			substVecFromVec(alphaAndDNow,alphaAndDNew);
			lambda /= lambdaDown;
			
			calculateAAndBParameters(AAndBNow,alphaAndDNow);
			
		/* unsuccesful step */
		} else {
			lambda *= lambdaUp;
			
			/* with very high values, the loop must be breaked due to numerical reasons */
			if (lambda > 1.0e10) {
				break;
			}
		}
			
		/* counting iterations */
		++i;
	}

	/* after succesful iterations, using present values as the result; otherwise, using initial guesses */
	if (isAlphaAndDWithinLimits(alphaAndDNow,alphaAndDMinInIterations,alphaAndDMaxInIterations) && (i!=iterMax)) {
		substVecFromVec(alphaAndD,alphaAndDNow);
		return;
		
	} else {
		substVecFromVec(alphaAndD,alphaAndDGuess);
		return;
	}
	
	return;
}

