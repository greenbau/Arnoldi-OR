# Arnoldi-OR
Matlab code for Arnoldi-OR
arnoldior_basic.m implements Arnoldi-OR solving a least squares problem from scratch at each iteration.  
arnoldior.m is the same as arnoldior_basic.m except that it reuses Givens rotations to solve the least squares problem.
In exact arithmetic, arnoldior_basic.m and arnoldior.m should give the same answers.  Both codes require the coefficients
  (in monomial basis) of the numerator and denominator polynomials.
arnoldior2_basic.m and arnoldior2.m are the same as the above codes except that the numerator and denominator
   are specified by their roots instead of their coefficients.
