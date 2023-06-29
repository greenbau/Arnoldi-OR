function [RofAb,Q,H,ls_resids,true_resids] = arnoldior_basic(A,N,D,b,kmax)
%
%  Given an mxm matrix A, an m-vector b, a numerator polynomial N, and
%  a denominator polynomial D, this routine finds the optimal approximation
%  (in D(A)'D(A)-norm) to inv(D(A))N(A)b from the kmax dimensional Krylov space 
%  span{b, Ab, ..., A^(kmax-1) b}.  It returns the result in RofAb,
%  the residual norms from from the least squares problem at each step in array ls_resids,
%  and the true residual norms ||N(A)b - D(A) x_k|| in array true_resids.
%
%  The coefficients of the polynomials N and D are stored from lowest to 
%  highest power:  N(A) = N(1) I + N(2) A + ... + N(end) A^(length(N)-1),
%  D(A) = D(1) I + D(2) A + ... + D(end) A^(length(D)-1).  See arnoldior2_basic.m
%  for a version where the roots of N and D are required instead.
%
%  kmax + max(deg(N),deg(D)) must be less than or equal to m.
%
[m,m] = size(A);

degN = length(N)-1; degD = length(D)-1;
deg = max([degN;degD]);

%  First form N(A) and D(A) so that we can compute true residuals.
DofA = D(1)*eye(m); NofA = N(1)*eye(m);
for j=1:degD, DofA = DofA + D(j+1)*A^j; end;
for l=1:degN, NofA = NofA + N(l+1)*A^l; end;

%  Run kmax+deg steps of Arnoldi.

Q = zeros(m,kmax+deg);  H = zeros(kmax+deg,kmax+deg-1);
ls_resids = []; true_resids = [];
normb = norm(b);
Q(:,1) = b/normb;

for k=1:kmax+deg,
  Q(:,k+1) = A*Q(:,k);
  for l=1:k,
    H(l,k) = Q(:,l)'*Q(:,k+1);
    Q(:,k+1) = Q(:,k+1) - H(l,k)*Q(:,l);
  end;
%    Orthogonalize twice.
  for l=1:k,
    Q(:,k+1) = Q(:,k+1) - (Q(:,l)'*Q(:,k+1))*Q(:,l);
  end;
  H(k+1,k) = norm(Q(:,k+1));
  Q(:,k+1) = Q(:,k+1)/H(k+1,k);

%    If k >= 1+deg, set up and solve the least squares problem.
  if k >= 1+deg,
    DMat = zeros(k,k-deg);
    NMat = zeros(k,k-deg);
    for j=0:deg,
      if j==0,  
        Hkpjk = eye(k-deg); 
        DMat(1:k-deg,1:k-deg) = D(1)*Hkpjk; 
        NMat(1:k-deg,1:k-deg) = N(1)*Hkpjk;
      else
        Hkpjk = H(1:k-deg+j,1:k-deg+j-1)*Hkpjm1k;
        if j <= degD,
          DMat(1:k-deg+j,1:k-deg) = DMat(1:k-deg+j,1:k-deg) + D(j+1)*Hkpjk;
        end;
        if j <= degN,
          NMat(1:k-deg+j,1:k-deg) = NMat(1:k-deg+j,1:k-deg) + N(j+1)*Hkpjk;
        end;
      end;
      Hkpjm1k = Hkpjk;
    end;

%      Form approximate solution.
    if k==1+deg,
      calD = DMat; eta = NMat(:,1);
    else
      calD = [calD, DMat(1:k-1,k-deg)]; calD = [calD; DMat(k,:)]; eta = [eta; 0];
    end;
    y = calD\eta;
    RofAb = Q(:,1:k-deg)*y*normb;

%      Compute residual norm in least squares problem and in RofAb.
    k, resnorm = norm(calD*y - eta)*normb, resnorm_true = norm(NofA*b - DofA*RofAb)
    if abs(resnorm - resnorm_true) > 1.e-8, pause, end;
    ls_resids = [ls_resids; resnorm];
    true_resids = [true_resids; resnorm_true];
  end;
end;
