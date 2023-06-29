function [RofAb,Q,H,ls_resids,true_resids] = arnoldior2_basic(A,N,D,b,kmax)
%
%  Given an mxm matrix A, an m-vector b, a numerator polynomial N, and
%  a denominator polynomial D, this routine finds the optimal approximation
%  (in D(A)'D(A)-norm) to inv(D(A))N(A)b from the kmax dimensional Krylov space 
%  span{b, Ab, ..., A^(kmax-1) b}.  It returns the result in RofAb,
%  the residual norms from from the least squares problem at each step in array ls_resids,
%  and the true residual norms ||N(A)b - D(A) x_k|| in array true_resids.
%
%  The roots of the polynomials N(z) and D(z) are stored in arrays N and D,
%  along with the coefficients of the highest power of z:
%  N(A) = N(1)*(A - N(2) I)*...*(A - N(end) I),
%  D(A) = D(1)*(A - D(2) I)*...*(A - D(end) I).
%
%  kmax + max(deg(N),deg(D)) must be less than or equal to m.
%
[m,m] = size(A);

degN = length(N)-1; degD = length(D)-1;
deg = max([degN;degD]);

%  First form N(A) and D(A) so that we can compute true residuals.
DofA = D(1)*eye(m); NofA = N(1)*eye(m);
for j=1:degD, DofA = DofA*(A - D(j+1)*eye(m)); end;
for l=1:degN, NofA = NofA*(A - N(l+1)*eye(m)); end;

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
    DMat = zeros(k,k-deg);   %% We need only the last row and last column of DMat.
    NMat = zeros(k,k-deg);   %% We need only the first column of NMat.
    DMat(1:k-deg,1:k-deg) = D(1)*eye(k-deg);
    NMat(1:k-deg,1:k-deg) = N(1)*eye(k-deg);
    for j=1:degD,
      DMat = (H(1:k,1:k) - D(j+1)*eye(k))*DMat;
    end;
    for l=1:degN,
      NMat = (H(1:k,1:k) - N(l+1)*eye(k))*NMat;
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
