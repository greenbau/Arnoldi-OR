function [RofAb,Q,H,ls_resids,true_resids] = arnoldior2(A,N,D,b,kmax)
%
%  Given an mxm matrix A, an m-vector b, a numerator polynomial N, and
%  a denominator polynomial D, this routine finds the optimal approximation
%  (in D(A)'D(A)-norm) to inv(D(A))N(A)b from the kmax dimensional Krylov space 
%  span{b, Ab, ..., A^(kmax-1) b}.  It returns the result in RofAb, 
%  the residual norms from from the least squares problem at each step in array ls_resids,
%  and the true residual norms ||N(A)b - D(A) x_k|| in array true_resids.  

%  The roots of the polynomials N(z) and D(z) are stored in arrays N and D,
%  along with the coefficient of the highest power of z:
%  N(A) = N(1)*(A - N(2) I)*...*(A - N(end) I),
%  D(A) = D(1)*(A - D(2) I)*...*(A-D(end) I).
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

cosines = ones(kmax,degD); sines = zeros(kmax,degD);    % Store rotations at each step.

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
    for j=1:degN,
      NMat = (H(1:k,1:k) - N(j+1)*eye(k))*NMat;
    end;

%      Store relevant pieces of DMat and NMat in calD and eta, which will be updated
%      at each step.

    if k==1+deg,
      calD = DMat; eta = NMat(:,1);
    else
      calD = [calD, DMat(1:k-1,k-deg)]; calD = [calD; DMat(k,:)]; eta = [eta; 0];
    end;

%      If k > 1+deg, apply previous rotations to last column of calD.

    for kk=1:k-deg-1,          % Rotations from steps deg+1 to k-1.
      for ii=kk+degD:-1:kk+1,  % Applied to rows kk+degD,...,kk+1
        tempkk = cosines(kk,ii-kk)*calD(kk,k-deg) + sines(kk,ii-kk)*calD(ii,k-deg);
        tempii = -conj(sines(kk,ii-kk))*calD(kk,k-deg) + cosines(kk,ii-kk)*calD(ii,k-deg);
        calD(kk,k-deg) = tempkk; calD(ii,k-deg) = tempii;
      end;
    end;

%      Now determine new rotations and apply to last col of calD and to eta.

    for ii=k-deg+degD:-1:k-deg+1,
%        There are more stable formulas for these rotations.
      if abs(calD(k-deg,k-deg)) ~= 0,
        denom = sqrt(abs(calD(k-deg,k-deg))^2 + abs(calD(ii,k-deg))^2);
        cosines(k-deg,ii-(k-deg)) = abs(calD(k-deg,k-deg))/denom;
        signf = calD(k-deg,k-deg)/abs(calD(k-deg,k-deg));
        sines(k-deg,ii-(k-deg)) = signf*conj(calD(ii,k-deg))/denom;
      end;
%      [G,Y] = planerot([calD(k-deg,k-deg); calD(ii,k-deg)]);
%      cosines(k-deg,ii-(k-deg)) = G(1,1); sines(k-deg,ii-(k-deg)) = G(1,2);
      tempk = cosines(k-deg,ii-(k-deg))*calD(k-deg,k-deg) + ...
              sines(k-deg,ii-(k-deg))*calD(ii,k-deg);
      tempii = -conj(sines(k-deg,ii-(k-deg)))*calD(k-deg,k-deg) + ...
              cosines(k-deg,ii-(k-deg))*calD(ii,k-deg);
      calD(k-deg,k-deg) = tempk; calD(ii,k-deg) = tempii;
      tempk = cosines(k-deg,ii-(k-deg))*eta(k-deg) + sines(k-deg,ii-(k-deg))*eta(ii);
      tempii = -conj(sines(k-deg,ii-(k-deg)))*eta(k-deg) + cosines(k-deg,ii-(k-deg))*eta(ii);
      eta(k-deg) = tempk; eta(ii) = tempii;
    end;

%      Solve the upper triangular linear system and form approx solution.

    y = calD(1:k-deg,1:k-deg)\eta(1:k-deg);
    RofAb = Q(:,1:k-deg)*y*normb;

%      Compute residual norm in least squares problem and in RofAb.
    k, resnorm = norm(calD*y - eta)*normb, resnorm_true = norm(NofA*b - DofA*RofAb);
    if abs(resnorm - resnorm_true) > 1.e-8, pause, end;
    ls_resids = [ls_resids; resnorm];
    true_resids = [true_resids; resnorm_true];
  end;
end;
