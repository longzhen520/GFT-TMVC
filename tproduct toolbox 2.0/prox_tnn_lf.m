function [X,tnn,trank] = prox_tnn_lf(Y,rho,rank_low)
[n1,n2,n3] = size(Y);
if rank_low<round(n3/2)
 r3=rank_low;
 if r3==0
     r3=1;
 end
else
    r3=round(n3/2);
end
max12 = max(n1,n2);
X = zeros(n1,n2,n3);
Y = fft(Y,[],3);
tnn = 0;
trank = 0;
        
% first frontal slice
[U,S,V] = svd(Y(:,:,1),'econ');
S = diag(S);
S = max(S-rho,0);
tol = max12*eps(max(S));
r = sum(S > tol);
S = S(1:r);
X(:,:,1) = U(:,1:r)*diag(S)*V(:,1:r)';
tnn = tnn+sum(S);


% i=2,...,halfn3
halfn3 = r3;
for i = 2 : halfn3
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    S = max(S-rho,0);
    tol = max12*eps(max(S));
    r = sum(S > tol);
    S = S(1:r);
    X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';    
    X(:,:,n3+2-i) = conj(X(:,:,i));
    tnn = tnn+sum(S)*2;
    trank = max(trank,r);
end

% if n3 is even
if mod(n3,2) == 0 && rank_low>round(n3/2)
    i = halfn3+1;
    [U,S,V] = svd(Y(:,:,i),'econ');
    S = diag(S);
    S = max(S-rho,0);
    tol = max12*eps(max(S));
    r = sum(S > tol);
    S = S(1:r);
    X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
    tnn = tnn+sum(S);
    trank = max(trank,r);
end
tnn = tnn/n3;
X = ifft(X,[],3);
