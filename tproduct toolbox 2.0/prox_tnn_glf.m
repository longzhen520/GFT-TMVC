function [X,tnn,trank] = prox_tnn_glf(Y,rho,rank_low,t)

%% constract L
      % choice 1
% %       options.Metric = 'Euclidean';
% %       options.NeighborMode = 'KNN';
% %       options.k = 5;
% % %       options.WeightMode = 'Binary'
% %         options.WeightMode = 'HeatKernel';

    options.KernelType='Gaussian';
%               'Gaussian'      - e^{-(|x-y|^2)/2t^2}
%               'Polynomial'    - (x'*y)^d
%               'PolyPlus'      - (x'*y+1)^d
%               'Linear'        -  x'*y
      options.t = t;
      Dim=size(Y); 
      Y_n=reshape(Y,[Dim(1)*Dim(2),Dim(3)]);
      W = constructKernel(Y_n',Y_n',options);
      L=diag(sum(W,2))-W;
     [Fin,S]=eigs(L, rank_low, 'sa');
      F=Fin';

      %% GFT
      Y = tmprod(Y,F,3);
      if rank_low<size(Y,3)
         r3=rank_low;
     else
         r3=size(Y,3);
      end
     X=zeros(size(Y));
     for i=1:r3
        [U,S,V] = svd(Y(:,:,i),'econ');
        S = diag(S);
        S = max(S-rho,0);
        tol =eps(max(S));
        r = sum(S > tol);
        S = S(1:r);
        X(:,:,i) = U(:,1:r)*diag(S)*V(:,1:r)';
     end
   X = tmprod(X,Fin,3);  
