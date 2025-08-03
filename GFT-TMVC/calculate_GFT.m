function kerN=calculate_GFT(Z,n)
  
        V=length(Z);
        S=0;
        for v=1:V
        S =S+ abs(Z{v})+abs(Z{v}');
        end
        S=S/(2*V);

        N = size(S,1);
         DN = diag( 1./sqrt(sum(S,2)) );
        LapN = speye(N)- DN * S * DN;
        %     LapN = DN * S * DN;
%         [uN,sN,vN] = svd(LapN);
%         kerN = uN(:,N-n:N);
         [kerN, ~] =  eig(LapN);
end