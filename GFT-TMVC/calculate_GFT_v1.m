function kerN=calculate_GFT_v1(Z,R)
  
        V=length(Z);
       
        Z =0.5*(abs(Z)+abs(Z'));
    

        N = size(Z,1);
        DN = diag( 1./sqrt(sum(Z)) );
        LapN = max(speye(N)- DN * Z * DN,0);
        [kerN, ~] =  eig(LapN); 
       % , R, 'smallestreal');
        %     LapN = DN * S * DN;
%         [uN,sN,vN] = svd(LapN);
%         kerN = vN(:,1:n);
end