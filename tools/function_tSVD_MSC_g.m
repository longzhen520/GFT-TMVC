function  [S]=function_tSVD_MSC_g(X,data,Fin)

         Lambda=[5e-4,5e-3,5e-2,0.5,1,1e1,1e2]; 

         switch data
             case "Yale"
                    lambda=Lambda(5);
             case "MSRC"
                    lambda=Lambda(5);
             case "EYaleB10_mtv"
                    lambda=Lambda(7);
             case "synthetic3d"
                    lambda=Lambda(3);
             case "NGs"
                    lambda=Lambda(5);
             case "WikipediaArticles"
                    lambda=Lambda(1);
             case "ORL"
                    lambda=Lambda(6);
         end
      
                  [S] = tSVD_MSC_g(X,lambda,Fin); % joint affinity matrix
       
end