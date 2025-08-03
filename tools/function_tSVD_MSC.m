function  [S]=function_tSVD_MSC(X,data)

         Lambda=[5e-4,5e-3,5e-2,0.5,1,1e1,1e2]; 
         switch data
             case "Yale"
                    lambda=Lambda(4);
             case "MSRC"
                    lambda=Lambda(6);
             case "EYaleB10_mtv"
                    lambda=Lambda(7);
             case "synthetic3d"
                    lambda=Lambda(1);
             case "synthetic3d1"
                    lambda=Lambda(1);
             case "NGs"
                    lambda=Lambda(2);
             case "WikipediaArticles"
                    lambda=Lambda(1);
             case "ORL"
                    lambda=Lambda(4);
             case "COIL20MV"
                    lambda=Lambda(6);
         end

                  [S] = tSVD_MSC(X,lambda); % joint affinity matrix

end