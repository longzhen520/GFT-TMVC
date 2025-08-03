function  [S,Fin]=function_GTNN_MVC(X,paras,data)

%     Lambda=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,1e2,1e3]; 
%     Gamma=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,1e2,1e3];  
         switch data
             case "Yale"
                   paras.lambda=0.01;
                   paras.gamma=1; 
             case "MSRC"
                    paras.lambda=0.1;
                    paras.gamma=0.01; 
             case "EYaleB10_mtv"
                    paras.lambda=1000;
                    paras.gamma=0.01; 
             case "synthetic3d"
                    paras.lambda=0.01; 
                    paras.gamma=0.01; 
             case "synthetic3d1"
                    paras.lambda=0.01; 
                    paras.gamma=0.01; 
             case "NGs"
                    paras.lambda=1e-4;
                    paras.gamma=1e-4;
             case "WikipediaArticles"
                    paras.lambda=1e-4;
                    paras.gamma=100; 
             case "ORL"
                    paras.lambda=1e-3;
                    paras.gamma=1e-3; 
             case "COIL20MV"
                    paras.lambda=1e-6;
                    paras.gamma=1e-5; 
         end
                     [S,Fin]=GTNN_MVC_v3(X,paras);
end