    clear all; 
    close all; 
    clc;

    addpath('Tensor-tensor-product-toolbox-master-2', 'mylib', 'datasets');
    %% load data "yale",  "yale","MSRC","EYaleB10_mtv",,"3sources","ORL","COIL20MV","3sources",
    dataset_names = ["synthetic3d","yale","MSRC","EYaleB10_mtv","NGs","ORL","WikipediaArticles",];

    for ds = 1:1:length(dataset_names)
    %% load datasets
        dataName = dataset_names{ds};   
        fprintf('\n Dataset:%s \n',dataName);
        data = dataName;
        load("E:/multi-view-complete/New Version1/Pdatasets/"+data)
    %       gt=Y;
        V=length(X);
        cls_num = length(unique(gt));
%data preparation...
        for v=1:V
            [X{v}]=NormalizeData(X{v});
        end
% %  %% parameter setting
        N = size(X{1},2); %sample number
        for v=1:V
            X{v}=X{v}(:,Nind);
        end
        gt=gt(Nind);
        
        %% parameter settings

        Lambda=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,1e2,1e3]; 
        Gamma=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,1e2,1e3];  
%         Beta=[1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1,10,1e2,1e3];  

        resultsAll=[];
        paras.gt=gt;
        paras.miu=1.5;
         for g=1:length(Gamma)
                for l=1:length(Lambda)
                    paras.lambda=Lambda(l);
                    paras.gamma=Gamma(g);
%                     paras.beta=Beta(b);
                    tic;
                    [S,max_ACC]=GTNN_MVC_v3(X,paras);
                        rng('default')
                        Z=0.5*(abs(S)+abs(S'));
                        C1 = SpectralClustering(Z,cls_num);   
                        ACC1(g,l)= Accuracy(C1,double(gt));
%                         max_ACC1(g,l)=max_ACC;
                        result = [Gamma(g),Lambda(l),ACC1(g,l)];
                        resultsAll = [resultsAll; result];    
%           end
                end
         end
          
save("clustering_results_without_GFT"+"_"+dataName+"_gtnn.mat",'resultsAll');        
 end



