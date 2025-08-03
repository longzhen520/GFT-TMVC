%listed as follows:
%1. t-SVD-MSC
%2. GFT-TMVC
%3. t-SVD-MSC-with-GFT

clc;
clear;
close all;

EN_GFT_TMVC      =1;
EN_t_SVD_MSC    =1;

addpath('tools','tproduct toolbox 2.0','ClusteringMeasure','Datasets');

dataset_names = ["Yale","MSRC","EYaleB10_mtv","NGs","WikipediaArticles","ORL"];
  
for ds = 1:1:length(dataset_names)
    
    %% load datasets
        dataName = dataset_names{ds};   
        fprintf('\n Dataset:%s \n',dataName);
        data = dataName;
         load("Datasets\"+data)
        V=length(X);
        cls_num = length(unique(gt));
%data preparation...
        for v=1:V
            [X{v}]=NormalizeData(X{v});
        end
% %  %% resuffer the samples index
   
        N = size(X{1},2); %sample number
%         rng('default');
%         Nind=randperm(N);
        for v=1:V
            X{v}=X{v}(:,Nind);
        end
        gt=gt(Nind);
%% Proform algorithm
methodname={'t-SVD-MSC','GFT-TMVC'};


i=0;
enList=[];
i=i+1;
%% t-SVD-MSC
if EN_t_SVD_MSC
    disp(['performing ', methodname{i},'...']);
    addpath(genpath('t-SVD-MSC'));
    
     tic;
     S{i}=function_tSVD_MSC(X,dataName);
     Time(i)=toc;
    C = SpectralClustering(S{i},cls_num);% C = kmeans(U,numClust,'EmptyAction','drop');
    ACC(i) = Accuracy(C,double(gt));
        disp('...')
    enList = [enList,i];
end
 i=i+1;
%% GFT-TMVC
if EN_GFT_TMVC
    disp(['performing ', methodname{i},'...']);
    addpath(genpath('GFT-TMVC'));
            %% parameter settings
       paras_gft.miu=1.5;
       paras_gft.gt=gt;  
       tic;
      [S{i},Fin]=function_GTNN_MVC(X,paras_gft,dataName);
        Time(i)=toc;       
    Z=0.5*(abs(S{i})+abs(S{i}'));
    C = SpectralClustering(Z,cls_num);
    ACC(i) = Accuracy(C,double(gt));
        disp('...')
    enList = [enList,i];
end

fprintf('================== Result =====================\n');
fprintf(' %6.18s \t   %5.8s \t %5.8s   \n','method','ACC','Time');
for i = 1:length(enList)
    fprintf(' %6.18s \t   %5.3f \t %5.3f    \n',...
    methodname{enList(i)},ACC(enList(i)), Time(enList(i)));
end
end




