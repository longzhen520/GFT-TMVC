     dataName = ["yale"];

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
        paras.miu=1.5;
        paras.gt=gt;  
        paras.lambda=0.01;
        paras.gamma=1; 
        [S,Fin]=GTNN_MVC_v3(X,paras);

    