clear all;
addpath([pwd, '/funs']);
addpath([pwd, '/datasets']);
% addpath([pwd, '/BKHK+GRAPH']);

%% load data
dataname='HandWritten4';
load(strcat(dataname,'.mat'));

nv = length(X);
nc = length(unique(Y));

%% Data pre-processing A
disp('------Data preprocessing------');
tic
for v = 1:nv
    a = max(X{v}(:));
    X{v} = double(X{v}./a);
end
toc

%% setting
% The following parameters and settings are used to reproduce the results in the paper 
% 'Multi-view Fuzzy Clustering based on Tensorized Label Learning'.

% MSRC     
% anchor_rate = 0.2;    p = 0.3;    lambda1= 10;    r = 0.9;
% HandWritten4    
anchor_rate = 0.1;    p = 0.4;    lambda1= 200;    r = 0.9;
% Mnist4   
% anchor_rate = 0.3;    p = 0.9;    lambda1= 49;    r = 0.5;
% scene15Big  
% anchor_rate = 0.8;    p = 0.7;    lambda1= 20;    r = 0.5;
% NUS
% anchor_rate = 0.002;    p = 0.3;    lambda1= 1000;    r = 0.5;
% reuters   
% anchor_rate = 0.001;    p = 0.2;    lambda1= 19;    r = 0.1;



%% main
IterMax = 160;
filename=['result-IRW-' dataname '.txt'];
fid = fopen(filename,'a');
for num1 = 1:length(anchor_rate)
    for num2 = 1:length(p)
        for num3 = 1:length(lambda1)
           [alpha,label] = My_main(X,Y,nv,nc,anchor_rate(num1),p(num2),lambda1(num3),r,IterMax);
           final_result = ClusteringMeasure1(Y,label);
           for n_result = 1:length(final_result)
                fprintf(fid, '%f ' ,final_result(n_result));
                fprintf('%f ' ,final_result(n_result));
           end
           fprintf('\n');
           fprintf('anchor_rate=%f_p=%f_lambda1=%f\n', anchor_rate(num1),p(num2),lambda1(num3));
           fprintf(fid, 'anchor_rate=%f_p=%f_lambda1=%f\n', anchor_rate(num1),p(num2),lambda1(num3));

        end 
    end
end

