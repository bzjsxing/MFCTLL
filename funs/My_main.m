function [alpha,label] = My_main(X,Y,nv,nc,anchor_rate,p,lambda1,r,IterMax)
N = size(X{1},1);
M = fix(N*anchor_rate);
alpha = repmat(1/nv, [1,nv]);
betaf = ones(nv, 1); 

%% initialize
for v = 1:nv
    H{v} = eye(N,nc);
    G{v} = ones(N,nc);
    J{v} = eye(N,nc);
    Y1{v} = zeros(N,nc);
    SumHB{v} = ones(1,1);
    OBJ{v} = ones(1,1);
end

mu1 = 1e-3;
max_mu = 1e9;

coe = 1.1;
final_result = zeros(1,7);
sX = [N, nc, nv];
epson = 1e-3;
iter = 0;
Isconverg = 0;
time_start = clock;
opt1. style = 1;
opt1. IterMax = IterMax;
opt1. toy = 0;

%% anchor grph
[S] = Construct_Graph(X,nc,anchor_rate, opt1,10);
clear X

while(Isconverg == 0) 
    %% update J
    for v =1:nv
        HY1{v} = H{v} + Y1{v}./mu1;
    end
    HY1_tensor = cat(3,HY1{:,:});
    [myj, ~] = wshrinkObj_weight_lp(HY1_tensor(:), lambda1*betaf./mu1,sX, 0,3,p);
    J_tensor = reshape(myj, sX);
    for k=1:nv
        J{k} = J_tensor(:,:,k);
    end
    clear J_tensor
    
    %% update H{v}
    for v = 1:nv
        [H_IRW, sum1,obj] = IRW(S{v}, J{v}-Y1{v}/mu1, H{v}, 1500, mu1/(2*alpha(v)^r),S{v}*S{v}');
        H{v} = H_IRW;
        SumHB{v} = sum1;
        OBJ{v} = obj;
    end


    %% update alpha
    sum_z = 0;
    sum_a = 0;
    for v = 1:nv       
        z{v} = SumHB{v}.^(1/(1-r));
        sum_z = sum_z+z{v};
    end
    
    for v = 1:nv
        alpha(v) = (z{v})/sum_z;
        sum_a = sum_a+alpha(v)^r;
    end

    %% update Y1{v},Y2{v}
    for v = 1:nv
        Y1{v} = Y1{v}+mu1*(H{v}-J{v});
    end
    
    %% update mu,pho
    mu1 = min(mu1*coe, max_mu);

    %% Clustering result
    
    H_sum = zeros(N,nc);
    HH_sum = zeros(M,nc);
    for v=1:nv
        H_sum = H_sum+H{v}*(alpha(v)^r);
        HH_sum = HH_sum+(S{v}'*H{v})*(alpha(v)^r);
    end

    [~, label] = max(H_sum, [], 2);

    
    %% converge
    for v=1:nv
        history.norm_H_J{v}= norm(H{v}-J{v},inf);
        fprintf('norm_H_J %7.10f \n', history.norm_H_J{v});

        Isconverg = 0;

    end

    for v=1:nv
        if (norm(H{v}-J{v},inf)>epson)  
            Isconverg = 0;
            break;
        else
            Isconverg  = 1;
        end
    end


    %%
    if (iter > IterMax)
        Isconverg  = 1;
    end
    fprintf('iter:%d\n',iter)
    iter = iter + 1;

end

time_end = clock;
fprintf('Time_all:%f s\n',etime(time_end,time_start))
fprintf('Time_average:%f s\n',etime(time_end,time_start)/iter)
fprintf('Final_iter:%d\n',iter)
    
end

