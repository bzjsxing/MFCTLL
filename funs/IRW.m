% B      n*m
% H      n*c
% U0      n*c initializaton 

% maxiter 
% lambda parameter
% BB=B'*B
function  [U,sum1,want_value] = IRW(B, H, U0, maxiter, lambda, BB)
% function [U,sum1,want_value] = PGD(B, H, U0, alpha, maxiter, lambda, BB)
    O = B';
    clear B
    B = O;
    clear O
    U=U0;
    %%

    %[m,c]=size(U);
    %%
    [n,c]=size(U);
    %%

    F=B*U;
    E=F'*F;%f'f

    %%
%     E=U'*U;

    d=diag(E).^(1/2);
    obj(1)=sum(d)-lambda*norm(U-H,'fro')^2;
for iter=1:maxiter    
    
%     A=BB*U*diag(1./(d+eps))+2*lambda*H;
    A=BB*U*diag(1./(d+eps))+2*lambda*H;
    AA=1/(2*lambda)*A;
   
   
     for i=1:n
     %   U(i,:)= EProjSimplex_new(AA(i,:));
      U(i,:)=  EProjSimplex(AA(i,:));
%         f(i)=trace(A'*U)-lambda*trace(A'*U);
     end

  %%

    F=B*U;
    E=F'*F;
    %%
    d=diag(E).^(1/2);
    sum1 = sum(d); 
    obj(iter+1)=sum(d)-lambda*norm(U-H,'fro')^2;
    want_value = sum(d)-lambda*norm(U-H,'fro')^2;

%     fprintf('obj(iter+1)-obj(iter) %7.10f \n', obj(iter+1)-obj(iter));

    if (iter>1 && abs(obj(iter+1)-obj(iter)) < 10^-6)
       break;
    end
    
end
%         O = B';
%         clear B
%         B = O;
%         clear O
%  plot(obj);

%  [~,Label]=max(F,[],2);
 


