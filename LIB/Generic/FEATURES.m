classdef FEATURES
 
 methods(Static = true)
  
  function X = selGrid(X,G,n)
   for k = 1:numel(X)
    idx = G{k}{2} == n;
    X{k}{2} = X{k}{2}(idx,:);
   end
  end
  
  function Y = subsample(X,S)
   L = cell2mat(X(:,3));
   Y = cell(L(end),1);
   for k = 1:L(end)
    x = cell2mat(X(k == L,2));
    N = size(x,1);
    
    if N > S
      idx = randperm(N,S);
    else
      warning('class with few samples')
      idx = 1:N;
    end
    Y{k} = x(idx,:);
   end
   
   idx = cell2mat(cellfun(@(x) ~ isempty(x),Y,'UniformOutput',false));
   
   Y = cell2mat(Y(idx));
  end
  
  
  function [X] = whiten(X,fudgefactor)
   X = bsxfun(@minus, X, mean(X));
   A = X'*X;
   [V,D] = eig(A);
   X = X*V*diag(1./(diag(D)+fudgefactor).^(1/2))*V';
  end
  
  function epsilon = mnv(X,Y)
   epsilon = min(min(X(:)),min(Y(:)));
   epsilon = abs(epsilon) + 0.5;
  end
  
  function X = mapwhite(X,V,D)
   X = X * V * diag(1./(diag(D)+eps));
  end
  
  function [V,D] = pcaWhite(X,k)
   
   A = bsxfun(@minus, X, mean(X));
   A = A'*A;
   
   [V,D] = eig(A);
   
   V = V(:,1:k);
   D = D(1:k,1:k);
  end
  
 end
end