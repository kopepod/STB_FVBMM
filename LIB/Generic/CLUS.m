classdef CLUS
 
 methods(Static = true)
  
  function ctd = binkmeans(X,K,maxit)
   
   [ctd,N] = CLUS.initRandSample(X,K);
   
   opctd = cell(maxit,1);
   opd   = zeros(maxit,1);
   
   for i = 1:maxit
    
    cls = cell(N,1);
    d   = cell(N,1);
    
    parfor j = 1:N
     [cls{j},d{j}] = CLUS.distHamming2(X(j,:),ctd);
    end
    
    cls = cell2mat(cls);
    d   = cell2mat(d);
    opd(i) = sum(d);
    
    fprintf('%d\n',opd(i));
    opctd{i} = ctd;
    
    ctd = CLUS.upctd(X,ctd,cls,K);
    
   end
   
   plot(opd)
   figure
   
   [~,opd] = min(opd);
   
   ctd = opctd{opd};
   
  end
  
  function ctd = upctd(X,ctd,cls,K)
   
   for k = 1:K
    idx = cls == k;
    if sum(idx) < 2
     continue
    end    
    y = X(idx,:);
    ctd(k,:) = mean(y) > 0.5;
   end
      
  end
  
  function [cls,d] = distHamming2(x,ctd)
   
   M = numel(x)/2;
   
   a = bsxfun(@xor,x(1:M),ctd(:,1:M));
   
   b = bsxfun(@and,x(M+1:end),ctd(:,M+1:end));
   
   c = and(a,b);
   
   d = sum(c,2);
   
   [d,cls] = min(d);
   
  end
  
  
  function [ctd,N] = initRandSample(X,K)
   N = size(X,1);
   idx = ceil(K * rand(N,1));
   ctd = cell(K,1);
   for k = K:-1:1
    y = X(idx == k,:);
    y = mean(y) > 0.5;
    if numel(y) < 2
     ctd(k) = [];
     continue
    end
    ctd{k} = y;
   end
   ctd = cell2mat(ctd);
  end
  
 end
end
  