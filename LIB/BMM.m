classdef BMM
 
 properties(Constant = true)
  alpha_ =  0.25;
 end
 
 methods(Static = true)

  
  %% Visual Map
  
  function [Bag,Lbl] = Bag2FisherVector(Bag,lambda,F)
   
   Lbl = cell(numel(Bag),1);
   for k = 1:numel(Bag)
    A = Bag{k};
    parfor j = 1:numel(Bag{k})
     G = BMM.FisherScore(A{j},lambda);
     A{j} = BMM.FisherVector(F,G);
     %G = BMM.binFisherScore(A{j},lambda);
     %A{j} = G;
    end
    Lbl{k} = ones(numel(Bag{k}),1) * k;
    Bag{k} = cell2mat(A);
   end
   Bag = cell2mat(Bag);
   Lbl = cell2mat(Lbl);
  end
  
  function [lambda,F] = FisherBMM(varargin)
   
   X  = varargin{1};
   N  = 5;
   it = 20;
   
   switch nargin
    case 2
     N = varargin{2};
    case 3
     N  = varargin{2};
     it = varargin{3};
   end
   
   lambda = struct;
   lambda = BMM.EMbin(X,N,lambda,it);   
   F = BMM.FisherImat(lambda,size(X,1));
  end
  
  
  %% Fisher Encoding
  
  
  function G = binFisherScore(X,lambda)
   T = size(X,1);
   N = numel(lambda.w);
   D = size(lambda.m,2);
   G = zeros(1,D*N);
   x = cell(1,N);
   y = zeros(1,N);
   for t = 1:T
    y = y + BMM.multibinomial(X(t,:),lambda.w,lambda.m)';
    for i = 1:N
     x{i} = BMM.dev(X(t,:),lambda.m(i,:));
    end
    G = G + cell2mat(x);
   end
   G = [y/T > 0.015,logical((sign(G)+1)/2)];
  end
  
  
  function G = FisherScore(X,lambda)
   T = size(X,1);
   N = numel(lambda.w);
   D = size(lambda.m,2);
   G = zeros(1,D*N);
   x = cell(1,N);
   for t = 1:T
    y = BMM.multibinomial(X(t,:),lambda.w,lambda.m);
    for i = 1:N
     x{i} =  y(i) * BMM.dev(X(t,:),lambda.m(i,:));
    end
    G = G + cell2mat(x);
   end
   G = G/T;
  end
  
  function z = FisherVector(F,G)
   z = (F.^(-1/2)) .* G;
   z = sign(z) .* (abs(z) .^ BMM.alpha_);
   z = z/norm(z);
  end
  
  
  %% Bernoulli Multivariate Mixture
  
  function F = FisherImat(lambda,T)
   N = numel(lambda.w);
   F = cell(1,N);
   for i = 1:N
    F{i} = T * lambda.w(i) * ( ...
     sum(bsxfun(@times,lambda.w,lambda.m)) ./ (lambda.m(i,:).^2) + ...
     sum(bsxfun(@times,lambda.w,1-lambda.m)) ./ ((1-lambda.m(i,:)).^2));
   end
   F = cell2mat(F);
  end
  
  function lambda = EMbin(Xs,N,lambda,it)
   [T,D]= size(Xs);
   lambda.w = ones(N,1)/N;
   lambda.m = rand(N,D)/2+0.25;   
   gamma = zeros(T,N);
   for k = 1:it
    % Expectation
    w = lambda.w;
    m = lambda.m;
    parfor t = 1:T 
     gamma(t,:) = BMM.multibinomial(Xs(t,:),w,m);
    end
    % Maximization
    S = sum(gamma);
    lambda.w = (S/T)';
    for i = 1:N
     lambda.m(i,:) = 1/S(i) * sum(bsxfun(@times,gamma(:,i),Xs));
    end
   end
  end
  
  %% kernel functions  
  
  function y = multibinomial(x,w,m)
   y = zeros(numel(w),1);
   for k = 1:numel(w)
    y(k) = w(k) * prod(binopdf(x,1,m(k,:)));
   end
   y = y/sum(y);
  end
  
  function y = dev(x,m)
   a = (-1) .^ (1-x);
   b = (m .^ x) .* ( (1-m) .^ (1-x)); 
   y = a./ (b+eps);
  end
  
  function cls = simFisher(x,y,l,N,k)
   
   op = bsxfun(@xor,x,y);
   
%    op(:,N+1:end) = ~ op(:,N+1:end);
   
   op = sum(op,2);
   
   [~,idx] = sort(op);
   
   cls = mode(l(idx(1:k)));
   
  end
  
 end
end