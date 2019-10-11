classdef BMMu64
 
 properties(Constant = true)
  alpha_ =  0.5;
  max_it = 500;
 end
 
 methods(Static = true)

   
  
  %% Visual Map
  
  function Y = encodeBMMfisher(X,lambda,F)
   
   if ischar(X)
    load(X)
   end
   
   Y = cell(size(X));
   Y(:,[1,3]) = X(:,[1,3]);
   
   fprintf('Encoding %d Vectors\n',size(X,1));
   
   w = lambda.w;
   m = lambda.m;
   m1 = 1./(lambda.m+realmin);
   m2 = -1./(1-lambda.m-realmin);
   tic
   parfor k = 1:size(X,1)
     %for k = 1:size(X,1)
%     try
       G = mexfisherscore(X{k,2},w,m,m1,m2);
 %    catch
 %      G = 0;
 %    end
     Y{k,2} = BMMu64.FisherVector(F,G);
   end
   toc
   
  end
  
  %% Estimate Distribution Parameters 
 
  function [lambda,F] = FisherBMM(X,idx,S,N)
   
   if ischar(X)
    load(X)
   end
   
   if ~ isempty(idx)
     X = X(idx,:);
   end
   
   X = FEATURES.subsample(X,S);
   
   fprintf('Fitting %d Features into %d-BMM\n',size(X,1),N);
   
   lambda = struct;
   
   it = BMMu64.max_it;
   tic
   lambda = BMMu64.EMbin(X,N,lambda,it);   
   F = BMMu64.FisherImat(lambda,size(X,1));
   toc
  end
  
  function lambda = EMbin(X64,N,lambda,it)
   [T,D]= size(X64);
   lambda.w = ones(N,1)/N;
   lambda.m = rand(64*D,N)/2+0.25;   
   %m = 0;
   %load('m.mat')
   %lambda.m = m';
   gamma = zeros(N,T);
   
   for k = 1:it
    % Expectation
    w = lambda.w;
    m = lambda.m;
    parfor t = 1:T 
    %for t = 1:T 
     gamma(:,t) = mexmultibinomial(X64(t,:),w,m);
    end
    % Maximization
    S = sum(gamma,2);
    lambda.w = (S/T);
    
    c = cell(N,1);
    
    for i = 1:N
     c{i} = lambda.m(:,i);
    end
    
    parfor i = 1:N
    %for i = 1:N
     c{i} = 1/S(i) * mexmeanbitwise64(X64,gamma(i,:));
    end
    
    for i = 1:N
     lambda.m(:,i) = c{i};
    end
    
    es = sum(abs(lambda.m(:) - m(:)));
    
    if es < 9
     fprintf('BMM early stop: exiting\n');
     return
    else
     if mod(k,20) == 0
      fprintf('iteration %d error: %f\n',k,es);
     end
    end
    
    %{
    if e0 > es
     e0 = es;
    else
     fprintf('BMM error minimum reached: exiting\n');
     return
    end
    %}
   end
   fprintf('BMM max number of iterations reached\n');
  end
  
  function F = FisherImat(lambda,T)
   N = numel(lambda.w);
   lambda.m = lambda.m';
   F = cell(1,N);
   for i = 1:N
    F{i} = T * lambda.w(i) * ( ...
     sum(bsxfun(@times,lambda.w,lambda.m)) ./ (lambda.m(i,:).^2+eps) + ...
     sum(bsxfun(@times,lambda.w,1-lambda.m)) ./ ((1-lambda.m(i,:)).^2+eps));
   end
   F = cell2mat(F);
  end
  
  
  %% Fisher Encoding
    
  function G = FisherScore(X64,lambda)
   [T,D] = size(X64);
   N = numel(lambda.w);
   G = zeros(64*D*N,1);
   x = cell(N,1);
   m = lambda.m;
   w = lambda.w;
   for t = 1:T
    z = X64(t,:);
    y = mexmultibinomial(z,w,m);
    for i = 1:N
     x{i} =  y(i) * mexmultinomaildev(z,m(:,i));
    end
    G = G + cell2mat(x);
   end
   G = G'/T;
  end
  
  function z = FisherVector(F,G)
   z = (F.^(-1/2)) .* G;
   z = sign(z) .* (abs(z) .^ BMMu64.alpha_);
   z = z/norm(z);
  end
  
  
  %% CELL processing
  
   
  function [His,Lbl] = mexBag2FV(dbname,N,sbj,prfx,inc_excl,lambda,F)
  
  sbj = cellfun(@(var) strcat(prfx,var),sbj,'UniformOutput',false);
  
  His = cell(N,1);
  Lbl = cell(N,1);
  c = cell(1);
  
  for k = 1:N
   load(strcat(dbname,num2str(k),'.mat'))
   fprintf('%d ',k);
   
   idx = BAG.findsubjects(c,sbj,inc_excl,prfx);
   c(idx) = [];
   for j = 1:numel(c)
    c{j} = c{j}{2};
   end
   
   %parfor j = 1:numel(c)
   for j = 1:numel(c)
    G = BMMu64.FisherScore(c{j},lambda);
    c{j} = BMMu64.FisherVector(F,G);
   end
   
   Lbl{k} = ones(numel(c),1) * k;
   
   His{k} = cell2mat(c);
   
  end
  fprintf('\n');
  Lbl = cell2mat(Lbl);
  
  His = cell2mat(His);
  
 end
  
  
  function [Bag,Lbl] = Bag2FisherVector(Bag,lambda,F)
   
   Lbl = cell(numel(Bag),1);
   for k = 1:numel(Bag)
    A = Bag{k};
    for j = 1:numel(Bag{k})
     G = BMMu64.FisherScore(A{j},lambda);
     A{j} = BMMu64.FisherVector(F,G);
     %G = BMMu64.binFisherScore(A{j},lambda);
     %A{j} = G;
    end
    Lbl{k} = ones(numel(Bag{k}),1) * k;
    Bag{k} = cell2mat(A);
   end
   Bag = cell2mat(Bag);
   Lbl = cell2mat(Lbl);
  end
  
  
  
  
 end
end
