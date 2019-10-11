classdef BVLAD
 properties(Constant = true)
  replicates = 16;
 end
 
 methods(Static = true)
  
  function [ctd,qk] = BoF(X,idx,S,K)
   
   if ischar(X)
    load(X)
   end
   
   X = X(idx,:);
   Xs = FEATURES.subsample(X,S);
   
   fprintf('Clustering %d Features into %d-cluster with %d replicates\n',size(Xs,1),K,BVLAD.replicates);
   
   ctd = cell(BVLAD.replicates,1);% centroids
   ccm = zeros(BVLAD.replicates,1);% cluster compactness
   %lbl = cell(BVLAD.replicates,1);% cluster compactness
   
   seed = 1e6*rand(BVLAD.replicates,1);
   
   parfor j = 1:BVLAD.replicates
    [~,bov,sumd] = mexbinkmeans(Xs,K,seed(j));
    ctd{j} = bov;
    ccm(j) = sum(sumd);
    %lbl{j} = l;
   end
   
   fprintf('Selecting best cluster (high compactness)\n');
   
   [~,ccm] = min(ccm);
   
   ctd = ctd{ccm};
   ctd(ctd(:,end) == 0,:) = [];
   
   K = size(ctd,1);
   qk = zeros(K,1);
   
   fprintf('Recomputing cluster ...\n');
   xm = cell(size(X,1),1);% x matching to centroid
   
   parfor j = 1:numel(xm)
    %fprintf('Matching vectors [%d/%d] ...\n',j,numel(xm));
    a = mexbinmatch(X{j,2},ctd);
    xm{j} = uint16(a(:,1));
   end
   ctd = cell(K,1);
   for k = 1:K
    xk = cell(size(xm));
    for j = 1:numel(xk)
     xk{j} = X{j,2}(k == xm{j},:);
    end
    xk = cell2mat(xk);
    qk(k) = size(xk,1);
    ctd{k} = mexbinmean(xk);
   end
   ctd = cell2mat(ctd);
  end
  
  
  function X = FeatureMap(X,ctd,qk)
   
   if ischar(X)
    load(X)
   end
   
   fprintf('Mapping %d Vectors ...\n',size(X,1));
   
   qk = qk/sum(qk);
   
   parfor j = 1:size(X,1)
    %X{j,2} = mexBVLADu64(X{j,2},ctd,qk);
    h = mexbinhist(X{j,2},ctd);
    X{j,2} = h'/sum(h);
   end
   
  end
  
  function X = mapvect(X,ctd)
   
   if ischar(X)
    load(X)
   end
   
   for j = 1:size(X,1)
    X{j,2} = BVLAD.BVLADvector(X{j,2},ctd);
   end
   
  end
  
  
  function y = BVLADvector(x,ctd)
   
   x = BVLAD.vect_u64_2_log(x);
   ctd = BVLAD.vect_u64_2_log(ctd);
   
   y = cell(size(x,1),1);
   
   for k = 1:size(x,1)
    u = xor(x(k,:),ctd);
    y{k} = u(:)';
   end
   
   y = cell2mat(y);
   
   y = sum(y);
   
   y = y/sum(y);
   
  end
  
  function y = vect_u64_2_log(x)
   map = 2.^uint64(63:-1:0);
   X = cell(size(x,1),1);
   for k = 1:size(x,1)
    X{k} = BVLAD.u64_2_log(x(k),map);
   end
   y = cell2mat(X);
  end
  
  function y = u64_2_log(x,map)
   y = false(1,64);
   for k = 1:64
    if x >= map(k)
     x = x - map(k);
     y(k) = true;
    end
   end
  end
  
 end
end