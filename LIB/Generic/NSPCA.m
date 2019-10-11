classdef NSPCA
 
 properties(Constant = true)
  L1_eps = 1e-4;    % iteration stopping criterion
  L1_optimize = 1;  % optimize weights for given sparsity pattern
  L1_NN_eps = 1e-4;    % iteration stopping criterion
  L1_NN_optimize = 1;  % optimize weights for given sparsity pattern
 end
 
 methods(Static = true)
  
  function T = nspca(X,K,op)
   
   X = NSPCA.zeromean(X);
   
   [~, w] = NSPCA.PCA(X,1);
   [~,d] = size(X);
   w = w * w';
   w  = - w;
   
   for j = 1:size(w,1)
    w(j,j) = w(j,j)+1;
   end
   
   %X = (w * X')';
   
   X = NSPCA.largeMatProd(w,X');
   
   X = X';
   
   T = zeros(d,K);
   
   if strcmp(op,'sparse')
    for k = 1:K
     T(:,k) = NSPCA.emPCA_L1(X,k);
    end
   else
    for k = 1:K
     T(:,k) = NSPCA.emPCA_L1_NN(X,k);
    end    
   end
  end
  
  function A = largeMatProd(A,B)
   
   [N,P] = size(A);
   [~,M] = size(B);
   
   if P > M
    R = zeros(1,M);
   else
    R = zeros(1,P);
    A = [A,zeros(N,M-P)];
   end
   
   for j = 1:N
    for i = 1:M
     R(i) = sum(A(j,1:P) * B(:,i));
    end
    
    if P > M
     A(j,1:M) = R;
    else
     A(j,:) = R;
    end
    
    
   end
   
   if P > M
    A = A(:,1:M);
   end
   
  end
  
  function Y = zeromean(X, dim)
   if nargin < 2 || dim == 1
    Y = X - repmat(mean(X),size(X,1),1);
   else
    Y = X - repmat(mean(X,2),1,size(X,2));
   end
  end
  
  function w = emPCA_L1(X, n_nz)
   % Expectation maximization for sparse PCA.
   %
   % INPUTS
   % X: (n x d) zero mean data matrix (samples are rows)
   % n_nz: number of nonzero coefficients in w
   %
   % OUTPUTS
   % w: sparse loadings of first principal component (normalized to unit length)
   %
   % Tuneable parameters are set in the global variable PRM.
   %
   % Copyright 2013 Christian Sigg (christian@sigg-iten.ch)
   % See license file for details.
   %
   % This code is no longer actively maintained. See my 'nsprcomp' R package
   % for a better documented and more feature complete implementation.

   epsilon_ = NSPCA.L1_eps;            % iteration stopping criterion
   optimize_ = NSPCA.L1_optimize;  % optimize weights for given sparsity pattern
   
   [~,d] = size(X);
   
   [~, w] = NSPCA.PCA(X,1);  % initialize with (full) loadings of first PC
   w_old = zeros(size(w));
   while(abs(w'*w_old) < 1 - epsilon_)
    w_old = w;
    
    % E step: orthogonally project on w
    y = X*w;
    
    % M step:
    w = NSPCA.isoQP_L1(y'*y, -(sum(X.*repmat(y,1,d)))', n_nz);
    w = w / norm(w);
   end
   
   if optimize_
    indx = abs(w) > 0;
    [~, v] = NSPCA.PCA(X(:,indx),1);
    w = zeros(d,1);
    w(indx) = v;
   end
  end

  function w = isoQP_L1(h, f, n_nz)
   % Solve min_w h/2*w'*w + f'*w, subject to norm(w,1) <= k, but instead of
   % specifying k, the number of nonzero elements n_nz is the criterion.
   %
   % The QP is of simple form due to the Hessian being a scaled identity
   % matrix, and therefore reduces to minimizing quadratic distance to the
   % unconstrained minimum.
   %
   % INPUTS
   % h: (non-negative) scalar of the Hessian
   % f: linear term
   % n_nz: number of nonzero coefficients in w
   %
   % OUTPUTS
   % w: optimum vector satisfying constraints
   %
   % Copyright 2013 Christian Sigg (christian@sigg-iten.ch)
   % See license file for details.
   %
   % This code is no longer actively maintained. See my 'nsprcomp' R package
   % for a better documented and more feature complete implementation.
   
   d = numel(f);
   w_star = -f/h;
   if (n_nz > d)
    n_nz = d;
   end
   
   % do NN monotone stagewise optimization
   [ms, mi] = sort(abs(w_star), 'descend');
   ms(d+1) = 0;
   
   w = zeros(d,1);
   w(1:n_nz) = ms(1:n_nz) - ms(n_nz+1);
   
   % backtransform solution w into correct orthant
   [~,mii] = sort(mi, 'ascend');
   w = w(mii);
   w = w.*sign(w_star);
  end
  
  function [Y,W,d,flag] = PCA(X,k)
   % PCA computes the k principal components with largest associated
   % eigenvalue.
   %
   % INPUTS
   % X: (n x d) zero mean data matrix (samples are rows)
   % k: number of principal components
   %
   % OUTPUTS
   % Y: (n x k) first k principal components
   % W: (d x k) respective loadings of PCs
   % d: eigenvalues in descending order
   % flag: convergence flag (see eigs)
   %
   % Copyright 2013 Christian Sigg (christian@sigg-iten.ch)
   % See license file for details.
   %
   % This code is no longer actively maintained. See my 'nsprcomp' R package
   % for a better documented and more feature complete implementation.
   
   % eigs options
   options.issym = 1;
   options.disp = 0;
   
   r = rank(X);
   if nargin < 2
    k = r;
   end
   
   if(size(X,1)>=size(X,2))
    if (k == r)
     [W,D] = eig(X'*X);
     flag = 1;
    else
     [W,D,flag] = eigs(X'*X,k,'lm',options);
    end
    [d,indx] = sort(diag(D),'descend');
    Y = X*W;
   else
    if (k == r)
     [W,D] = eig(X*X');
     flag = 1;
    else
     [W,D,flag] = eigs(X*X',k,'lm',options);
    end
    [d,indx] = sort(diag(D),'descend');
    W = X'*W;
    for j=1:size(W,2)
     W(:,j) = W(:,j)/norm(W(:,j));
    end
    Y = X*W;
   end
   
   Y = Y(:,indx(1:k));
   W = W(:,indx(1:k));
   d = d(1:k);
  end
  
  function w = emPCA_L1_NN(X, n_nz)
   % Expectation maximization for non-negative sparse PCA. Profits from
   % random restarts.
   %
   % INPUTS
   % X: (n x d) zero mean data matrix (samples are rows)
   % n_nz: number of nonzero coefficients in w (this is only an upper limit)
   %
   % OUTPUTS
   % w: principal component (normalized to unit length)
   %
   % Tuneable parameters are set in the global variable PRM.
   %
   % Copyright 2013 Christian Sigg (christian@sigg-iten.ch)
   % See license file for details.
   %
   % This code is no longer actively maintained. See my 'nsprcomp' R package
   % for a better documented and more feature complete implementation.
   
   epsilon_ = NSPCA.L1_NN_eps;        % iteration stopping criterion
   optimize_ = NSPCA.L1_NN_optimize;  % perform optimization for given sparsity pattern
   
   [~,d] = size(X);
   
   w = NSPCA.iter(X, n_nz, epsilon_);
   if (optimize_ && n_nz < d)
    indx = w > 0;
    v = NSPCA.iter(X(:,indx), sum(indx), epsilon_);
    w = zeros(d,1);
    w(indx) = v;
   end
  end
    
  function w = iter(X, n_nz, eps)
   
   [~,d] = size(X);
   
   % initialize with random unit vector in nonnegative orthant
   w = abs(randn(d,1));
   w = w/norm(w);
   
   w_old = zeros(size(w));
   while(abs(w'*w_old) < 1 - eps)
    w_old = w;
    
    % E step: orthogonally project on w
    y = X*w;
    
    % M step:
    w = NSPCA.isoQP_L1_NN(y'*y, -(sum(X.*repmat(y,1,d)))', n_nz);
    w = w/norm(w);
   end
  end

  function w = isoQP_L1_NN(h, f, n_nz)
   % Solve min_w h/2*w'*w + w'*f, subject to norm(w,1) <= k and w_i >= 0.
   % Instead of specifying k, the number of nonzero elements n_nz is the
   % stopping criterion.
   %
   % The optimal solution sets w_i = 0 for f_i < 0 and calls isoQP_L1 (without
   % the non-negativity constraint) on the remaining variables.
   %
   % INPUTS
   % h: scalar of the Hessian
   % f: linear term
   % n_nz: number of nonzero coefficients in w (this is only an upper limit)
   %
   % OUTPUTS
   % w: optimum vector satisfying constraints
   %
   % Copyright 2013 Christian Sigg (christian@sigg-iten.ch)
   % See license file for details.
   %
   % This code is no longer actively maintained. See my 'nsprcomp' R package
   % for a better documented and more feature complete implementation.
   
   d = numel(f);
   nni = f<0;
   if (sum(nni) == 0)
    error('no non-negative components left');
   end
   
   w = zeros(d,1);
   w(nni) = NSPCA.isoQP_L1(h, f(nni), n_nz);
  end
  
 end
end