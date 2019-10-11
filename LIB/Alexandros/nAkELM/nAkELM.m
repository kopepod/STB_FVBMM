% Nystrom-based approximate kernel ELM
function [labels, A,Ot] = nAkELM(train_data, test_data, train_lbls, C, n, clst_bool)
% train_data --> DxN matrix (N training samples)
% test_data --> DxNt matrix (Nt test samples)
% train_lbls --> Nx1 training label vector
% C --> regularization parameter
% n --> number of reference vector (for the reduced kernel)
% clst_bool --> boolean for using clustering approach or randomized

calc_alpha = 1;

if nargin<6,  clst_bool = 1;  end
N = size(train_data,2);  M = size(test_data,2);

n = round(n * N);

% calculate target vectors
noOfClasses = length(unique(train_lbls));
%T = -ones(noOfClasses,noOfTrainingData);    
T = -ones(noOfClasses,size(train_data,2));
for i=1:N, T(train_lbls(i),i) = 1.0; end

% use single precision
%train_data = single(train_data);  test_data = single(test_data);  T = single(T);

%%%%% training
if clst_bool~=1
    pp = randperm(N);  ref_data = train_data(:,pp(1:n));
else
    options.max_iter = 5;
    if N > 100000
        pp = randperm(N);        
        [tmp_vec, ref_data] = fkmeans(double(train_data(:,pp(1:50000))'), n, options);
    else
        [tmp_vec, ref_data] = fkmeans(double(train_data'), n, options.max_iter);
    end
    clear tmp_vec;    ref_data = ref_data';
    % use single precision
    %ref_data = single(ref_data);
end

K_Nn = ((sum(train_data'.^2,2)*ones(1,n))+(sum(ref_data'.^2,2)*ones(1,N))'-(2*(train_data'*ref_data)));
K_nn = ((sum(ref_data'.^2,2)*ones(1,n))+(sum(ref_data'.^2,2)*ones(1,n))'-(2*(ref_data'*ref_data)));
if calc_alpha==1,  A = 2 * mean(mean([K_Nn; K_nn]));  end
K_Nn = exp(-K_Nn/A);   [V,D] = eig(exp(-K_nn/A));

K_Nn = (K_Nn * V * D^(-0.5))';
W = pinv( K_Nn*K_Nn' + eye(n)/C ) *K_Nn * T';

%%%%% testing
if N > 100000
    Ktest = [];
    start_ind = 1;  NN = 10000;  
    while start_ind<N
        end_ind = start_ind + NN -1;
        if end_ind > N,  end_ind = N;  end
        MM = end_ind - start_ind+ 1;
        curr_Ktest = ((sum(train_data(:,start_ind:end_ind)'.^2,2)*ones(1,M))+(sum(test_data'.^2,2)*ones(1,MM))'-(2*(train_data(:,start_ind:end_ind)'*test_data)));
        start_ind = end_ind + 1;
        Ktest = [Ktest; curr_Ktest];
        disp(['end: ',num2str(end_ind)])
    end
    
else
    Ktest = ((sum(train_data'.^2,2)*ones(1,M))+(sum(test_data'.^2,2)*ones(1,N))'-(2*(train_data'*test_data)));
end
Ktest = exp(-Ktest/A);
clear train_data;  clear test_data;
Ktest = pinv(K_Nn') * Ktest;
Ot = W' * Ktest;  
[tmp_vec labels] = max(Ot);     labels = labels';
