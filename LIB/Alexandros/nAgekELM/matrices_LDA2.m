function [Slda_w, Slda_b] = matrices_LDA2(data, lbls)

noOfClasses = length(unique(lbls));  [n,N] = size(data);
Slda_w = zeros(n);      Slda_b = zeros(n);
mu = mean(data,2);
for k=1:noOfClasses
    curr_ind = find(lbls==k); Ni = length(curr_ind); 
    curr_mu = mean(data(:,curr_ind),2);
    diff1 = curr_mu*ones(1,Ni) - data(:,curr_ind);
    diff2 = curr_mu - mu;
    Slda_w = Slda_w + (diff1*diff1');
    Slda_b = Slda_b + Ni*(diff2*diff2');
end

Slda_w = Slda_w / N;
Slda_b = Slda_b / N;