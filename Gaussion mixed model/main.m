%suppose the datasets are D1 and D2. Each row represents an observation or datapoint

D = readmatrix("clean_data_GSE39582_all_male.csv", NumHeaderLines=1);
%D2 = readtable("clean_data_GSE39582_all_female.csv");
D1 = D(1:86,11:419);
D2 = D(87:224, 11:419);
%D2= D2(:,11:419)
[n1,p] = size(D1);  
[n2,p] = size(D2);
D = [D1; D2];  
n = n1 + n2;
omega = n2/n;
mu1 = mean(D1, 1); 
mu2 = mean(D2, 1); 
delta = mu1 - mu2;
muhat = (mu1 + mu2)/2; %muhat is a p-dimensional row vector
ecov = cov(D);
while cond(ecov) > 1e+6
    ecov = ecov + 0.2*sqrt(log(p)/n)*diag(ones(1,p));
end
beta_init = ecov \ delta'; %beta_init = inv(ecov)*delta';

%labels = [ones(1,n1) repelem(2,n2)]; 
%labelshat = zeros(1,n);
%lam = 0.1*0.9.^(1:100);
%error_rate = zeros(1, 100);
%for k = 1:100
%    beta = clime(beta_init, ecov, delta', lam(k));
%    for i = 1:n
%        temp = (D(i,:) - muhat)*beta;
%        if  temp >= log(omega/(1-omega))
%            labelshat(i) = 1;
%        else
%            labelshat(i) = 2;
%        end
%    end
%    error_rate(k) = numel(find(labels ~= labelshat))/n;
%end
%[minim, I] = min(error_rate);
%opt_lam = 0.1*0.9.^I;

%opt_lam = 0.1*0.9.^30;
opt_lam = 0.005;
beta = clime(beta_init, ecov, delta', opt_lam);
sig_count = sum(abs(beta) > 0.2);
sig_value = maxk(abs(beta), sig_count);
[~, sig_index] = maxk(abs(beta), sig_count);

