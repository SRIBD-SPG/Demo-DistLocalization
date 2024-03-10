%% eigen value decomposition of first r principle eigenvalues in descending order
%  author: Rui Zhou
function [Vs, Ds] = eigDesOrder(S, r)
[V, D] = eig(S);
[d,ind] = sort(diag(D), 'descend');
Ds = D(ind,ind);
Vs = V(:,ind);

Ds = Ds(1:r, 1:r);
Vs = Vs(:, 1:r);
