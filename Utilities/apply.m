%% execute a cell of functions on a given dataset
% Author: Rui ZHOU

function res = apply(funs, Xs,r,gamma)

if (~iscell(funs))
    funs = {funs};
end
if (~iscell(Xs))
    Xs = {Xs};
end

n_f = length(funs);
n_X = length(Xs);

res = zeros(n_f, n_X);
for i = 1:n_f
    for j = 1:n_X
        res(i, j) = funs{i}(Xs{j},r,gamma);
    end
end

end
