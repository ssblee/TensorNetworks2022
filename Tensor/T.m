% Just exchange first two tensor indecies
function result = T(matrices)
    result = permute(matrices, [2, 1, (3:ndims(matrices))]);
end
