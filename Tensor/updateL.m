function result = updateL(C, B, X, A)

    if (isempty(X))
        % assume B's 3rd leg has same dimension with A's 3rd
        X = eye(size(B, 3));
        % No error catch now
    end

    % same for C
    if (isempty(C))
        C = eye(size(B, 1));
        % No error catch now
    end

    % contract X and A first
    XA = contract(X, ndims(X), 2, A, ndims(A), 3);

    % contract C and XA and so on
    CX_bondLegCount = max(0, min(ndims(C) - 2, ndims(X) - 2));
    % [2] for C-A bond, [3 ~ 2 + CX_bondLegCount] for C-X bond
    legsC = [2, [3:3 - 1 + CX_bondLegCount]];
    % [ndims(X)] for C-A bond, [2 ~ 1 + CX_bondLegCount] for C-X bond
    legsXA = [ndims(X), [2:2 - 1 + CX_bondLegCount]];
    CXA = contract(C, ndims(C), legsC, XA, ndims(XA), legsXA);

    B_c = conj(B);
    legC = 1 + ((ndims(C) - 2) - CX_bondLegCount) + 1;
    BCXA = contract(B_c, ndims(B_c), [1, 3], CXA, ndims(CXA), [1, legC]);
    result = permute(BCXA, [1, ndims(BCXA), [2:ndims(BCXA) - 1]]);
end
