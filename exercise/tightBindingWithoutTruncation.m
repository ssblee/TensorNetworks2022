clear

constants.N = 5;
constants.t = ones(1, constants.N - 1);

[F, Z, S, I] = getLocalSpace('FermionS');

dummy = 1;
H_prev = dummy;
I_prev = dummy;
H_next = I;
H = zeros(size(H_next, 1));
ZF = contract(Z, 2, 2, F, 3, 1);

for i = 2:constants.N
    %  STEP 1 Expand I_prev/next
    %   ╔════════H_next════════╗
    %  ╔══H_prev══╗
    %   ┌────────   ──I_prev──   ──I_next──
    %   |               |             |
    %   |
    %   |               |             |
    %   └────────   ──I_prev──   ──I_next──
    %  ╚══════════╝
    %   ╚══════════════════════╝
    %  STEP 2 Update Hamiltonian
    %   ╔═════════════H_hopping═════════════╗
    %  ╔════════H_accum════════╗
    %   ┌────────   ──I_prev──   ──I_next──
    %   |               |             |
    %   |               F─────   ────ZF
    %   |               |             |
    %   └────────   ──I_prev──   ──I_next──
    %  ╚═══════════════════════╝
    %   ╚═══════════════════════════════════╝
    %  STEP 3 update H_prev/next
    %  ╔═══════════════H_next═══════════════╗
    %   ╔════════H_prev════════╗
    %   ┌────────   ──I_prev──   ──I_next──
    %   |               |             |
    %   |
    %   |               |             |
    %   └────────   ──I_prev──   ──I_next──
    %   ╚══════════════════════╝
    %  ╚════════════════════════════════════╝
    I_prev = getIdentity(H_prev, 2, I, 2, [1 3 2]);
    I_next = getIdentity(H_next, 2, I, 2, [1 3 2]);

    H_accum = updateL(H_prev, I_prev, F, I_prev);
    H_hopping = updateL(H_accum, I_next, dagger(ZF), I_next);

    H_prev = H_next;
    H_next = updateL(H_next, I_next, [], I_next);

    H = updateL(H, I_next, [], I_next);
    H = H + H_hopping + dagger(H_hopping);
end

mine = ground(H);

% Compare with nonIntTB.m
answer = nonIntTB(constants.t);
gwayeon(mine, 2 * answer);
