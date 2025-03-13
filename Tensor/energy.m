function result = energy(hamiltonian)
    result = sort(eig(Hermitian(hamiltonian)));
end
