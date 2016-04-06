function y = getit(A, x, k)
    [S,D] = eig(A);
    y = S*diag(diag(D).^k)* (S\x);
end
