function ET = getinfo(T,E)
L = size(E,1);
N = max(max(T));
A = sparse(E(:,1),E(:,2),(1:L)',N,N);
ET = zeros(N,3);
i = 1;
for tri=T'
    Eloc = full(A(tri,tri)); Eloc = Eloc + Eloc';
    ET(i,:) = Eloc([8 7 4]);
    i = i + 1;
end
end
