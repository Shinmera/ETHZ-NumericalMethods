function [E,Eb] = processmesh(T)
    N = max(max(T)); 
    M = size(T,1);   
    T = sort(T,2);
    C = [T(:,1) T(:,2); T(:,2) T(:,3); T(:,1) T(:,3)];
    A = sparse(C(:,1), C(:,2), ones(3*M,1), N, N);
    [I,J] = find(A > 0); E = [I,J];
    [I,J] = find(A == 1); Eb = [I,J];
end
