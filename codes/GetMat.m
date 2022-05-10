function [A]=GetMat(a,b,c,n)
    A(1,1:2)=[b,c];
    for i = 2 : n - 1
        A(i, i - 1 : i + 1) = [a, b, c];
    end
    A(n, n-1:n)=[a,b];
 