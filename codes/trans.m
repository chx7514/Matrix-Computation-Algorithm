function [U]=trans(x,n)
    k=1;
    U=zeros(n,n);
    for i=1:n
        U(i,1:n)=x(k:k+n-1);
        k=k+n;
    end