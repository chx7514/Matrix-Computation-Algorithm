function [x]=real_eigen(A,l)
    n=size(A,1);
    x=rand(n,1);
    k=0;
    B=A-l*eye(n);
    while k<2  
        z=B\x;
        x=z/norm(z);
        k=k+1;
    end