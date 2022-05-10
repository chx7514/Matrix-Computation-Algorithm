function [x,times]=Gauss(A,b)
    flag=1;
    times=0;
    tol=1e-5;
    n=size(A,1);
    x=rand(n,1);
    DL=tril(A);
    M=eye(n)-DL\A;
    g=DL\b;
    while (flag)
        x=M*x+g;
        R=abs(A*x-b);
        times=times+1;
        if all(R<tol)
            flag=0;
        end
    end