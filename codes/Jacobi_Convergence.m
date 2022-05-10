function [x,times]=Jacobi_Convergence(A,b)
    Dn=diag(1./diag(A));
    B=eye(size(A))-Dn*A;
    f=Dn*b;
    x=rand(size(A,1),1);
    flag=1;
    tol=1e-5;
    times=0;
    while(flag)
        x=B*x+f;
        R=abs(A*x-b);
        times=times+1;
        if all(R<tol)
            flag=0;
        end
    end
    