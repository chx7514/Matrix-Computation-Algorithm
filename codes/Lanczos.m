function []=Lanczos()
    n=1000;
    kmax=10;
%     A=rand(n);
%     A=A+A';
%     d=sort(rand(n,1));
    d=zeros(n,1);
    d(1:n-1)=5*ones(n-1,1);
    d(n)=10;
    d(n-1)=d(n-1)+1e-3;
    D=diag(d);
    U=orth(rand(n));
    A=U*D*U';
    q=rand(n,1);
    Lanczos_for_eigen(A,q,kmax);