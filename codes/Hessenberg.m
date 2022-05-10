%input: an n*n matrix
%output: H(a Hessenberg matrix from A), an orthogonal matrix Q that Q'*A*Q=H
function [A,Q]=Hessenberg(A)
    n=size(A,1);
    Q=eye(n);
    for i=1:n-2
        [v,b]=House(A(i+1:n,i));
        H=b*(v*v');
        T=eye(n-i)-H;
        Q=blkdiag(eye(i),T)*Q;
        A(i+1:n,i:n)=A(i+1:n,i:n)-H*A(i+1:n,i:n);
        A(1:n,i+1:n)=A(1:n,i+1:n)-A(1:n,i+1:n)*H;
    end
    Q=Q';
    