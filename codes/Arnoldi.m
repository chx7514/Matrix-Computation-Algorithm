% Arnoldi procedure
function [Q,H,t]=Arnoldi(A,q)
    n=size(A,1);
    k=0;
    r=q;
    tol=1e-16;
    Q=zeros(n,1);
    H=[];
    while k<n && (k==0 || abs(H(k+1,k))>tol)
        if k==0
            Q(:,1)=r/norm(r);
        else
            Q=[Q,r/H(k+1,k)];
        end
        k=k+1;
        r=A*Q(:,k);
        H=[H,zeros(k,1)];
        for i=1:k
            H(i,k)=Q(:,i)'*r;
            r=r-H(i,k)*Q(:,i);
        end 
        H(k+1,k)=norm(r);
    end
    t=k;
    H=H(1:n,:);