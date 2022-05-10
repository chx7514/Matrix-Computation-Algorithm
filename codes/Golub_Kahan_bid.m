function [U,V,a,b]=Golub_Kahan_bid(A,v)
    [m,n]=size(A);
    U=zeros(m);
    V=zeros(n);
    a=zeros(n,1);
    b=zeros(n-1,1);
    tol=1e-15;
    %k=1
    k=1;
    v=v/norm(v);
    V(:,1)=v;
    r=A*V(:,1);
    a(1)=norm(r);
    U(:,1)=r/a(1);
    p=A'*U(:,1)-V(:,1)*a(1);
    b(1)=norm(p);
    while k<n-1 && b(k)>tol
        V(:,k+1)=p/b(k);
        k=k+1;
        r=A*V(:,k)-b(k-1)*U(:,k-1);
        a(k)=norm(r);
        U(:,k)=r/a(k);
        p=A'*U(:,k)-a(k)*V(:,k);         
        b(k)=norm(p);
    end
    if b(k)>tol
        V(:,k+1)=p/b(k);
        k=k+1;
        r=A*V(:,k)-b(k-1)*U(:,k-1);
        a(k)=norm(r);
        U(:,k)=r/a(k);
        if norm(A'*U(:,k)-a(n)*V(:,k))>tol
            b(k-1)=-b(k-1);
            V(:,k)=p/b(k-1);
            r=A*V(:,k)-b(k-1)*U(:,k-1);
            a(k)=norm(r);
            U(:,k)=r/a(k);
        end
    end