% Lanczos algorithm for computing the matrix functional v^*f(A)v, where A is Hermitian
function [answer]=Lanczos_for_function(A,u,kmax)
% f(x)=exp(x)
    nor=norm(u)^2;
    Ans=zeros(kmax,1);
    n=size(u,1);
    T=zeros(n);
    q1=u/norm(u);
    w=A*q1;
    T(1,1)=q1'*w;
    r=w-q1*T(1,1);
    T(2,1)=norm(r);
    T(1,2)=T(2,1);
    k=1;
    F=T(1,1);
    Ans(k)=exp(F)*nor;
    tol=1e-16;
    while (k<kmax) && T(k+1,k)>tol && (k==1 || (k>1 && abs(Ans(k)-Ans(k-1))/Ans(k-1)>tol))
        k=k+1;
        q2=r/T(k,k-1);
        w=A*q2;
        T(k,k)=q2'*w;
        if k<n
            r=w-q2*T(k,k)-q1*T(k,k-1);
            T(k+1,k)=norm(r);
            T(k,k+1)=T(k+1,k);
        end
        F=T(1:k,1:k);
        [U,D]=eig(F);
        for i=1:k
            D(i,i)=exp(D(i,i));
        end
        F=U*D*U';
        Ans(k)=F(1,1)*nor;
        q1=q2;
    end
    plot([1:k]',Ans(1:k));
    answer=Ans(k);
end