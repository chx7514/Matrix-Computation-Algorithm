% 求解 min||Ax-b||, s.t. Cx=d
function [x]=solve2(A,b,C,d)
    [m,n]=size(A); %A:m*n
    p=size(C,1); %C:p*n
    C=C'; %C':n*p
    s=zeros(p,1); %储存Householder变换
    for j=1:p
        if j<n
            [v,r]=House(C(j:n,j));
            H=r*(v*v');
            C(j:n,j:p)=C(j:n,j:p)-H*C(j:n,j:p);
            A(:,j:n)=A(:,j:n)-A(:,j:n)*H;
            s(j)=r;
            C(j+1:n,j)=v(2:n-j+1);
        end
    end % 对C'进行QR分解 同时A变成[AQ1,AQ2]
    for i=1:p-1
        d(i)=d(i)/C(i,i);
        d(i+1:p)=d(i+1:p)-d(i)*C(i,i+1:p)';
    end
    d(p)=d(p)/C(p,p); % d=L^-1*d
    b=b-A(:,1:p)*d; % b=b-A*Q1*L^-1*d
    A=A(:,p+1:n); % A=A*Q2
    n=n-p;
    for j=1:n
        if j<m
            [v,r]=House(A(j:m,j));
            H=r*(v*v');
            A(j:m,j:n)=A(j:m,j:n)-H*A(j:m,j:n);
            b(j:m)=b(j:m)-H*b(j:m);
        end
    end % 对AQ2进行QR分解
    for i=n:-1:2
        b(i)=b(i)/A(i,i);
        b(1:i-1)=b(1:i-1)-A(1:i-1,i)*b(i);
    end
    b(1)=b(1)/A(1,1);
    x=[d;b(1:n)];
    n=n+p;
    for j=p:-1:1
        v=zeros(n-j+1,1);
        v(2:n-j+1)=C(j+1:n,j);
        v(1)=1;
        x(j:n)=x(j:n)-s(j)*(v*v')*x(j:n);
    end % x=Qx