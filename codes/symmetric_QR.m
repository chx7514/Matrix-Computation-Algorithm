% 对称矩阵QR算法求解特征值
function [A,Q,times]=symmetric_QR(A)
    n=size(A,1);
    tol=1e-16;
    times=0;
    Q=eye(n);
    for k=1:n-2
        [v,b]=House(A(k+1:n,k));
        Q=Q*blkdiag(eye(k),(eye(n-k)-b*(v*v')));
        u=b*A(k+1:n,k+1:n)*v;
        w=u-b*(v'*u)*v/2;
        A(k+1:n,k+1:n)=A(k+1:n,k+1:n)-v*w'-w*v';
        A(k+1,k)=norm(A(k+1:n,k));
        A(k,k+1)=A(k+1,k);
        A(k+2:n,k)=zeros(n-k-1,1);
        A(k,k+2:n)=zeros(1,n-k-1);
    end
    while 1
        for i=2:n
            if abs(A(i,i-1))<=(abs(A(i-1,i-1))+abs(A(i,i)))*tol
                A(i,i-1)=0;
            end
        end %置零
        m=n+1;
        for i=n:-1:2
            if A(i,i-1)==0
                m=i; 
                if i==2
                    m=1;
                end %是1*1
            else
               break
            end
        end %设置m
        t=m;
        for i=m-1:-1:2
            if A(i,i-1)~=0
                t=i-1;
            else
                break
            end
        end %设置l
        t=t-1;
        m=n-m+1;
        if m==n
            return
        end
        times=times+1;
        [A(t+1:n-m,t+1:n-m),P]=sym_QR_conv(A(t+1:n-m,t+1:n-m));
        Q=Q*blkdiag(eye(t),P,eye(m));
    end
    
function [A,Q]=sym_QR_conv(A)
    n=size(A,1);
    Q=eye(n);
    d=(A(n-1,n-1)-A(n,n))/2;
    u=A(n,n)+d-sign(d)*sqrt(d^2+A(n-1,n)^2);
    x=A(1,1)-u;
    z=A(2,1);
    for k=1:n-1
        [c,s]=Givens(x,z);
        G=Get_G(c,s,n,k);
        A=G*A*G';
        Q=Q*G';
        if k<n-1
            x=A(k+1,k);
            z=A(k+2,k);
        end
    end
        
function [G]=Get_G(c,s,n,k)
    G=eye(n);
    G(k,k)=c;
    G(k+1,k+1)=c;
    G(k+1,k)=-s;
    G(k,k+1)=s;
        