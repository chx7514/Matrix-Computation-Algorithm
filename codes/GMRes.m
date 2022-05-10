% GMRES algorithm
function [x]=GMRes(A,b,x0,m)
    k=0;
    n=size(A,1);
    r=b-A*x0;
    tol=1e-16;
    Q=zeros(n,1);
    H=[];
    G=[];
    tmp=norm(r);
    p=tmp*[1;zeros(m-1,1)];
    res=[]; %保存每一次的r
    m=min(m,n); 
    while k<m && tmp>tol
        if k==0
            Q(:,1)=r/tmp;
        else
            Q=[Q,r/tmp];
        end
        k=k+1;
        r=A*Q(:,k);
        H=[H,zeros(m,1)];
        for i=1:k
            H(i,k)=Q(:,i)'*r;
            r=r-H(i,k)*Q(:,i);
        end 
        tmp=norm(r);
        H(k+1,k)=tmp;
        for i=1:k-1
            H(i:i+1,k)=[G(1,i),G(2,i);-G(2,i),G(1,i)]*H(i:i+1,k);
        end
        if k==m
            break
        end
        [c,s]=Givens(H(k,k),H(k+1,k));
        G=[G,[c;s]];
        H(k:k+1,k)=[c,s;-s,c]*H(k:k+1,k);
        p(k:k+1)=[c,s;-s,c]*p(k:k+1);
        rk=p;
        for i=k:-1:2
            rk(i)=rk(i)/H(i,i);
            rk(1:i-1)=rk(1:i-1)-H(1:i-1,i)*rk(i);
        end
        rk(1)=rk(1)/H(1,1);
        x=x0+Q(:,1:k)*rk(1:k);
        res=[res;norm(A*x-b)];
    end   
    for i=k:-1:2
        p(i)=p(i)/H(i,i);
        p(1:i-1)=p(1:i-1)-H(1:i-1,i)*p(i);
    end
    p(1)=p(1)/H(1,1);
    x=x0+Q(:,1:k)*p(1:k);
    res=[res;norm(A*x-b)];
    time=[1:k]';
    plot(time,res);
    
function [c,s]=Givens(a,b)
    if b==0
        c=1;
        s=0;
    else
        if abs(b)>abs(a)
            t=a/b;
            s=1/sqrt(1+t^2);
            c=s*t;
        else
            t=b/a;
            c=1/sqrt(1+t^2);
            s=c*t;
        end
    end