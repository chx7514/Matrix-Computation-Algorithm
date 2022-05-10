%% Function nonsymmetric_QR
%input:an n*n matrix
%output: T(real schur form of A), an orthogonal matrix Q that Q'*A*Q=T, interation number
function [A,Q,time]=nonsymmetric_QR(A)
    n=size(A,1);
    [A,Q]=Hessenberg(A);
    u=1e-16;
    time=0;
    %partition the matrix
    while 1
        for i=2:n
            if abs(A(i,i-1))<=(abs(A(i-1,i-1))+abs(A(i,i)))*u
                A(i,i-1)=0;
            end
        end % set to zero
        m=n+1;
        for i=n:-1:2
            if A(i,i-1)==0 %then is 1*1
                m=i; 
                if i==2
                    m=1;
                end 
            elseif i==2 || A(i-1,i-2)==0 %then is 2*2
                b=trace(A(i-1:i,i-1:i));
                c=det(A(i-1:i,i-1:i));
                if b*b-4*c < 0
                    m=i-1;
                else
                    break
                end 
            else
               break
            end
        end %get m
        t=m;
        for i=m-1:-1:2
            if A(i,i-1)~=0
                t=i-1;
            else
                break
            end
        end %get t
        t=t-1;
        m=n-m+1;
        if m==n
            return
        end
        [A(t+1:n-m,t+1:n-m),P]=QR_Convergence(A(t+1:n-m,t+1:n-m)); %perform a Francis QR step on H22
        Q=Q*blkdiag(eye(t),P,eye(m)); %calculate the orthogonal matrix Q
        A(1:t,t+1:n-m)=A(1:t,t+1:n-m)*P; 
        A(t+1:n-m,n-m+1:n)=P'*A(t+1:n-m,n-m+1:n); %calculate other blocks of A
        time=time+1;
    end
end
%% Function Hessenberg
%input: an n*n matrix
%output: H(a Hessenberg matrix from A), an orthogonal matrix Q that Q'*A*Q=H
function [A,Q]=Hessenberg(A)
    n=size(A,1);
    Q=eye(n);
    for i=1:n-2
        [v,b]=House(A(i+1:n,i));
        H=b*(v*v');
        T=eye(n-i)-H;
        Q=Q*blkdiag(eye(i),T);
        A(i+1:n,i:n)=A(i+1:n,i:n)-H*A(i+1:n,i:n);
        A(1:n,i+1:n)=A(1:n,i+1:n)-A(1:n,i+1:n)*H;
    end
end
%% Function QR_Convergence
%input: an unreduced Hessenberg matrix H
%output: an matrix after one Francis QR step on H, an orthogonal matrix Q that Q'*H*Q=H1
function [H,P]=QR_Convergence(H)
    n=size(H,1);
    if n>2 %perform Francis QR step
        m=n-1;
        s=H(m,m)+H(n,n);
        t=H(m,m)*H(n,n)-H(m,n)*H(n,m);
        x=H(1,1)^2+H(1,2)*H(2,1)-s*H(1,1)+t;
        y=H(2,1)*(H(1,1)+H(2,2)-s);
        z=H(2,1)*H(3,2);
        P=eye(n);
        for i=0:n-3
            [v,b]=House([x,y,z]');
            T=b*(v*v');
            P=P*blkdiag(eye(i),eye(3)-T,eye(n-i-3));
            p=max(1,i);
            H(i+1:i+3,p:n)=H(i+1:i+3,p:n)-T*H(i+1:i+3,p:n);
            q=min(n,i+4);
            H(1:q,i+1:i+3)=H(1:q,i+1:i+3)-H(1:q,i+1:i+3)*T;
            x=H(i+2,i+1);
            y=H(i+3,i+1);
            if i<n-3
                z=H(i+4,i+1);
            end
        end
        [v,b]=House([x,y]');
        T=b*(v*v');
        P=P*blkdiag(eye(n-2),eye(2)-T);
        H(n-1:n,n-2:n)=H(n-1:n,n-2:n)-T*H(n-1:n,n-2:n);
        H(1:n,n-1:n)=H(1:n,n-1:n)-H(1:n,n-1:n)*T;
    else % if H is 2*2 and has real eigenvalue, use single shift
        t=H(2,2);
        [c,s]=Givens(H(1,1)-t,H(2,1));
        H=[c,s;-s,c]*H*[c,-s;s,c];
        P=[c,-s;s,c];
    end
end
%% Function House
%perform Householder transformation
function [v,b]= House (x)
	n=length(x);
	x=x/norm(x,inf);
	sigma=x(2:n)'*x(2:n);
    v=zeros(n,1);
	v(2:n)=x(2:n);
	if sigma==0
		b=0;
	else
		a=sqrt((x(1))^2+sigma); 
		if x(1)<=0
			v(1)=x(1)-a;
		else
			v(1)=-sigma/(x(1)+a);
        end
		b=2*(v(1)^2)/(sigma+(v(1)^2));
		v=v/v(1);
    end
end
%% Function Givens
%perform Givens transformation
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
end
%% Function Get_eigenpairs
%input: a real schur form matrix A
%output: all eigenvalues of A
function [Q,x]=Get_eigenpairs(T,A)
    n=size(T,1);
    x=zeros(n,1);
    i=1;
    Q=zeros(n);
    R=zeros(n,1);
    while 1
        if T(i+1,i)==0
            x(i)=T(i,i);
            Q(:,i)=real_eigen(A,x(i));
            R(i)=norm(A*Q(:,i)-Q(:,i)*x(i));
            i=i+1;
        else 
            b=trace(T(i:i+1,i:i+1));
            c=det(T(i:i+1,i:i+1));
            delta=b*b-4*c;
            x(i)=(b+sqrt(delta))/2;
            Q(:,i)=im_eigen(A,x(i));
            R(i)=norm(A*Q(:,i)-Q(:,i)*x(i));
            R(i+1)=R(i);
            Q(:,i+1)=conj(Q(:,i));
            x(i+1)=(b-sqrt(delta))/2;
            i=i+2;
        end
        if i>=n
            break
        end
    end
    if i==n
        x(n)=T(n,n);
        Q(:,i)=real_eigen(A,x(i));
        R(i)=norm(A*Q(:,i)-Q(:,i)*x(i));
    end
    semilogy(1:n,R);
end