%% Function mysvd
%  Input: an m*n matrix A
%  Output: U,D,V; U and V are orthogonal, D is diagonal,
%          and A=U*D*V'
function [U,A,V,t]=mysvd(A)
    [m,n]=size(A);
    % mysvd_1 can only solve the case that m>=n
    if m>=n
        [U,A,V,t]=mysvd_1(A);
    else
        [U,A,V,t]=mysvd_1(A');
        A=A';
        T=U;
        U=V;
        V=T;
    end    
end
%% Function mysvd_1
%  Input: an m*n matrix A, while m>=n
%  Output: U,D,V; U and V are orthogonal, D is diagonal,
%          and A=U*D*V'
function [U,A,V,time]=mysvd_1(A) 
    n=size(A,2);
    [U,A,V]=bidiagonalize(A); 
    % U'*A*V=T; T is bidiagonal
    tol=1e-12;
    time=0;
    while 1
        % set to zero
        norm_A=norm(A,'inf');
        for k=1:n
            if k<n && abs(A(k,k+1))<=(abs(A(k,k))+abs(A(k+1,k+1)))*tol
                A(k,k+1)=0;
            end
            if abs(A(k,k))<=norm_A*tol
                A(k,k)=0;
            end
        end
        % compute p and q, while A(p+1:n-q,p+1:n-q) has a non zero
        % superdiagnoal, A(n-q+1:n,n-q+1:n) is diagonal
        q=n+1;
        for i=n:-1:2
            if A(i-1,i)==0 
                q=i;
                if i==2
                    q=1;
                end
            else
                break
            end
        end
        p=q;
        for i=q-1:-1:2
            if A(i-1,i)~=0
                p=i-1;
            else
                break
            end
        end
        p=p-1;
        q=n-q+1;
        % A is diagonal
        if q==n
            break
        end
        time=time+1;
        % perform a Golub-Kahan SVD step on A(p+1:n-q,p+1:n-q)
        [U(:,p+1:n-q),A(p+1:n-q,p+1:n-q),V(:,p+1:n-q)]=SVD_conver(U(:,p+1:n-q),A(p+1:n-q,p+1:n-q),V(:,p+1:n-q));
    end
    % Get the common form of SVD
    [U,A,V]=common_form(U,A,V); 
end
%% Function bidiagonalize
%  Input: an m*n matrix A, while m>=n
%  Output: U,B,V; U and V are orthogonal, B is bidiagonal,
function [U,A,V]=bidiagonalize(A)
    [m,n]=size(A);
    U=eye(m);
    V=eye(n);
    for k=1:n
        if k<m
            % Householder transformation on column vectors
            [v,b]=House(A(k:m,k));
            T=b*(v*v');
            A(k:m,k:n)=A(k:m,k:n)-T*A(k:m,k:n);
            U(:,k:m)=U(:,k:m)-U(:,k:m)*T;
            if k<n-1
                % Householder transformation on row vectors
                [v,b]=House(A(k,k+1:n)');
                T=b*(v*v');
                A(k:m,k+1:n)=A(k:m,k+1:n)-A(k:m,k+1:n)*T;
                V(:,k+1:n)=V(:,k+1:n)-V(:,k+1:n)*T;
            end
        end
    end
    A=A(1:n,:);
end
%% Function SVD_conver
%  Input: a bidiagonal matrix B, which has a non zero superdiagnoal
%  Output: the result of a Golub-Kahan SVD step on B
function [P,B,Q]=SVD_conver(P,B,Q)
    n=size(B,1);
    if n==1 % in case B is 1*1
        return
    end
    for i=1:n-1
        if B(i,i)==0
            % perform a Golub-Kahan SVD step on B(1:i,1:i)
            [P(:,1:i),B(1:i,1:i),Q(:,1:i)]=SVD_conver(P(:,1:i),B(1:i,1:i),Q(:,1:i)); 
            % zero the superdiagonal entry in the ith row
            for j=i+1:n
                [c,s]=Givens(B(j,j),B(i,j));
                B(j,j)=c*B(j,j)+s*B(i,j);
                B(i,j)=0;
                Tmp=[P(:,i),P(:,j)]*[c,s;-s,c];
                P(:,i)=Tmp(:,1);
                P(:,j)=Tmp(:,2);
                if j<n
                    Tmp=[c,-s;s,c]*[B(i,j+1);B(j,j+1)];
                    B(i,j+1)=Tmp(1);
                    B(j,j+1)=Tmp(2);
                end
            end
            % perform a Golub-Kahan SVD step on B(i+1:n,i+1:n)
            [P(:,i+1:n),B(i+1:n,i+1:n),Q(:,i+1:n)]=SVD_conver(P(:,i+1:n),B(i+1:n,i+1:n),Q(:,i+1:n));
            return
        end
    end
    % another case is that B is irreducible
    a=B(n,n)*B(n,n)+B(n-1,n)*B(n-1,n);
    % calculate Wilkinson shift of B'*B
    if n>2
        d=(B(n-1,n-1)*B(n-1,n-1)+B(n-2,n-1)*B(n-2,n-1)-a)/2;
    else
        d=(B(1,1)*B(1,1)-a)/2;
    end
    b=B(n-1,n-1)*B(n-1,n);
    u=a-b*b/(d+sign(d)*sqrt(d*d+b*b));
    y=B(1,1)*B(1,1)-u;
    z=B(1,1)*B(1,2);
    for k=1:n-1
        [c,s]=Givens(y,z); % Get Q
        if k>1
            B(k-1,k)=c*y+s*z;
        end
        Tmp=B(k:k+1,k:k+1)*[c,-s;s,c];
        B(k:k+1,k+1)=Tmp(:,2);
        y=Tmp(1,1);
        z=Tmp(2,1);
        Q(:,k:k+1)=Q(:,k:k+1)*[c,-s;s,c];
        [c,s]=Givens(y,z); % Get P
        B(k,k)=c*y+s*z;
        if k<n-1
            Tmp=[c,s;-s,c]*B(k:k+1,k+1:k+2);
            B(k+1,k+1:k+2)=Tmp(2,:);
            y=Tmp(1,1);
            z=Tmp(1,2);
        else
            B(k:k+1,k+1)=[c,s;-s,c]*B(k:k+1,k+1);
        end
        P(:,k:k+1)=P(:,k:k+1)*[c,-s;s,c];
    end
end
%% Function common_form
%  Input: original matrix of SVD
%  Output: common form of SVD, while the sigular values are sorted
function [U,D,V]=common_form(U,A,V)
    m=size(U,1);
    n=size(V,1);
    % turn the negetive sigular values to the positive ones
    for i=1:n
        if A(i,i)<0
            A(i,i)=-A(i,i);
            U(:,i)=-U(:,i);
        end
    end
    % sort U,A,V by sigular values
    [~,I]=sort(diag(A),'descend');    
    for i=1:m
        D=U(i,1:n);
        D=D(I);
        U(i,1:n)=D;
    end
    for i=1:n
        D=V(i,:);
        D=D(I);
        V(i,:)=D;
    end
    D=zeros(n);
    for i=1:n
        D(i,i)=A(I(i),I(i));
    end
    D=[D;zeros(m-n,n)];
end
%% Function Givens
%  calculate a Givens transformation
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
%% Function House
%  calculate a Householder transformation
function [v,b]= House (x)
	n=length(x);
    v=zeros(n,1);
    if x==0
        b=0;
        return;
    end
	x=x/norm(x,inf);
	sigma=x(2:n)'*x(2:n);
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