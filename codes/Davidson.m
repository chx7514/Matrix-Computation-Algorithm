% Davidson Algorithm
function [x,lambda]=Davidson(A,m) % m is the number of the eigenpairs to get
kmax=25;
n=size(A,1);
V=zeros(n,kmax*m);
E=[];
V(:,1:m)=orth(rand(n,m));
t=zeros(n,m);
H=V(:,1:m)'*A*V(:,1:m);
H=(H+H')/2;
[Eig_vecs,lambda]=eigs(H,m,'smallestabs');
x=V(:,1:m)*Eig_vecs;
r=A*x-x*lambda;
E=[E,norm(r,'fro')];
tol=1e-6;
D=diag(A);
times=0;
flag=0;
tic
while 1
    while 1
        for k=2:kmax
            times=times+1;
            for i=1:m
                t(:,i)=-r(:,i)./(D-lambda(i,i));
            end
            V(:,(k-1)*m+1:k*m)=t;
            for i=1:k*m
                V(:,i)=V(:,i)/norm(V(:,i));
                for j=i+1:k*m
                    V(:,j)=V(:,j)-V(:,i)'*V(:,j)*V(:,i);
                end
            end
            H=V(:,1:k*m)'*A*V(:,1:k*m);
            H=(H+H')/2;
            [s,lambda]=eigs(H,m,'smallestabs');
            x=V(:,1:k*m)*s;
            r=A*x-x*lambda;
            E=[E,norm(r,'fro')];
            if norm(r,'fro')<tol
                flag=1;
                break;
            end
        end
        if flag
            break;
        end
        V(:,1:m)=x;
        if times>500
            break;
        end
    end
end
    times
    semilogy(E);
    title('convergence history');
    xlabel('times');
    tlabel('lg(||r||)');
end