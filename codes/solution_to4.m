% solution to poison equation
function []=solution_to4(n)
    T=2*eye(n)+(-1)*diag(ones(n-1,1),1)+(-1)*diag(ones(n-1,1),-1);
    I=eye(n);
    A=kron(T,I)+kron(I,T);
    b=zeros(n*n,1);
    for i=1:n
        b(i)=sin(i*pi/(n+1));
    end
    for i=1:n
        x(n*(i-1)+1:n*i,1)=i*ones(n,1);
    end
    x=x/n;
    y0=[1:n]';
    for i=1:n
        y(n*(i-1)+1:n*i,1)=y0;
    end
    y=y/n;
    z=CG(A,b,zeros(n*n,1), 1000);
    figure(2);
    scatter3(x,y,z,'.');
    title('solution by CG')