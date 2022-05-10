% use random projection to compute full rank decomposition
function []=fullrank()
    r=10;
    m=400;
    n=400;
    t=50;
    B=rand(m,r);
    C=rand(n,r);
    A=B*C';
    omega=rand(t,m);
    T=omega*A;
    [~,~,V]=svd(T);
    V=V(:,1:r);
    U=A*V;
    R=U*V'-A;
    rank(U)
    rank(V)
    norm(R,'fro')
    imagesc(log10(abs(R)))
    colorbar
    title("log|U*V'-A|")
    