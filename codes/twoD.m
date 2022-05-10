% 求解二维谐振子系统
function [y,lambda]=twoD(m,omega,N,r)
    n=N^2;
    dx = 2 * r / (N-1);
    x = linspace(-r, r, N);
    D = sparse(1:N, 1:N, 2*ones(1, N), N, N);
    E = sparse(2:N, 1:N-1, -1*ones(1, N-1), N, N);
    T = D + E + E';
    I = sparse(1:N, 1:N, ones(1, N), N, N);
    A = kron(T, I) + kron(I, T);
    A = A / (2 * m * dx^2);
    d = zeros(n, 1);
    for i = 1 : N
        d(i * N - N + 1:i * N) = x(i)^2 * ones(1, N)+ x.^2;
    end
    d = m * omega^2 * d / 2;
    A = A + sparse(1:n, 1:n, d, n, n);
    [y,lambda]=Davidson(A,k);