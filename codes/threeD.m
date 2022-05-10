% 求解三维谐振子系统
function [y,lambda]=threeD(m,omega,N,r,k)
    dx = 2 * r / (N-1);
    x = linspace(-r, r, N);
    D = sparse(1:N, 1:N, 2*ones(1, N), N, N);
    C = sparse(2:N, 1:N-1, -1*ones(1, N-1), N, N);
    T = D + C + C';
    I = sparse(1:N, 1:N, ones(1, N), N, N); 
    B = kron(T, I);
    A = B + kron(I, T);
    A = kron(B, I) + kron(I, A);
    A = A / (2 * m * dx^2);
    d1 = zeros(N^2, 1);
    for i = 1 : N
        d1(i * N - N + 1 : i * N) = x(i)^2 * ones(1, N)+ x.^2;
    end
    d = zeros(N^3, 1);
    n = N^2;
    for i = 1 : N
        d(i * n - n + 1 : i * n) = d1 + x(i)^2;
    end
    d = m * omega^2 * d / 2;
    N = N^3;
    A = A + sparse(1:N, 1:N, d, N, N);
    [y,lambda]=Davidson(A,k);
end