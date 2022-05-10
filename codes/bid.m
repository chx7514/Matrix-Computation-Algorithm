%  Input: an m*n matrix A, while m>=n
%  Output: U,B,V; U and V are orthogonal, B is bidiagonal,
function [U,A,V]=bid(A)
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