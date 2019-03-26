function X=chebyshev_basis(L, x, degree)
X=x;
X(:,2)=L*x;
for i=2:degree
    X(:,i+1)=2*L*X(:,i)-X(:,i-1);
end