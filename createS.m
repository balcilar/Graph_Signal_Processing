function Sx=createS(u,s,basis)
n=size(u,1);
nf=size(s,2);
K=size(basis,2);

Sx=zeros(size(basis,1),nf*K);
if nf>1
    for i =1:size(s,2)
        Sx(:,(i-1)*K+1:i*K)=createS(u,s(:,i),basis);
    end
    return
end



U=zeros(n*n,n);
S=zeros(n,n*n);
it=0;
for i=1:n:n*n
    it=it+1;
    for j=1:n
        U(i:i+n-1,j)=u(:,j)*u(it,j);
    end
    S(it,i:i+n-1)=s';
end

Sx=S*U*basis;