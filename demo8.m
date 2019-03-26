clear all
close all
clc

% create time series grid as graph
n=101;
W=zeros(n,n);
for i=1:n-1  
    
    W(i,i+1)=1;
    W(i+1,i)=1;  
    
end

W(1,n)=1;
W(n,1)=1;


% calculate combinatorial Laplacian Matrix

d = sum(W,2);
L = diag(d)-W;


% calculate basis
[u v]=eig(L);


% make eignevalue as vector
v=diag(v);
v(v<0)=0;
% get maximum eigenvalue
lmax=max(v);

nv=linspace(0,2,21);
basis=bspline_basis(21, nv,v, 3);

nu=u*basis;


figure;subplot(2,3,1);plot(u(:,5));ylim([-0.15 0.15]);xlim([1 101]);
title(['5th eigenvector, value=' num2str(v(5))]);
subplot(2,3,2);plot(u(:,11));ylim([-0.15 0.15]);xlim([1 101]);
title(['11th eigenvector, value=' num2str(v(11))]);
subplot(2,3,3);plot(u(:,end));ylim([-0.15 0.15]);xlim([1 101]);
title(['last eigenvector, value=' num2str(v(end))]);

subplot(2,3,4);plot(nu(:,5));xlim([1 101]); %ylim([-0.15 0.15]);
title(['5th eigenvector, value=' num2str(nv(5))]);
subplot(2,3,5);plot(nu(:,11));xlim([1 101]); %ylim([-0.15 0.15]);
title(['11th eigenvector, value=' num2str(nv(11))]);
subplot(2,3,6);plot(nu(:,end));xlim([1 101]); %ylim([-0.15 0.15]);
title(['last eigenvector, value=' num2str(nv(end))]);
