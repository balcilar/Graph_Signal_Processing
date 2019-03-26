clear all
close all
clc

% load predefined W matrix for 100 nodes
load mydata


% calculate combinatorial Laplacian Matrix
d = sum(W,2);
L = diag(d)-W;
% calculate  Laplacian Matrix

% find eigenvector and eigenvalues of combinatorial Laplacian
[u v]=eig(L);


% make eignevalue as vector
v=diag(v);
% get maximum eigenvalue
lmax=max(v);
v(v<0)=0;

% create signal where first node is 1 rest of them zero
s=zeros(size(W,1),1);
s(1)=1;

% determine filter
flt =exp(-20*v);
% apply that filter on to graph signal
sf=u*(flt.*(u'*s));



% filter over new eigenvalue basis
K=50;
nv=linspace(0,8,K)';
basis=bspline_basis(K, nv,v, 3);
alpha=exp(-20*nv);
flt=basis*alpha;


% apply that filter on to graph signal

sf2=u*diag(flt)*u'*s;

figure;plot(sf,'b--','linewidth',2)
hold on;plot(sf2,'r-')
xlabel('node id')
ylabel('node signal value')
legend({'original basis','spline basis'})

