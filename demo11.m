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

% determine filter by chebyshev
K=120;
nL=2*L/lmax-eye(size(L));
nu=chebyshev_basis(nL, s, K);
alpha=pinv(nu)*sf;


figure;plot(alpha)
xlabel('coeff id');
title('learned chebyshev filter coefficient')

% apply that filter on to graph signal
sf2=nu*alpha;


% visualize input and result
run gspbox/gsp_start

G=gsp_graph(W,coord);
figure;gsp_plot_signal(G,sf2)
title('Filtered signal on first Graph');


load data2
     

% calculate combinatorial Laplacian Matrix
d = sum(WW,2);
L = diag(d)-WW;
% calculate  Laplacian Matrix

% find eigenvector and eigenvalues of combinatorial Laplacian
[u v]=eig(L);


% make eignevalue as vector
v=diag(v);
% get maximum eigenvalue
lmax=max(v);
v(v<0)=0;

% create signal where first node is 1 rest of them zero
s=zeros(size(WW,1),1);
s(1)=1;

% determine filter
nL=2*L/lmax-eye(size(L));
nu=chebyshev_basis(nL, s, K);
% apply that filter on to graph signal
sf2=nu*alpha;

% determine filter
flt =exp(-20*v);
% apply that filter on to graph signal
sf=u*(flt.*(u'*s));


G=gsp_graph(WW,coord2);
figure;gsp_plot_signal(G,sf2)
title('Filtered signal on second graph by learned coeff');
figure;gsp_plot_signal(G,sf)
title('Filtered signal on second graph by standart filter');

figure;plot(sf2);hold on;plot(sf,'r-')
xlabel('node id')
legend({'filter result by learned coeff','filter result by standart filter'})


