clear all
close all
clc

% load predefined W matrix for 100 nodes
load mydata


% calculate combinatorial Laplacian Matrix
d = sum(W,2);
L = diag(d)-W;


% find eigenvector and eigenvalues of combinatorial Laplacian
[u v]=eig(L);


% make eignevalue as vector
v=diag(v);
% get maximum eigenvalue
lmax=max(v);


% create signal where first node is 1 rest of them zero
s=zeros(size(W,1),1);
s(1)=1;


% determine filter
flt =exp(-100*v);

% apply that filter on to graph signal
sf=u*(flt.*(u'*s));


% visualize input and result
%coord=u(:,2:4);
G=gsp_graph(W,coord);
figure;gsp_plot_signal(G,s)
title('Input signal');
figure;gsp_plot_signal(G,sf)
title('Filtered signal');


