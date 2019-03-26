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
v(v<0)=0;
% get maximum eigenvalue
lmax=max(v);



% create signal where first node is 1 rest of them zero
s=zeros(size(W,1),1);
s(1)=1;

% determine filter
K=15;
alpha=rand(K,1)-0.5;
C=ones(size(W,1),1);
for i=1:K-1
    C=[C v.^(i/5)];
end
C(isinf(C))=0;
flt=C*alpha;

figure;plot(alpha);
xlabel('coeff id');
title('applied filter coeffs  \alpha');


% apply that filter on to graph signal
sf=u*diag(flt)*u'*s;


% visualize input and result
run gspbox/gsp_start
%coord=u(:,2:4);
G=gsp_graph(W,coord);
%figure;gsp_plot_signal(G,s)
%title('Input signal');
figure;gsp_plot_signal(G,sf)
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

C=ones(size(WW,1),1);
for i=1:K-1
    C=[C v.^(i/5)];
end
C(isinf(C))=0;
flt=C*alpha;


% apply that filter on to graph signal
sf2=u*diag(flt)*u'*s;


G=gsp_graph(WW,coord2);
%figure;gsp_plot_signal(G,s)
%title('Input signal');
figure;gsp_plot_signal(G,sf2)
title('Filtered signal on second graph');


