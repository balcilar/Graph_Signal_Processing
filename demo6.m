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


% visualize input and result
run gspbox/gsp_start
%coord=u(:,2:4);
G=gsp_graph(W,coord);
%figure;gsp_plot_signal(G,s)
%title('Input signal');
figure;gsp_plot_signal(G,sf)
title('Filtered signal on original Graph');


W(200,200)=0;
p=randperm(200);
WW=zeros(200,200);
coord(200,1:2)=0;
corr=zeros(200,2);

for i =1:200
    for j=1:200
        WW(i,j)=W(p(i),p(j));
    end
    coor(i,:)=coord(p(i),:);
end

        

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
s(find(p==1))=1;

% determine filter
flt =exp(-20*v);
% apply that filter on to graph signal
sf2=u*(flt.*(u'*s));


G=gsp_graph(WW,coor);
%figure;gsp_plot_signal(G,s)
%title('Input signal');
figure;gsp_plot_signal(G,sf2)
title('Filtered signal on Graph which includes 100 dummy nodes');


pp=[p' [1:200]'];
pp=sortrows(pp,1);

figure;plot(sf2(pp(:,2)),'b--','linewidth',2);
hold on;plot(sf,'r-');
title('fitered nodes values');
xlabel('Node id');
legend({'node values on Graph added 100 dummy node and shuffled','node values on original graph'})




