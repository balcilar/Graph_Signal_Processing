clear all
close all
clc

% create time series grid as graph
n=101;
W=zeros(n,n);
ang=linspace(0,2*pi,1010);
p=randperm(length(ang));
p=sort(p(1:101));
p=[p p(1)];

for i=1:n+1   
    loc(i,:)=[cos(ang(p(i))) sin(ang(p(i)))];
end

for i=1:n-1
    d=norm(loc(i,:)-loc(i+1,:));
    W(i,i+1)=(1/d);
    W(i+1,i)=(1/d);
end
d=norm(loc(1,:)-loc(n,:));
W(1,n)=(1/d);
W(n,1)=(1/d);

W=W/max(W(:));




% calculate combinatorial Laplacian Matrix
d = sum(W,2);
L = diag(d)-W;


% calculate basis
[u v]=eig(L);
% make eignevalue as vector
v=diag(v);
% get maximum eigenvalue
lmax=max(v);



% create arbitrary signal 
s=1:n;




%% show graph and signal
run gspbox/gsp_start



G=gsp_graph(W,loc(1:end-1,:));
figure;gsp_plot_signal(G,s)
title('original coordinates ');

G=gsp_graph(W,u(:,2:3));
figure;gsp_plot_signal(G,s)
title('esimated coordinates');



du=diff(u(:,2:3));
dt=diff(loc(1:end-1,:));
figure;plot(diag(du*du'),diag(dt*dt'),'r.')
xlabel('estimated '); 
ylabel('actual ');
title('consecutive node distance');

