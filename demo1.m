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
% s=sign(u(1,:));
% s(s==0)=1;
% u=u.*s; 

figure;subplot(1,3,1);plot(u(:,1));ylim([-0.15 0.15]);xlim([1 101]);
title(['1st eigenvector, value=' num2str(v(1,1))]);
subplot(1,3,2);plot(u(:,2));ylim([-0.15 0.15]);xlim([1 101]);
title(['2nd eigenvector, value=' num2str(v(2,2))]);
subplot(1,3,3);plot(u(:,3));ylim([-0.15 0.15]);xlim([1 101]);
title(['3th eigenvector, value=' num2str(v(3,3))]);


% make eignevalue as vector
v=diag(v);
% get maximum eigenvalue
lmax=max(v);


% create arbitrary signal 
s=randn(size(W,1),1);
%s(50)=1;

% calcualte graph frequency profile
f=u'*s;

% calculate classical frequency profile
cf=abs(fft(s));
cf=cf(1:51);

% F=[0;f];
% F=reshape(F,2,51);
% f=sum(abs(F));

figure;subplot(1,2,1);plot(v,f,'r*-');%plot(v(1:2:end),f,'r*-');
xlabel('eigenvalues');ylabel('coefficient');
title('Graph Frequency profile');

subplot(1,2,2);plot(0:50,cf,'b*-');
xlabel('Frequencies');ylabel('magnitude');
title('Frequency profile');


% create filter for classical signal processing
fs=fftshift(fft(s));
f=linspace(-(n-1)/2,(n-1)/2,n)';
flt = exp(-abs(f)*0.1);
% apply that filter on to time series signal
ot=ifft(ifftshift(flt.*fs));


% create the same filter for graph signal processing
flt=[1:(n-1)/2]';
flt=[flt flt]';
flt=[0; flt(:)];
flt = exp(-abs(flt)*0.1);

% apply that signal on to graph signal
sf=u*(flt.*(u'*s));

figure;
hold on;plot(s,'k--');plot(ot,'b-','linewidth',2);xlim([1 n]);
plot(sf,'r-','linewidth',1);xlabel('node order or time');
plot(ot-sf,'g-')
legend({'original signal','graph filtering','time series filtering','differences'});





%% show graph and signal
run gspbox/gsp_start


coor=u(:,2:3);
G=gsp_graph(W,coor);
figure;gsp_plot_signal(G,s)

title('Input signal');
figure;gsp_plot_signal(G,sf)
title('Filtered signal');




figure;
frm=0;
for i=1:size(u,2) 
    frm=frm+1;
    ff=abs(fft(u(:,i)));
    [a b]=max(ff(1:(length(ff)+1)/2));
    msg=[num2str(i) 'th eigenvector measured freq:' num2str(b-1) 'Hz'];
    
    h=figure(5);
    
    plot(u(:,i));
    xlabel('sample suppose total 1 sec');
    ylim([-0.2 0.2]);
    xlim([1 n])
    title(msg);
    %pause(0.1);
    frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        %Write to the GIF File
        if frm == 1        
            imwrite(imind,cm,'eigenvectors11.gif','gif', 'Loopcount',inf);
        else        
            imwrite(imind,cm,'eigenvectors11.gif','gif','WriteMode','append');
        end
        
end







