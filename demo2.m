clear all
close all
clc

% create 2d euclidien mesh-grid  as graph
n=21;
W=zeros(n*n,n*n);
for i=1:n 
    for j=1:n
        p=sub2ind([n n],i,j);
        try
            p1=sub2ind([n n],i+1,j);
            W(p,p1)=1;
            W(p1,p)=1;
        catch
        end
        
        try
            p2=sub2ind([n n],i-1,j);
            W(p,p2)=1;
            W(p2,p)=1;
        catch
        end
        
        try
            p3=sub2ind([n n],i,j+1);
            W(p,p3)=1;
            W(p3,p)=1;
        catch
        end
        
        try
            p4=sub2ind([n n],i,j-1);   
            W(p,p4)=1;
            W(p4,p)=1;
        catch
        end
    end
end
for i=1:n 
    p1=sub2ind([n n],1,i);
    p2=sub2ind([n n],n,i);
    W(p1,p2)=1;
    W(p2,p1)=1;
    
    p1=sub2ind([n n],i,1);
    p2=sub2ind([n n],i,n);
    W(p1,p2)=1;
    W(p2,p1)=1;    
end

        
coor=zeros(n*n,2);
for i=1:n
    for j=1:n
        p=sub2ind([n n],i,j);
        coor(p,:)=[j i];
    end
end

% calculate Laplacian Matrix
d = sum(W,2);
L = diag(d)-W;


% calculate basis
[u v]=eig(L);


% make eignevalue as vector
v=diag(v);




% get maximum eigenvalue
lmax=max(v);

%figure;
GG=[];
frm=0;
for i=1:size(u,2) 
    frm=frm+1;
    uu=reshape(u(:,i),[n n]);
    ff=abs(fft2(uu));
    ff=ff(1:(n+1)/2,1:(n+1)/2);
    [a b]=max(ff(:));
    [p o]=ind2sub([(n+1)/2 (n+1)/2],b);
    
    GG=[GG;i o-1 p-1];

    
    msg=[num2str(i) 'th eigenvector measured freq: u=' num2str(o-1),' v=' num2str(p-1)];
    
    h=figure(5);
    
    imagesc(uu);
    colormap(jet);
    axis image
    axis off
    xlabel('sample suppose total 1 sec');
    %ylim([-0.2 0.2]);
    %xlim([1 n])
    title(msg);
    pause(0.1);
    frame = getframe(h);
        im = frame2im(frame);
        [imind,cm] = rgb2ind(im,256);
        %Write to the GIF File
        if frm == 1        
            imwrite(imind,cm,'eigenvectors2.gif','gif', 'Loopcount',inf);
        else   % imwrite(imind,cm,'graphembed.gif','gif', 'Loopcount',inf);     
            imwrite(imind,cm,'eigenvectors2.gif','gif','WriteMode','append');
        end        
end


% generate random iamge
s=randn(n,n);

% calculate centered frequcy of given image
F=fftshift(fft2(s));

[U V]=meshgrid(-fix(size(F,2)/2):fix(size(F,2)/2),-fix(size(F,1)/2):fix(size(F,1)/2));

UV=sqrt(U.^2+V.^2);
flt=exp(-1*UV);
% filter image in frequency domain
FD=F.*flt;
% transform the filtered image in to spatial domain
ot=ifft2(ifftshift(FD));

fltg=zeros(n*n,1);
for i=1:n*n
    i1=GG(i,2);
    j1=GG(i,3);
    fltg(i,1)=flt((n+1)/2+j1,(n+1)/2+i1);
end



sf=u*(fltg.*(u'*s(:)));
figure;
subplot(2,2,1);imagesc(s);axis image;title('Input image');
subplot(2,2,2);imagesc(ot,[-0.05 0.2]);axis image;title('Classic Filtering Result');
subplot(2,2,3);imagesc(reshape(sf,[n n ]),[-0.05 0.2]);axis image;title('Graph Signal Filtering Result');
subplot(2,2,4);imagesc(reshape(sf,[n n ])-ot);axis image;title('Differences');colorbar

%% show graph and signal
run gspbox/gsp_start
coor2=u(:,2:4);
G=gsp_graph(W,coor2);

figure;
gsp_plot_signal(G,s(:));title('input signal');
figure;
gsp_plot_signal(G,sf);title('filtered signal');



   