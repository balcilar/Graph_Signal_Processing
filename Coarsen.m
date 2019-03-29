function [G, parents]= Coarsen(W,level)

% Coarsen a graph multiple times using the Heavy Edge Matching (HEM).

N=size(W,1);

ss=full(sum(W))';
ss(:,2)=1:N;

sss=sortrows(ss,1);
rid=sss(:,2);
ss(:,2)=[];

degree=full(sum(W)-diag(W)')';
G{1}=W;
NN=N;

for i=1:level
    [cc,rr,vv] = find(W);
    
    nnz = length(cc);
    marked = zeros(NN, 1);
    rowstart = ones(NN, 1);
    rowlength = zeros(NN, 1);
    cluster_id = zeros(NN, 1);
    
    oldval = rr(1);
    count = 1;
    clustercount = 1;
    
    for ii= 1:nnz
        rowlength(count) = rowlength(count) + 1;
        if rr(ii) > oldval
            oldval = rr(ii);
            rowstart(count+1) = ii;
            count = count + 1;
        end
    end
    
    for ii =1:NN
        tid = rid(ii);
        if marked(tid)==0
            wmax = 0.0;
            rs = rowstart(tid);
            marked(tid) = 1;
            bestneighbor = -1;
            for jj =1:rowlength(tid)
                nid = cc(rs+jj);
                if marked(nid)==1
                    tval = 0.0;
                else
                    Wij = vv(rs+jj);
                    Wii = vv(rowstart(tid));
                    Wjj = vv(rowstart(nid));
                    di = degree(tid);
                    dj = degree(nid);
                    tval = (2.*Wij + Wii + Wjj) * 1./(di+dj+1e-9);
                end
                
                if tval > wmax
                    wmax = tval;
                    bestneighbor = nid;
                end
            end
            
            cluster_id(tid) = clustercount;
            
            if bestneighbor > -1
                cluster_id(bestneighbor) = clustercount;
                marked(bestneighbor) = 1;
            end
            
            clustercount = clustercount+1;
        end
    end
    parents{i}=cluster_id;
    
    nrr = cluster_id(rr);
    ncc = cluster_id(cc);
    nvv = vv;
    Nnew = max(cluster_id) + 1;
    W=sparse(nrr,ncc,nvv);
    G{i+1}=W;
    
    NN=size(W,1);
    ss=full(sum(W))';
    ss(:,2)=1:NN;
    sss=sortrows(ss,1);
    rid=sss(:,2);
    ss(:,2)=[];

    
    degree=full(sum(W)-diag(W)')';
        
    
    a=1;
    
end
