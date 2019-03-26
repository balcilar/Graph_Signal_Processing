function basis=bspline_basis(K, v, x, degree)

kv1 = min(v) * ones(degree,1);
kv2 = linspace(min(v), max(v), K-degree+1)';
kv3 = max(v) * ones(degree,1);
kv = [kv1 ;kv2 ;kv3];

for k =1:K
    basis(:,k)=cox_deboor(k, degree);
end
%basis(end,end)=1;



    function ret= cox_deboor(k, d)
        
        % Test for end conditions, the rectangular degree zero spline.
        if (d == 0)
            ret= double(((x - kv(k) >= 0) & (x - kv(k + 1) < 0)));
        else
            
            denom1 = kv(k + d) - kv(k);
            term1 = 0;
            if denom1 > 0
                term1 = ((x - kv(k)) / denom1) .* cox_deboor(k, d - 1);
            end
            
            denom2 = kv(k + d + 1) - kv(k + 1);
            term2 = 0;
            if denom2 > 0
                term2 = ((-(x - kv(k + d + 1)) / denom2) .* cox_deboor(k + 1, d - 1));
            end
            ret= term1 + term2;
        end
        
    end

end