function Ahat= DHsuSymm_recur(y,P) % dimension d, num trails
% second version of recursion P'MP
n = size(P,1);
d = size(P,2);
if (d==0) 
    Ahat= [];
end
if (d==1) 
    Ahat= P;
end

if(d>1)
% project signals to the subspace P
    V = sqrt(-1);
    while ~isreal(V)
        phi = randn(n,1);
        phi = phi/norm(phi,2);
        G1_phi = G1(phi,y);
        G2_phi = G2(phi,y);
        G3_phi = G3(phi,y);
        M1 = G1_phi-G2_phi-2*G3_phi;
        M1 = P'*M1*P;
        
        psi = randn(n,1);
        psi = psi/norm(psi,2);
        G1_psi = G1(psi,y);
        G2_psi = G2(psi,y);
        G3_psi = G3(psi,y);
        M2 = G1_psi-G2_psi-2*G3_psi;
        M2 = P'*M2*P;
        if(sum(sum(M2~=M2)) + sum(sum(M2== Inf)) + sum(sum(M2== -Inf))==0)
            [U,D,~] = svd(M2);
            X = U*D.^(1/2);

            M = X\M1;
            M = M/X';  
            if(sum(sum(M~=M)) + sum(sum(M== Inf)) + sum(sum(M== -Inf))==0)  
               [V,sigma] = svd(M);    
               V = X*V;
            else
                V = sqrt(-1);
            end
        else
            V = sqrt(-1);    
        end 
    end
    
    
    for i=1:d
        V(:,i) = V(:,i)/norm(V(:,i),2); 
    end
    
    sigma = diag(sigma);
    [sigma, I] = sort(sigma);
    V = V(:,I);
    gap = sigma(2:d) - sigma(1:d-1);
    [~,k] = max(gap);
    
    P1 = V(:,1:k);
    P2 = V(:,k+1:d);
    W1 = DICASymm_recur(y,P*P1);
    W2 = DICASymm_recur(y,P*P2);
    
    Ahat = [W1,W2];
end