function Ahat= DICA_recur(y,P) % dimension d, num trails
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

    V = sqrt(-1);
    while ~isreal(V)
        phi = randn(n,1);
        G1_phi = G1(phi,y);
        G2_phi = G2(phi,y);
        G3_phi = G3(phi,y);
        M = G1_phi-G2_phi-2*G3_phi;
        if(sum(sum(M~=M)) + sum(sum(M== Inf)) + sum(sum(M== -Inf))==0)
            
            [U,D,~] = svd(M);
            X = U*D.^(1/2);

            psi1 = randn(n,1);       
            G1_psi1 = G1(X'\psi1,y);
            G2_psi1 = G2(X'\psi1,y);
            G3_psi1 = G3(X'\psi1,y);
            M1 = G1_psi1-G2_psi1-2*G3_psi1;
            %T1 = X\M1;
            %T1 = T1/X';

            psi2 = randn(n,1);
            G1_psi2 = G1(X'\psi2,y);
            G2_psi2 = G2(X'\psi2,y);
            G3_psi2 = G3(X'\psi2,y);
            M2 = G1_psi2-G2_psi2-2*G3_psi2;
            %T2 = X\M2;
            %T2 = T2/X';

            M3 = M1/M2;
            M3 = P'*M3*P;
            if(sum(sum(M3~=M3)) + sum(sum(M3 == Inf)) + sum(sum(M3 == -Inf))==0)  
                [V,sigma] = eig(M3);
            else
                V = sqrt(-1);
            end
        else
            V = sqrt(-1);
        end
    end
   
    sigma = diag(sigma);
    [sigma, I] = sort(sigma);
    V = V(:,I);
    for i=1:d
        V(:,i) = V(:,i)/norm(V(:,i),2); 
    end
    
    gap = sigma(2:d) - sigma(1:d-1);
    [~,k] = max(gap);
    
    P1 = V(:,1:k);
    P2 = V(:,k+1:d);
    W1 = DICA_recur(y,P*P1);
    W2 = DICA_recur(y,P*P2);
    
    Ahat = [W1,W2];
   

end
  