function Ahat= DHsuPlus(y) % dimension d, num trails

d=size(y,1);
Ahat = zeros(d);

    
    phi = randn(d,1);
    phi = phi/norm(phi,2);    
    
    G1_phi = G1(phi,y);
    G2_phi = G2(phi,y);
    G3_phi = G3(phi,y);
    M1 = G1_phi-G2_phi-2*G3_phi;
    
    for i = 1:d
        psi = randn(d,1);
        psi = psi/norm(psi,2);  
        G1_psi = G1(psi,y);
        G2_psi = G2(psi,y);
        G3_psi = G3(psi,y);
        M2 = G1_psi-G2_psi-2*G3_psi;
        M = M2/M1;  
        F = Ahat(:,1:i-1);
        M = M - F*diag(diag(F'*M*F))*F';
        %M = M*(eye(d) - F*F');
        
        lambda = eig(M);

        while ~isreal(lambda)
            psi = randn(d,1);
            psi = psi/norm(psi,2);  
            G1_psi = G1(psi,y);
            G2_psi = G2(psi,y);
            G3_psi = G3(psi,y);
            M2 = G1_psi-G2_psi-2*G3_psi;
            M = M2/M1;  
            
            F = Ahat(:,1:i-1);
            M = M - F*diag(diag(F'*M*F))*F';
            lambda = eig(M);        
        end

        x = randn(d,1);
        eps = 0.1;
        while eps>0.0001
            x1 = M*x;
            x1 = x1/norm(x1,2);
            eps = norm(x-x1,2);
            x = x1;
        end
        Ahat(:,i) = x;
    end    
    
    