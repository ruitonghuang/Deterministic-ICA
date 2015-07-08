function [V,y_white] = whitening(y)

d=size(y,1);
V = sqrt(-1);
while ~isreal(V)
    phi = randn(d,1);
    G1_phi = G1(phi,y);
        G2_phi = G2(phi,y);
        G3_phi = G3(phi,y);
        M1 = G1_phi-G2_phi-2*G3_phi;
        %[U,D,~] = svd(M1);
        %X = U*D.^(1/2);

        psi = randn(d,1);
        G1_psi = G1(psi,y);
        G2_psi = G2(psi,y);
        G3_psi = G3(psi,y);
        M2 = G1_psi-G2_psi-2*G3_psi;

        M = M2/M1;
        %M = M/X';    
        [V,~] = eig(M);
end
y_white  = V\y;