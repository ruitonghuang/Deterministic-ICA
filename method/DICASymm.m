function Ahat= DICASymm(y) % dimension d, num trails

d=size(y,1);
%V = sqrt(-1);

%while ~isreal(V)
    phi = randn(d,1);
    G1_phi = G1(phi,y);
    G2_phi = G2(phi,y);
    G3_phi = G3(phi,y);
    M = G1_phi-G2_phi-2*G3_phi;
    [U,D,~] = svd(M);
    X = U*D.^(1/2);

    psi1 = randn(d,1);
    G1_psi1 = G1(X'\psi1,y);
    G2_psi1 = G2(X'\psi1,y);
    G3_psi1 = G3(X'\psi1,y);
    M1 = G1_psi1-G2_psi1-2*G3_psi1;
    
    psi2 = randn(d,1);
    G1_psi2 = G1(X'\psi2,y);
    G2_psi2 = G2(X'\psi2,y);
    G3_psi2 = G3(X'\psi2,y);
    M2 = G1_psi2-G2_psi2-2*G3_psi2;
    [U,D,~] = svd(M2);
    X = U*D.^(1/2);

    M3 = X\M1;
    M3 = M3/X';
    [V,~] = eig(M3);
%end
V = X*V;
if isreal(V)
    Ahat = V;
else
    disp({'DICASymm',10000});
    Ahat = eye(d);
end
  