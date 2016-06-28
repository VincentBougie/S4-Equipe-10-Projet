clear all

if exist('mS') == 0
    run('bancEssaiConstantes')
end

b_E1 = 13.029359254409743;
positionPlaqueZ =  0.015;
phi = 0;
theta = 0;
positionSphere = [1, 1]; %[x, y]
masseSphere = mS;

positionActionneurs = [xvec_ABC, yvec_ABC, zeros(3,1)];
positionCapteurs = [xvec_DEF, yvec_DEF, zeros(3,1)];

Zk_eq = positionPlaqueZ - positionActionneurs(:, 1) * theta + positionActionneurs(:, 2) * phi;

coefFE = [1, 1, 1, 1];
coefFS = [1, 1, 1, 1];

Ck = [(XB * g * masseSphere * positionSphere(2) - XC * g * masseSphere * positionSphere(2) - YB * g * masseSphere * positionSphere(1) + YC * g * masseSphere * positionSphere(1) - XB * YC * g * mP + XC * YB * g * mP - XB * YC * g * masseSphere + XC * YB * g * masseSphere) / (XA * YB - XB * YA - XA * YC + XC * YA + XB * YC - XC * YB);
      -(XA * g * masseSphere * positionSphere(2) - XC * g * masseSphere * positionSphere(2) - YA * g * masseSphere * positionSphere(1) + YC * g * masseSphere * positionSphere(1) - XA * YC * g * mP + XC * YA * g * mP - XA * YC * g * masseSphere + XC * YA * g * masseSphere) / (XA * YB - XB * YA - XA * YC + XC * YA + XB * YC - XC * YB);
      (XA * g * masseSphere * positionSphere(2) - XB * g * masseSphere * positionSphere(2) - YA * g * masseSphere * positionSphere(1) + YB * g * masseSphere * positionSphere(1) - XA * YB * g * mP + XB * YA * g * mP - XA * YB * g * masseSphere + XB * YA * g * masseSphere) / (XA * YB - XB * YA - XA * YC + XC * YA + XB * YC - XC * YB)];

I_eq = b_E1 - sqrt(b_E1 ^ 2 - 4 * sign(Ck(:, 1)) .* ((coefFS(1) * Zk_eq .^ 0 + coefFS(2) * Zk_eq + coefFS(3) * Zk_eq .^ 2 + coefFS(4) * Zk_eq .^ 3)) .* (Ck + 1 ./ (coefFS(1) * Zk_eq .^ 0 + coefFS(2) * Zk_eq + coefFS(3) * Zk_eq .^ 2 + coefFS(4) * Zk_eq .^ 3))) ./ (2 * sign(Ck));

u_eq = [I_eq(1) / RA;
        I_eq(2) / RB;
        I_eq(3) / RC];

%[phi, theta, z, omega_phi, omega_theta, vitesse_z, x_sphere, y_sphere, v_x_sphere, v_y_sphere, courant_A, courant_B, courant_C]
conditionsInitiales_eq = [phi, theta, positionPlaqueZ, 0, 0, 0, positionSphere(1), positionSphere(2), 0, 0, I_eq'];

derivePartForceCourant = (2 * abs(I_eq) + b_E1) ./ (coefFE(1) + coefFE(2) * Zk_eq + coefFE(3) * Zk_eq .^ 2 + coefFE(4) * Zk_eq .^ 3);
derivePartForceZ = -((I_eq .^ 2 + b_E1 * abs(I_eq)) .* sign(I_eq)) .* (coefFE(2) * Zk_eq .^ 0 + 2 * coefFE(3) * Zk_eq + 3 * coefFE(4) * Zk_eq .^ 2) ./ ((coefFE(1) * Zk_eq .^ 0 + coefFE(2) * Zk_eq + coefFE(3) * Zk_eq .^ 2 + coefFE(4) * Zk_eq .^ 3) .^ 2) + ...
                   (coefFS(2) * Zk_eq .^ 0 + 2 * coefFS(3) * Zk_eq + 3 * coefFS(4) * Zk_eq .^ 2) ./ ((coefFS(1) * Zk_eq .^ 0 + coefFS(2) * Zk_eq + coefFS(3) * Zk_eq .^ 2 + coefFS(4) * Zk_eq .^ 3) .^ 2);
derivePartForceTheta = -positionActionneurs(:, 1) .* derivePartForceZ;
derivePartForcePhi = positionActionneurs(:, 2) .* derivePartForceZ;

%Sous-matrices
PS = [0,(-mS*g)/Jx;
      (mS*g)/Jx, 0;
      0, 0];
SP = [0, (-5/7)*g,0;
      (5/7)*g,0,0];
CC = [-RA/LA, 0, 0;
      0, -RB/LB, 0;
      0, 0, -RC/LC];
CV = [1/LA, 0, 0;
      0, 1/LB, 0;
      0, 0, 1/LC];
PP = [-sum(positionActionneurs(2, :) * derivePartForcePhi) / Jx, -sum(positionActionneurs(2, :) * derivePartForceTheta) / Jx, -sum(positionActionneurs(2, :) * derivePartForceZ) / Jx;
       sum(positionActionneurs(1, :) * derivePartForcePhi) / Jx, sum(positionActionneurs(1, :) * derivePartForceTheta) / Jx, sum(positionActionneurs(1, :) * derivePartForceZ) / Jx;
       sum(derivePartForcePhi) / (mP + mS), sum(derivePartForceTheta) / (mP + mS), sum(derivePartForceZ) / (mP + mS)];
PC = [-(positionActionneurs(2, 1) * derivePartForceZ(1, 1)) / Jx, -(positionActionneurs(2, 2) * derivePartForceZ(2, 1)) / Jx, -(positionActionneurs(2, 3) * derivePartForceZ(3, 1)) / Jx;
       (positionActionneurs(1, 1) * derivePartForceZ(1, 1)) / Jx, (positionActionneurs(1, 2) * derivePartForceZ(2, 1)) / Jx, (positionActionneurs(1, 3) * derivePartForceZ(3, 1)) / Jx;
       derivePartForceZ(1) / (mS + mP), derivePartForceZ(2) / (mS + mP), derivePartForceZ(3) / (mS + mP)];

%Matrices de linearisation
A = zeros(13, 13);
A(1:3, 1:3) = zeros(3, 3);
A(1:3, 4:6) = eye(3, 3);
A(1:3, 7:13) = zeros(3, 7);
A(4:6, 1:3) = PP;
A(4:6, 4:6) = zeros(3, 3);
A(4:6, 7:8) = PS;
A(4:6, 9:10) = zeros(3, 2);
A(4:6, 11:13) = PC;
A(7:8, 1:8) = zeros(2, 8);
A(7:8, 9:10) = eye(2, 2);
A(7:8, 11:13) = zeros(2, 3);
A(9:10, 1:3) = SP;
A(9:10, 4:13) = zeros(2, 10);
A(11:13, 1:10) = zeros(3, 10);
A(11:13, 11:13) = CC;

B = zeros(13, 3);
B(1:10, 1:3) = zeros(10, 3);
B(11:13, 1:3) = CV;

C = zeros(7, 13);
C(1:3, 1:3) = TDEF;
C(1:3, 4:13) = zeros(3, 10);
C(4:7, 1:6) = zeros(4, 6);
C(4:7, 7:10) = eye(4, 4);
C(4:7, 11:13) = zeros(4, 3);

D = zeros(7, 3);








