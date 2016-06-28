if exist('Fs') == 0
    load('DonneesActionneurs\Fe_attraction')
    load('DonneesActionneurs\Fs')
end

i = [-1, -2];
be1 = 13.029359254409743;
num = (i .^ 2 + be1 .* abs(i)) .* sign(i);

fs_ajustee = Fs .^ -1;
fe_m1A_ajustee = (Fe_m1A / num(1)) .^ -1;
fe_m2A_ajustee = (Fe_m2A / num(2)) .^ -1;

A_s = [z_pos .^ 0, z_pos, z_pos .^ 2, z_pos .^ 3];
A_mA1 = [z_m1A .^ 0, z_m1A, z_m1A .^ 2, z_m1A .^ 3];
A_mA2 = [z_m2A .^ 0, z_m2A, z_m2A .^ 2, z_m2A .^ 3];

X_s = pinv(A_s) * fs_ajustee;
X_mA1 = pinv(A_mA1) * fe_m1A_ajustee;
X_mA2 = pinv(A_mA2) * fe_m2A_ajustee;

fs_calculee = -1 ./ (X_s(4) * z_pos .^ 3 + ... 
                     X_s(3) * z_pos .^ 2 + ...
                     X_s(2) * z_pos + ...
                     X_s(1));

fe_mA1_calculee = num(1) ./ (X_mA1(4) * z_m1A .^ 3 + ... 
                             X_mA1(3) * z_m1A .^ 2 + ...
                             X_mA1(2) * z_m1A + ...
                             X_mA1(1));

fe_mA2_calculee = num(2) ./ (X_mA2(4) * z_m2A .^ 3 + ... 
                             X_mA2(3) * z_m2A .^ 2 + ...
                             X_mA2(2) * z_m2A + ...
                             X_mA2(1));
                         
figure
subplot(3, 1, 1)
plot(z_pos, fs_calculee, 'b', z_pos, Fs, 'r.')
xlabel('Distance en z (m)')
ylabel('Force (N)')
title('Force naturelle de l''aimant par rapport à la distance de la plaque')

subplot(3, 1, 2)
plot(z_m1A, fe_calculee_mA1, 'b', z_m1A, Fe_m1A, 'r.')
xlabel('Distance en z (m)')
ylabel('Force (N)')
title('Force de l''électroaimant avec un courant d''un ampère par rapport à la distance de la plaque')

subplot(3, 1, 3)
plot(z_m2A, fe_calculee_mA2, 'b', z_m2A, Fe_m2A, 'r.')
xlabel('Distance en z (m)')
ylabel('Force (N)')
title('Force de l''électroaimant avec un courant de deux ampères par rapport à la distance de la plaque')






