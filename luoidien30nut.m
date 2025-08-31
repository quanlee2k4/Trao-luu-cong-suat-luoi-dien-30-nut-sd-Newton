clear;
clc;
close all;

%% 1. Cac bien ban dau
basemva = 100;
accuracy = 5e-5;
maxiter = 50;
% Du lieu nut: so nut - Loai nut(1Slack,2PV,0PQ) - Dien ap - goc - Tai MW,
% Mvar - Phat MW, Mvar
busdata = [
 1  1  1.06   0   0.0   0.0     0.0    0.0 ;
 2  2  1.043  0  21.7  12.7    40.0    0.0 ;
 3  0  1.00   0   2.4   1.2     0.0    0.0 ;
 4  0  1.00   0   7.6   1.6     0.0    0.0 ;
 5  2  1.01   0  94.2  19.0     0.0    0.0 ;
 6  0  1.00   0   0.0   0.0     0.0    0.0 ;
 7  0  1.00   0  22.8  10.9     0.0    0.0 ;
 8  2  1.01   0  30.0  30.0     0.0    0.0 ;
 9  0  1.00   0   0.0   0.0     0.0    0.0 ;
10  0  1.00   0   5.8   2.0     0.0    0.0 ;
11  2  1.082  0   0.0   0.0     0.0    0.0 ;
12  0  1.00   0  11.2   7.5     0.0    0.0 ;
13  2  1.071  0   0.0   0.0     0.0    0.0 ;
14  0  1.00   0   6.2   1.6     0.0    0.0 ;
15  0  1.00   0   8.2   2.5     0.0    0.0 ;
16  0  1.00   0   3.5   1.8     0.0    0.0 ;
17  0  1.00   0   9.0   5.8     0.0    0.0 ;
18  0  1.00   0   3.2   0.9     0.0    0.0 ;
19  0  1.00   0   9.5   3.4     0.0    0.0 ;
20  0  1.00   0   2.2   0.7     0.0    0.0 ;
21  0  1.00   0  17.5  11.2     0.0    0.0 ;
22  0  1.00   0   0.0   0.0     0.0    0.0 ;
23  0  1.00   0   3.2   1.6     0.0    0.0 ;
24  0  1.00   0   8.7   6.7     0.0    0.0 ;
25  0  1.00   0   0.0   0.0     0.0    0.0 ;
26  0  1.00   0   3.5   2.3     0.0    0.0 ;
27  0  1.00   0   0.0   0.0     0.0    0.0 ;
28  0  1.00   0   0.0   0.0     0.0    0.0 ;
29  0  1.00   0   2.4   0.9     0.0    0.0 ;
30  0  1.00   0  10.6   1.9     0.0    0.0 ];
% Du lieu duong day: x-y : R : X : B/2 : Tap
linedata = [
 1   2  0.0192  0.0575  0.0000  1 ;
 1   3  0.0452  0.1852  0.0000  1 ;
 2   4  0.0570  0.1737  0.0000  1 ;
 3   4  0.0132  0.0379  0.0000  1 ;
 2   5  0.0472  0.1983  0.0000  1 ;
 2   6  0.0581  0.1763  0.0000  1 ;
 4   6  0.0119  0.0414  0.0000  1 ;
 5   7  0.0460  0.1160  0.0000  1 ;
 6   7  0.0267  0.0820  0.0000  1 ;
 6   8  0.0120  0.0420  0.0000  1 ;
 6   9  0.0000  0.2080  0.0000  1 ;
 6  10  0.0000  0.5560  0.0000  1 ;
 9  11  0.0000  0.2080  0.0000  1 ;
 9  10  0.0000  0.1100  0.0000  1 ;
 4  12  0.0000  0.2560  0.0000  1 ;
12  13  0.0000  0.1400  0.0000  1 ;
12  14  0.1231  0.2559  0.0000  1 ;
12  15  0.0662  0.1304  0.0000  1 ;
12  16  0.0945  0.1987  0.0000  1 ;
14  15  0.2210  0.1997  0.0000  1 ;
16  17  0.0824  0.1923  0.0000  1 ;
15  18  0.1073  0.2185  0.0000  1 ;
18  19  0.0639  0.1292  0.0000  1 ;
19  20  0.0340  0.0680  0.0000  1 ;
10  20  0.0936  0.2090  0.0000  1 ;
10  17  0.0324  0.0845  0.0000  1 ;
10  21  0.0348  0.0749  0.0000  1 ;
10  22  0.0727  0.1499  0.0000  1 ;
21  22  0.0116  0.0236  0.0000  1 ;
15  23  0.1000  0.2020  0.0000  1 ;
22  24  0.1150  0.1790  0.0000  1 ;
23  24  0.1320  0.2700  0.0000  1 ;
24  25  0.1885  0.3292  0.0000  1 ;
25  26  0.2544  0.3800  0.0000  1 ;
25  27  0.1093  0.2087  0.0000  1 ;
28  27  0.0000  0.3960  0.0000  1 ;
27  29  0.2198  0.4153  0.0000  1 ;
27  30  0.3202  0.6027  0.0000  1 ;
29  30  0.2399  0.4533  0.0000  1 ;
 8  28  0.0636  0.2000  0.0000  1 ;
 6  28  0.0169  0.0599  0.0000  1 ];


%% 2. Tao ma tran Ybus
nbus = max(max(linedata(:,1)), max(linedata(:,2)));
Ybus = zeros(nbus, nbus);
for i = 1:size(linedata, 1)
    from = linedata(i, 1);
    to = linedata(i, 2);
    R = linedata(i, 3);
    X = linedata(i, 4);
    y = 1 / (R + 1j * X);

    Ybus(from, from) = Ybus(from, from) + y;
    Ybus(to, to) = Ybus(to, to) + y;
    Ybus(from, to) = Ybus(from, to) - y;
    Ybus(to, from) = Ybus(to, from) - y;
end
disp('Ma tran Ybus:');
disp(Ybus);
G = real(Ybus);
B = imag(Ybus);

%% 3. Khoi tao
bus_code = busdata(:, 2);
V_mag = busdata(:,3);
V_ang = busdata(:,4) * pi/180;
P_sch = (busdata(:,7) - busdata(:,5)) / basemva;
Q_sch = (busdata(:,8) - busdata(:,6)) / basemva;
non_slack_buses = find(bus_code == 0 | bus_code == 2); %là PV+PQ
pq_buses = find(bus_code == 0);
V = V_mag .* exp(1j * V_ang);

%% 4. Vong lap NEWTON-RAPHSON
iter = 0;
converged = false;
fprintf('\n--- BAT DAU VONG LAP NEWTON-RAPHSON ---\n');
fprintf('Lap\tPhan du toi da\n');
fprintf('------------------------\n');

while ~converged && iter < maxiter
    iter = iter + 1;

    % a. Tinh toan cong suat va phan du
    I_inj = Ybus * V;
    S_inj = V .* conj(I_inj);
    P_calc = real(S_inj);
    Q_calc = imag(S_inj);
    delta_P = P_sch - P_calc;
    delta_Q = Q_sch - Q_calc;
    % Vector delta cho cac nut PQ va PV 
    mismatch = [delta_P(non_slack_buses); delta_Q(pq_buses)];
    % b. Kiem tra hoi tu
    max_mismatch = max(abs(mismatch));
    fprintf('%d\t\t%.6f\n', iter, max_mismatch);
    if max_mismatch < accuracy
        converged = true;
        continue;
    end

    % c. Xay dung ma tran Jacobian
    % Phan tu cua J1 (dP/d_ang)
    J1 = zeros(length(non_slack_buses), length(non_slack_buses));
    for i = 1:length(non_slack_buses)
        m = non_slack_buses(i);
        for k = 1:length(non_slack_buses)
            n = non_slack_buses(k);
            if n == m  % Phan tu duong cheo
                J1(i,k) = -Q_calc(m) - B(m,m)*(V_mag(m)^2);
            else       % Phan tu ngoai duong cheo
                J1(i,k) = V_mag(m)*V_mag(n)*(G(m,n)*sin(V_ang(m)-V_ang(n)) - B(m,n)*cos(V_ang(m)-V_ang(n)));
            end
        end
    end
    
    % Phan tu cua J2 (dP/d_mag)
    J2 = zeros(length(non_slack_buses), length(pq_buses));
     for i = 1:length(non_slack_buses)
        m = non_slack_buses(i);
        for k = 1:length(pq_buses)
            n = pq_buses(k);
            if n == m   % Phan tu duong cheo
                J2(i,k) = P_calc(m)/V_mag(m) + G(m,m)*V_mag(m);
            else        % Phan tu ngoai duong cheo
                J2(i,k) = V_mag(m)*(G(m,n)*cos(V_ang(m)-V_ang(n)) + B(m,n)*sin(V_ang(m)-V_ang(n)));
            end
        end
    end

    % Phan tu cua J3 (dQ/d_ang)
    J3 = zeros(length(pq_buses), length(non_slack_buses));
    for i = 1:length(pq_buses)
        m = pq_buses(i);
        for k = 1:length(non_slack_buses)
            n = non_slack_buses(k);
             if n == m   % Phan tu duong cheo
                J3(i,k) = P_calc(m) - G(m,m)*(V_mag(m)^2);
             else        % Phan tu ngoai duong cheo
                J3(i,k) = -V_mag(m)*V_mag(n)*(G(m,n)*cos(V_ang(m)-V_ang(n)) + B(m,n)*sin(V_ang(m)-V_ang(n)));
            end
        end
    end
    
    % Phan tu cua J4 (dQ/d_mag)
    J4 = zeros(length(pq_buses), length(pq_buses));
    for i = 1:length(pq_buses)
        m = pq_buses(i);
        for k = 1:length(pq_buses)
            n = pq_buses(k);
             if n == m   % Phan tu duong cheo
                J4(i,k) = Q_calc(m)/V_mag(m) - B(m,m)*V_mag(m);
             else        % Phan tu ngoai duong cheo
                J4(i,k) = V_mag(m)*(G(m,n)*sin(V_ang(m)-V_ang(n)) - B(m,n)*cos(V_ang(m)-V_ang(n)));
            end
        end
    end
    J = [J1, J2; J3, J4];
    vectonghiemhpt = J \ mismatch;

    % e. Cap nhat dien ap
    num_non_slack = length(non_slack_buses);
    delta_ang = vectonghiemhpt(1:num_non_slack);
    delta_mag = vectonghiemhpt(num_non_slack+1:end);
    
    V_ang(non_slack_buses) = V_ang(non_slack_buses) + delta_ang;
    V_mag(pq_buses) = V_mag(pq_buses) + delta_mag;
    
    V = V_mag .* exp(1j * V_ang);
end

%% 5. HIEN THI KET QUA
if converged
    fprintf('\nNewton-Raphson hoi tu sau %d lan lap.\n', iter);
else
    fprintf('\nNewton-Raphson khong hoi tu sau %d lan lap.\n', maxiter);
end

% Tinh toan cong suat nut bu (slack bus)
I_inj = Ybus * V;
S_inj = V .* conj(I_inj);
V_ang_deg = V_ang * 180 / pi;

fprintf('\n================================================================');
fprintf('\n                         KET QUA TRAO LUU CONG SUAT              ');
fprintf('\n================================================================');
fprintf('\nBus |  Voltage   |  Angle   | ---Load---- | ---Gen---- | ---Injected--');
fprintf('\nNo. |  Mag(pu)   |  (deg)   | MW    Mvar  | MW    Mvar  | MW      Mvar');
fprintf('\n----------------------------------------------------------------');
for i = 1:nbus
    if bus_code(i) == 1 % Nut bu
        Gen_MW = real(S_inj(i))*basemva;
        Gen_Mvar = imag(S_inj(i))*basemva;
    elseif bus_code(i) == 2 % Nut PV
        Gen_MW = busdata(i,7);
        Gen_Mvar = imag(S_inj(i))*basemva + busdata(i,6);
    else % Nut PQ
        Gen_MW = 0;
        Gen_Mvar = 0;
    end
    
    fprintf('\n %d  |  %.4f    | %8.3f | %6.1f %6.1f | %6.1f %6.1f | %8.2f %8.2f', ...
        i, V_mag(i), V_ang_deg(i), busdata(i,5), busdata(i,6), Gen_MW, Gen_Mvar, real(S_inj(i))*basemva, imag(S_inj(i))*basemva);
end
fprintf('\n\n');

% Tinh toan dong chay tren duong day va ton that
fprintf('\n================================================================');
fprintf('\n                         DONG CONG SUAT VA TON THAT              ');
fprintf('\n================================================================');
fprintf('\n Line  | ---Power at bus & line flow--- | ---Line loss--- |');
fprintf('\nfrom to|       MW              Mvar          |    MW      Mvar  |');
fprintf('\n----------------------------------------------------------------');
total_loss_P = 0;
total_loss_Q = 0;
for i = 1:size(linedata, 1)
    from = linedata(i,1);
    to = linedata(i,2);
    R = linedata(i,3);
    X = linedata(i,4);
    y_line = 1/(R + 1j*X);
    I_from_to = y_line * (V(from) - V(to));
    I_to_from = -I_from_to;
    S_from_to = V(from) * conj(I_from_to) * basemva;
    S_to_from = V(to) * conj(I_to_from) * basemva;
    
    line_loss = S_from_to + S_to_from;
    total_loss_P = total_loss_P + real(line_loss);
    total_loss_Q = total_loss_Q + imag(line_loss);
    
    fprintf('\n %2d - %2d | %12.4f %12.4f | %8.4f %8.4f', from, to, real(S_from_to), imag(S_from_to), real(line_loss), imag(line_loss));
    fprintf('\n %2d - %2d | %12.4f %12.4f |', to, from, real(S_to_from), imag(S_to_from));
    fprintf('\n----------------------------------------------------------------');
end
fprintf('\n Tong ton that                                   | %8.4f %8.4f', total_loss_P, total_loss_Q);

fprintf('\n================================================================\n');
