function [x_est,b_est,Ax_f_instance,Ay_f_instance,Az_f_instance,p_f_instance,q_f_instance,r_f_instance, AoA_f_instance] = SID230118278(c_k, d_k, t, dt)

z_k = d_k;
u_k = c_k;
Ts = dt; %time step 
N = length(t); %total number of steps

stdw = [0.01 0.01 0.01 deg2rad(0.01) deg2rad(0.01) deg2rad(0.01) 1 1 1 deg2rad(1) deg2rad(1) deg2rad(1)]; %standard deviation of w
stdv = [15 15 20 1 1 1 deg2rad(0.9) deg2rad(0.9) deg2rad(0.9) 1 deg2rad(0.9) deg2rad(0.9)];
Ex_0=[z_k(1,1) z_k(1,2) z_k(1,3) 90 0 0 z_k(1,7) z_k(1,8) z_k(1,9) 15 5 15 5 15 105 0 0 0];      %expected value of x_0
stdx_0=[5 5 5 25 25 25 5 5 5 150 150 150 150 150 150 50 50 50];  %standard deviation of x_0


xhat_km1_km1 = Ex_0; % x(0|0) = E{x_0}
P_km1_km1 = diag(stdx_0.^2); % P(0|0) = P(0)

Q = diag(stdw.^2);
R = diag(stdv.^2);

n = length(xhat_km1_km1); % n: state dimension
m = size(u_k, 2); % m: observation dimension
p = size(z_k, 2); % m: observation dimension
u_km1 = [zeros(1,m); u_k]; % shifted to have the right indices

% Preallocate storage
stdx_cor = zeros(N, n); % \sigma(k-1|k-1), standard deviation of state estimation error (hint: diagonal elements of P(k-1|k-1))
x_cor = zeros(N, n); % \hat{x}(k-1|k-1), previous estimation
K = cell(N, 1); % K(k) Kalman Gain
innov = zeros(N, p); % y(k)-y(k|k-1), innovation, with y(k|k-1)=h(\hat{x}(k|k-1),u(k|k-1),k);


% Assuming initialization and pre-allocation have been done

for k=1:N
    % Step 1: Prediction
    [t_nonlin, x_nonlin] = ode45(@(t,x) funcf(x, u_km1(k,:), t), [0 Ts], xhat_km1_km1);
    xhat_k_km1 = x_nonlin(end,:); % x(k|k-1) (prediction)

    % Step 2: Covariance matrix of state prediction error / Minimum
    % prediction MSE
    [Phi_km1, Gamma_km1] = funcLinDisDyn(xhat_km1_km1, u_km1(k,:), Ts); % Phi(k,k-1), Gamma(k,k-1)
    P_k_km1 = Phi_km1 * P_km1_km1 * Phi_km1' + Gamma_km1 * Q * Gamma_km1'; % P(k|k-1) (prediction)

    % Step 3: Kalman Gain
    H_k = funcLinDisObs(xhat_k_km1, u_km1(k,:), []);
    Ve = (H_k * P_k_km1 * H_k' + R); % Pz(k|k-1) (prediction)
    K_k = P_k_km1 * H_k' / Ve; % K(k) (gain)
    
    % Step 4: Measurement Update (Correction)
    z_k_km1 = funch(xhat_k_km1,u_km1(k,:),[]); % z(k|k-1) (prediction of output)
    innov(k,:)= z_k(k,:) - z_k_km1;
    sk(k,:)=innov(k,:)/(Ve).^(0.5);
    if mean(sk(k:11))>0.15
       innov(k,11)=0; 
       R(11,11)=deg2rad(0.01*1000)^2;
    end 
    xhat_k_k = xhat_k_km1 + ( innov(k,:))*K_k'; % x(k|k) (correction)
     
    % Step 5: Correction for Covariance matrix of state Estimate error /
    % Minimum MSE
    I_KH = eye(n) - K_k * H_k;
    P_k_k = I_KH * P_k_km1 * I_KH' + K_k * R * K_k'; % P(k|k) (correction)

    % Save data: State estimate and std dev
    stdx_cor(k,:) = sqrt(diag(P_km1_km1)); % \sigma(k-1|k-1) Standard deviation of state estimation error
    x_cor(k,:) = xhat_km1_km1; % \hat{x}(k-1|k-1), estimated state
    K{k,1} = K_k; % K(k) (gain)
    
    d_k(k,:) = z_k_km1;

    % Recursive step
    xhat_km1_km1 = xhat_k_k; 
    P_km1_km1 = P_k_k;
end




% Assuming the Kalman filter section has been executed and you have 'innov' available

% Extract relevant portions from the Kalman filter output for CUSUM analysis
x_est = [x_cor(:,1:9) x_cor(:, 16:18)];
b_est = x_cor(:, 10:15);

% Display baseline estimation plot
figure;
plot(b_est(:, 6));  % Example plot of one set of baseline estimates
title('Baseline Estimates for Sensor 6');
xlabel('Time Step');
ylabel('Estimate Value');

%% Initialize CUSUM variables
ops = cell(1, 6);  % To store indices of alarms
modified_leakage = [0.9 0.9 0.9 0.1 0.2 0.5];  % Adjusted leakage factors

% Loop through each measurement for CUSUM testing
for num = 1:6
    straingauge = b_est(500:end, num);
    baseline = 0;
    sigma = 1;

    % Dynamically calculate thresholds based on real-time variability
    threshold = 2 * sigma;  % Dynamically set threshold based on real-time std deviation
    threshold_pos = baseline + threshold;
    threshold_neg = baseline - threshold;

    g_pos = zeros(size(straingauge));
    g_neg = zeros(size(straingauge));

    % Perform CUSUM test
    for k = 2:length(straingauge)
        deviation = straingauge(k) - baseline;
        g_pos(k) = max(0, g_pos(k-1) + deviation - modified_leakage(num));
        g_neg(k) = min(0, g_neg(k-1) + deviation + modified_leakage(num));

        % Check thresholds and reset
        if g_pos(k) > threshold_pos
            ops{num} = [ops{num}, k];  % Log the index at which the threshold was exceeded
            g_pos(k) = 0;  % Reset
        end
        if g_neg(k) < threshold_neg
            ops{num} = [ops{num}, k];
            g_neg(k) = 0;  % Reset
        end
    end

    % Plotting the CUSUM results dynamically for each sensor
    figure;
    hold on;
    plot(g_pos, 'b-', 'LineWidth', 2);
    plot(g_neg, 'r-', 'LineWidth', 2);
    title(sprintf('CUSUM Test Results for Sensor %d', num));
    xlabel('Time Step');
    ylabel('CUSUM Value');
    legend('Positive CUSUM', 'Negative CUSUM');
    yline(threshold_pos, 'b--', 'Pos Threshold');
    yline(threshold_neg, 'r--', 'Neg Threshold');
    hold off;
end

%% Assign fault instance outputs based on CUSUM results
Ax_f_instance = ops{1};
Ay_f_instance = ops{2};
Az_f_instance = ops{3};
p_f_instance = ops{4};
q_f_instance = ops{5};
r_f_instance = ops{6};

%% Angle of Attack (AoA) CUSUM Test

% Extract Angle of Attack data from innovations
AoA_data = innov(:, 11); % Assuming column 11 is AoA in 'innov' matrix

% Initialize CUSUM variables
g_pos_AoA = zeros(size(AoA_data));
g_neg_AoA = zeros(size(AoA_data));
AOA_op = [];

% Determine baselines and thresholds
theta_0 = mean(AoA_data(1:100));  % Use early data for baseline calculation
sigma_0 = std(AoA_data(1:100));   % Use early data for standard deviation calculation

% Define CUSUM thresholds
threshold_pos = theta_0 + 1.5 * sigma_0;  % Set positive threshold
threshold_neg = -theta_0 - 1.5 * sigma_0; % Set negative threshold

% Perform CUSUM calculation starting from an appropriate point
for k = 101:length(AoA_data)  % Start after 100 samples to ensure data stability
    s_k = (AoA_data(k) - theta_0) / sigma_0;
    g_pos_AoA(k) = max(0, g_pos_AoA(k-1) + s_k - 0.5);  % Leakage adjusted
    g_neg_AoA(k) = min(0, g_neg_AoA(k-1) + s_k + 0.5);  % Leakage adjusted

    % Check thresholds and log alarms
    if g_pos_AoA(k) > threshold_pos
        AOA_op = [AOA_op, k];  % Log the index
        g_pos_AoA(k) = 0;       % Reset after detection
    end
    if g_neg_AoA(k) < threshold_neg
        AOA_op = [AOA_op, k];  % Log the index
        g_neg_AoA(k) = 0;       % Reset after detection
    end
end

% Assign Angle of Attack fault instances based on detections
AoA_f_instance = AOA_op(1);

% Plotting the results for Angle of Attack CUSUM
figure;
hold on;
plot(g_pos_AoA, 'b-', 'LineWidth', 2);
plot(g_neg_AoA, 'r-', 'LineWidth', 2);
title('CUSUM Analysis for Angle of Attack');
xlabel('Sample Index');
ylabel('CUSUM Value');
legend('Positive CUSUM', 'Negative CUSUM');
yline(threshold_pos, 'b--', 'Threshold Positive');
yline(threshold_neg, 'r--', 'Threshold Negative');
hold off;

   %% Functions
function x_dot = funcf(x1,u1,t1)
    %for comparison of bias and non bias and tracked bias 
    x=x1(1);y=x1(2); z=x1(3); u=x1(4); v=x1(5); w=x1(6);
    phi=x1(7); theta=x1(8); psii=x1(9); 
    b_A_x=x1(10); b_A_y=x1(11); b_A_z=x1(12); b_p=x1(13);
    b_q=x1(14); b_r=x1(15);
    v_wxE=x1(16); v_wyE=x1(17); v_wzE=x1(18);
    g=9.81;
    A_x=u1(1); A_y=u1(2); A_z=u1(3); p=u1(4);  q=u1(5);  r=u1(6); 
   
    %Paste from Lab B
    x_dott =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psii) - (v*cos(phi)-w*sin(phi))*sin(psii)+v_wxE;
    y_dot =(u*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psii) + (v*cos(phi)-w*sin(phi))*cos(psii)+v_wyE;
    z_dot =-u*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
    u_dot =(A_x-b_A_x)-g*sin(theta)+(r-b_r)*v-(q-b_q)*w;
    v_dot =(A_y-b_A_y)+g*cos(theta)*sin(phi)+(p-b_p)*w-(r-b_r)*u;
    w_dot =(A_z-b_A_z)+g*cos(theta)*cos(phi)+(q-b_q)*u-(p-b_p)*v;
    phi_dot =(p-b_p)+(q-b_q)*sin(phi)*tan(theta)+(r-b_r)*cos(phi)*tan(theta);
    theta_dot =(q-b_q)*cos(phi)-(r-b_r)*sin(phi);
    psii_dot =(q-b_q)*sin(phi)/cos(theta)+(r-b_r)*cos(phi)/cos(theta);
    b_A_x_dot =0;
    b_A_y_dot =0;
    b_A_z_dot =0;
    b_p_dot =0;
    b_q_dot =0;
    b_r_dot =0;
    v_wxE_dot =0;
    v_wyE_dot =0;
    v_wzE_dot =0;
    
    % Compute derivatives (note only here you need ; as divider between variables)
    x_dot =[x_dott; y_dot; z_dot ;u_dot; v_dot;w_dot; phi_dot; theta_dot; psii_dot; b_A_x_dot; b_A_y_dot; b_A_z_dot; ...
        b_p_dot; b_q_dot; b_r_dot; v_wxE_dot; v_wyE_dot; v_wzE_dot] ;
end

function y = funch(x1,u1,t1)
      g=9.81;
    %Paste from Lab B
    x=x1(1); y=x1(2); z=x1(3); u1=x1(4); v=x1(5); w=x1(6);
    phi=x1(7); theta=x1(8); psii=x1(9);
    
    
    v_wxE=x1(16); v_wyE=x1(17); v_wzE=x1(18);
    g=9.81;
        x_GPS=x;
        y_GPS=y;
        z_GPS=z;
        u_GPS=(u1*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*cos(psii)-(v*cos(phi)-w*sin(phi))*sin(psii)+v_wxE;
        v_GPS=(u1*cos(theta)+(v*sin(phi)+w*cos(phi))*sin(theta))*sin(psii)+(v*cos(phi)-w*sin(phi))*cos(psii)+v_wyE;
        w_GPS=-u1*sin(theta)+(v*sin(phi)+w*cos(phi))*cos(theta)+v_wzE;
        phi_GPS=phi;
        theta_GPS=theta;
        psii_GPS=psii;
        Vtas=sqrt(u1^2+v^2+w^2);
        alpha=atan(w/u1);
        beta=atan(v/sqrt(u1^2+w^2));
    % Compute output
    y = [x_GPS y_GPS z_GPS u_GPS v_GPS w_GPS phi_GPS theta_GPS psii_GPS Vtas alpha beta];
end


function [Phi,Gamma] = funcLinDisDyn(x_linpt,u_linpt,Ts)
     %for comparison of bias and non bias and tracked bias 
    % Extract values
    g=9.81;
    x=x_linpt(1); y=x_linpt(2); z=x_linpt(3); u=x_linpt(4); v=x_linpt(5); w=x_linpt(6);
    phi=x_linpt(7); theta=x_linpt(8); psii=x_linpt(9);
    b_A_x=x_linpt(10); b_A_y=x_linpt(11); b_A_z=x_linpt(12);
    b_p=x_linpt(13);
    b_q=x_linpt(14); b_r=x_linpt(15);
    v_wxE=x_linpt(16); v_wyE=x_linpt(17); v_wzE=x_linpt(18);
    A_xm=u_linpt(1); A_ym=u_linpt(2); A_zm=u_linpt(3); 
    pm=u_linpt(4); qm=u_linpt(5); rm=u_linpt(6);
   %[omega_A_x ,omega_A_y, omega_A_z, omega_p ,omega_q, omega_r]=deal(0);
    
    
  
    
    % Numerical evaluation of continuous - time dynamics
    % Note these matrices could be different from your Lab B ones due to
    % input measurement noises
 F =[...
0, 0, 0, cos(psii)*cos(theta), cos(psii)*sin(phi)*sin(theta) - cos(phi)*sin(psii), sin(phi)*sin(psii) + cos(phi)*cos(psii)*sin(theta), sin(psii)*(w*cos(phi) + v*sin(phi)) + cos(psii)*sin(theta)*(v*cos(phi) - w*sin(phi)),                                -cos(psii)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psii)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psii)*(v*cos(phi) - w*sin(phi)),  0,  0,  0,  0,                    0,                    0, 1, 0, 0;
0, 0, 0, cos(theta)*sin(psii), cos(phi)*cos(psii) + sin(phi)*sin(psii)*sin(theta), cos(phi)*sin(psii)*sin(theta) - cos(psii)*sin(phi), sin(psii)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psii)*(w*cos(phi) + v*sin(phi)),                                -sin(psii)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psii)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psii)*(v*cos(phi) - w*sin(phi)),  0,  0,  0,  0,                    0,                    0, 0, 1, 0;
0, 0, 0,          -sin(theta),                                cos(theta)*sin(phi),                                cos(phi)*cos(theta),                                                 cos(theta)*(v*cos(phi) - w*sin(phi)),                                           - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 1;
0, 0, 0,                    0,                                           rm - b_r,                                           b_q - qm,                                                                                    0,                                                                                   -g*cos(theta),                                                                                                       0, -1,  0,  0,  0,                    w,                   -v, 0, 0, 0;
0, 0, 0,             b_r - rm,                                                  0,                                           pm - b_p,                                                                g*cos(phi)*cos(theta),                                                                          -g*sin(phi)*sin(theta),                                                                                                       0,  0, -1,  0, -w,                    0,                    u, 0, 0, 0;
0, 0, 0,             qm - b_q,                                           b_p - pm,                                                  0,                                                               -g*cos(theta)*sin(phi),                                                                          -g*cos(phi)*sin(theta),                                                                                                       0,  0,  0, -1,  v,                   -u,                    0, 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                      sin(phi)*tan(theta)*(b_r - rm) - cos(phi)*tan(theta)*(b_q - qm),               - cos(phi)*(b_r - rm)*(tan(theta)^2 + 1) - sin(phi)*(b_q - qm)*(tan(theta)^2 + 1),                                                                                                       0,  0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta), 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                            cos(phi)*(b_r - rm) + sin(phi)*(b_q - qm),                                                                                               0,                                                                                                       0,  0,  0,  0,  0,            -cos(phi),             sin(phi), 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                  (sin(phi)*(b_r - rm))/cos(theta) - (cos(phi)*(b_q - qm))/cos(theta), - (cos(phi)*sin(theta)*(b_r - rm))/cos(theta)^2 - (sin(phi)*sin(theta)*(b_q - qm))/cos(theta)^2,                                                                                                       0,  0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta), 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                                                                    0,                                                                                               0,                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                                                                    0,                                                                                               0,                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                                                                    0,                                                                                               0,                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                                                                    0,                                                                                               0,                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                                                                    0,                                                                                               0,                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                                                                    0,                                                                                               0,                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                                                                    0,                                                                                               0,                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                                                                    0,                                                                                               0,                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0;
0, 0, 0,                    0,                                                  0,                                                  0,                                                                                    0,                                                                                               0,                                                                                                       0,  0,  0,  0,  0,                    0,                    0, 0, 0, 0];
G = [0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,    0,    0;
 0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,    0,    0;
 0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,    0,    0;
-1,  0,  0,  0,                    w,                   -v,    0,    0,    0,    0,    0,    0;
 0, -1,  0, -w,                    0,                    u,    0,    0,    0,    0,    0,    0;
 0,  0, -1,  v,                   -u,                    0,    0,    0,    0,    0,    0,    0;
 0,  0,  0, -1, -sin(phi)*tan(theta), -cos(phi)*tan(theta),    0,    0,    0,    0,    0,    0;
 0,  0,  0,  0,            -cos(phi),             sin(phi),    0,    0,    0,    0,    0,    0;
 0,  0,  0,  0, -sin(phi)/cos(theta), -cos(phi)/cos(theta),    0,    0,    0,    0,    0,    0;
 0,  0,  0,  0,                    0,                    0, 1/Ts,    0,    0,    0,    0,    0;
 0,  0,  0,  0,                    0,                    0,    0, 1/Ts,    0,    0,    0,    0;
 0,  0,  0,  0,                    0,                    0,    0,    0, 1/Ts,    0,    0,    0;
 0,  0,  0,  0,                    0,                    0,    0,    0,    0, 1/Ts,    0,    0;
 0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0, 1/Ts,    0;
 0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,    0, 1/Ts;
 0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,    0,    0;
 0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,    0,    0;
 0,  0,  0,  0,                    0,                    0,    0,    0,    0,    0,    0,    0];
 
 
    
    % Discretisation of dynamics
    [Phi , Gamma ]= c2d (F ,G , Ts ) ;
end


    function H = funcLinDisObs(x_linpt,u_linpt,t)
%for comparison of bias and non bias and tracked bias 
    g=9.81;
    x=x_linpt(1); y=x_linpt(2); z=x_linpt(3); u=x_linpt(4); v=x_linpt(5); w=x_linpt(6);
    phi=x_linpt(7); theta=x_linpt(8); psii=x_linpt(9);
    b_A_x=x_linpt(10); b_A_y=x_linpt(11); b_A_z=x_linpt(12);
    b_p=x_linpt(13);
    b_q=x_linpt(14); b_r=x_linpt(15);
    v_wxE=x_linpt(16); v_wyE=x_linpt(17); v_wzE=x_linpt(18);
    A_xm=u_linpt(1); A_ym=u_linpt(2); A_zm=u_linpt(3); 
    pm=u_linpt(4); qm=u_linpt(5); rm=u_linpt(6);
    
    % Numerical evaluation
    % Note these matrices could be different from your Lab B ones due to
    % input measurement noises
    H=[...
    1, 0, 0,                                                0,                                                  0,                                                  0,                                                                                    0,                                                                0,                                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 1, 0,                                                0,                                                  0,                                                  0,                                                                                    0,                                                                0,                                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 1,                                                0,                                                  0,                                                  0,                                                                                    0,                                                                0,                                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                             cos(psii)*cos(theta), cos(psii)*sin(phi)*sin(theta) - cos(phi)*sin(psii), sin(phi)*sin(psii) + cos(phi)*cos(psii)*sin(theta), sin(psii)*(w*cos(phi) + v*sin(phi)) + cos(psii)*sin(theta)*(v*cos(phi) - w*sin(phi)), -cos(psii)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))), - sin(psii)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - cos(psii)*(v*cos(phi) - w*sin(phi)), 0, 0, 0, 0, 0, 0, 1, 0, 0;
    0, 0, 0,                             cos(theta)*sin(psii), cos(phi)*cos(psii) + sin(phi)*sin(psii)*sin(theta), cos(phi)*sin(psii)*sin(theta) - cos(psii)*sin(phi), sin(psii)*sin(theta)*(v*cos(phi) - w*sin(phi)) - cos(psii)*(w*cos(phi) + v*sin(phi)), -sin(psii)*(u*sin(theta) - cos(theta)*(w*cos(phi) + v*sin(phi))),   cos(psii)*(sin(theta)*(w*cos(phi) + v*sin(phi)) + u*cos(theta)) - sin(psii)*(v*cos(phi) - w*sin(phi)), 0, 0, 0, 0, 0, 0, 0, 1, 0;
    0, 0, 0,                                      -sin(theta),                                cos(theta)*sin(phi),                                cos(phi)*cos(theta),                                                 cos(theta)*(v*cos(phi) - w*sin(phi)),            - sin(theta)*(w*cos(phi) + v*sin(phi)) - u*cos(theta),                                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 1;
    0, 0, 0,                                                0,                                                  0,                                                  0,                                                                                    1,                                                                0,                                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                                                0,                                                  0,                                                  0,                                                                                    0,                                                                1,                                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                                                0,                                                  0,                                                  0,                                                                                    0,                                                                0,                                                                                                       1, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                        u/(u^2 + v^2 + w^2)^(1/2),                          v/(u^2 + v^2 + w^2)^(1/2),                          w/(u^2 + v^2 + w^2)^(1/2),                                                                                    0,                                                                0,                                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0,                           -w/(u^2*(w^2/u^2 + 1)),                                                  0,                                1/(u*(w^2/u^2 + 1)),                                                                                    0,                                                                0,                                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0;
    0, 0, 0, -(u*v)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),        1/((u^2 + w^2)^(1/2)*(v^2/(u^2 + w^2) + 1)),   -(v*w)/((u^2 + w^2)^(3/2)*(v^2/(u^2 + w^2) + 1)),                                                                                    0,                                                                0,                                                                                                       0, 0, 0, 0, 0, 0, 0, 0, 0, 0];

    end    


end
