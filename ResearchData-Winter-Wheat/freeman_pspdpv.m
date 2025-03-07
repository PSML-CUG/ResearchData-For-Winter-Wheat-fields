%%% Freeman-Durden three-Component Decomposition
%%% Author:   Wenxin Xue;
%%% Date:     2024.10;
%% Setup
clear;clc;

%% Read Data
  load T3_0509.mat T3;
  T11 = T3(:,:,1);
  T12 = T3(:,:,2);
  T13 = T3(:,:,3);
  T22 = T3(:,:,4);
  T23 = T3(:,:,5);
  T33 = T3(:,:,6);
  [line,col,bs] = size(T3);
 
 %% Decomposition
% Allocate memory for parameter matrices
ps = zeros(line, col);
pd = zeros(line, col);
pv = zeros(line, col);
alpha = zeros(line, col);
beta = zeros(line, col);

% Calculate total number of pixels
N = line * col;

% Parallel processing loop
parfor jj = 1:N
    % Reconstruct T3 matrix for current pixel
    T3_matrix = [T11(jj) T12(jj) T13(jj);
                conj(T12(jj)) T22(jj) T23(jj);
                conj(T13(jj)) conj(T23(jj)) T33(jj)];
    
    % Orientation Angle Compensation
    poa0 = 0.25 * atan2(real(T23(jj)), (T22(jj) - T33(jj))/2);
    flag = 8*(T22(jj) - T33(jj))*cos(4*poa0) + 16*real(T23(jj))*sin(4*poa0);
    
    % Determine final orientation angle
    if flag > 0
        poa = poa0;
    else
        poa = poa0 + pi/4;
        flag = 8*(T22(jj) - T33(jj))*cos(4*poa) + 16*real(T23(jj))*sin(4*poa);
        if flag < 0
            poa = poa0 - pi/4;
        end
    end
    
    % Rotation matrix construction
    R3 = [1 0 0;
          0 cos(2*poa) sin(2*poa);
          0 -sin(2*poa) cos(2*poa)];
    
    % Apply rotation to T3 matrix
    T3_rot = R3 * T3_matrix * R3';
    T11(jj) = T3_rot(1,1);
    T12(jj) = T3_rot(1,2);
    T13(jj) = T3_rot(1,3);
    T22(jj) = T3_rot(2,2);
    T23(jj) = T3_rot(2,3);
    T33(jj) = T3_rot(3,3);
    
    % Check for valid matrix
    if (T3_rot(1,1) + T3_rot(2,2) + T3_rot(3,3)) < eps
        ps(jj) = 0;
        pd(jj) = 0;
        pv(jj) = 0;
        alpha(jj) = 0;
        beta(jj) = 0;
        continue;
    end
    
    % Volume scattering coefficients
    Tv11 = 0.5;   % 2.0/4.0
    Tv12 = 0.0;
    Tv22 = 0.25;  % 1.0/4.0
    Tv33 = 0.25;  % 1.0/4.0
    
    % Calculate volume scattering component
    fv = CalVolumeCoe(T11(jj), T12(jj), T13(jj), T22(jj), T23(jj), T33(jj),...
                     Tv11, Tv12, 0.0, Tv22, 0.0, Tv33);
    pv(jj) = fv;
    S = T3_rot(1,1) - pv(jj)/2;
    D = T3_rot(2,2) - pv(jj)/4;
    C = T3_rot(1,2);
    
    % Determine dominant scattering mechanism
    C0 = S - D;
    
    % Surface scattering dominant case
    if C0 > 0
        alpha(jj) = 0;
        beta(jj) = conj(C/S);
        ps(jj) = S + abs(C)^2/S;
        pd(jj) = D - abs(C)^2/S;
    end
    
    % Double-bounce scattering dominant case
    if C0 < 0
        beta(jj) = 0;
        alpha(jj) = C/D;
        ps(jj) = S - abs(C)^2/D;
        pd(jj) = D + abs(C)^2/D;
    end
end

%% Save Results
freeman(:,:,1) = pd;
freeman(:,:,2) = pv;
freeman(:,:,3) = ps;
save freeman_pspdpv_20190509.mat freeman;


%% Volume scattering coefficient
function fv = CalVolumeCoe(T11,T12,~,T22,~,T33,...
               Tv11, Tv12,~,Tv22,~,Tv33)
    
    %Initial parameters of Coherency matrix
    A11 = T11;
    A12 = T12;
    A22 = T22;
    A33 = T33;
    %Inital parameters ofcoherency matrix of Volume scattering
    B11 = Tv11;
    B12 = Tv12;
    B22 = Tv22;
    B33 = Tv33;
    
    AT = A11^2*B22^2-4*A11*real(A12)*B12*B22-2*A11*A22*B11*B22+4*A11*A22*B12^2+4*imag(A12^2)*B11*B22-4*imag(A12^2)*B12^2+4*real(A12^2)*B11*B22-4*real(A12^2)*B12^2+A22^2*B11^2;
    x1=(A11*B22-2*real(A12)*B12+A22*B11+sqrt(AT))/(2*(B11*B22-B12^2));
    x2=(A11*B22-2*real(A12)*B12+A22*B11-sqrt(AT))/(2*(B11*B22-B12^2));
    x3=A33/B33;
    
    %The Volume scattering coefficient
    fv = min([real(x1),real(x2),real(x3)]);
end

