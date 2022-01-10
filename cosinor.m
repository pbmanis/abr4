function [Amp, theta, M, p_3a] = cosinor(t,y,w,alpha, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% COSINOR 	[]=cosinor(t,y,w,alpha)
%
% Description:
%   Cosinor analysis uses the least squares method to fit a sine wave to a
%   time series. Cosinor analysis is often used in the analysis
%   of biologic time series that demonstrate predictible rhythms. This
%   method can be used with an unequally spaced time series.
%
%   Follows cosinor analysis of a time series as outlined by
%   Nelson et al. "Methods for Cosinor-Rhythmometry" Chronobiologica.
%   1979. Please consult reference.
%
% Input:
%   t - time series
%   y - value of series at time t
%   w - cycle length, defined by user based on prior knowledge of time
%       series
%   alpha - type I error used for cofidence interval calculations. Usually
%       set to be 0.05 which corresponds with 95% cofidence intervals
%
% Define Variables:
%   M - Mesor, the average cycle value
%   Amp - Amplitude, half the distance between peaks of the fitted
%       waveform
%   phi - Acrophase, time point in the cycle of highest amplitude (in
%       radians)
%   RSS - Residual Sum of Squares, a measure of the deviation of the
%       cosinor fit from the original waveform
%
% Subfunctions:
%   'CIcalc.m'
%
% Example:
%   Define time series:
%       y = [102,96.8,97,92.5,95,93,99.4,99.8,105.5];
%       t = [97,130,167.5,187.5,218,247.5,285,315,337.5]/360;
%   Define cycle length and alpha:
%       w = 2*pi;
%       alpha = .05;
%   Run Code:
%       cosinor(t,y,w,alpha)

% Record of revisions:
%     Date           Programmmer        Description of change
%     =====          ===========        ======================
%     5/16/08        Casey Cox          Original Code
%     6/24/08        Casey Cox          Revisions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 4
    error('Incorrect number of inputs');
end
if nargin > 4
    plotflag = varargin{1};
else
    plotflag = 0;
end;

if length(t) < 4
    error('There must be at least four time measurements')
end
% I. Parameter Estimation

n = length(t);

% Substituition
x = cos(w.*t);
z = sin(w.*t);

% Set up and solve the normal equations simultaneously
NE = [  n        sum(x)       sum(z)     sum(y);
    sum(x)   sum(x.^2)    sum(x.*z)    sum(x.*y);
    sum(z)   sum(x.*z)    sum(z.^2)    sum(z.*y);];

RNE = rref(NE);
M = RNE(1,4); beta = RNE(2,4); gamma = RNE(3,4); % "Mesor"

%Calculate amplitude and acrophase from beta and gamma
Amp = sqrt(beta^2 + gamma^2);
theta = atan(abs(gamma/beta));

% Calculate acrophase (phi) and convert from radians to degrees
a = sign(beta);
b = sign(gamma);
phi = 0;
if (a == 1 || a == 0) && b == 1
    phi = -theta;
elseif a == -1 && (b == 1 || b == 0)
    phi = -pi + theta;
elseif (a == -1 || a == 0) && b == -1
    phi = -pi - theta;
elseif a == 1 && (b == -1 || b == 0)
    phi = -2*pi + theta;
end
% the next part requires statistics toolbox.
% II. Confidence Limtes for Single Cosinor

%Residual sum of errors
RSS = sum((y - (M + beta.*x + gamma.*z)).^2);

%Residual varience estimation
sigma = sqrt(RSS/(n-3));

%Find confidence interval for mesor
X = 1/n * sum((x - mean(x)).^2);
Z = 1/n * sum((z - mean(z)).^2);
T = 1/n * sum((x - mean(x)).*(z - mean(z)));

%Confidence interval for the mesor
CI_M = tinv(1-alpha/2,n-3)*sigma^2*sqrt(((sum(x.^2))*(sum(z.^2)) - (sum(x.*z))^2)/(n^3*(X*Z - T^2))); %#ok<NASGU>

%Find confidence intervals for the amplitude and acrophase
[CI_Amp_min, CI_Amp_max, CI_phi_min, CI_phi_max] = CIcalc(X,T,Z,beta,gamma,n,sigma,Amp,phi,alpha, plotflag); %#ok<NASGU,NASGU>

p_3a = fpdf((n*(X*beta^2 + 2*T*beta*gamma + Z*gamma^2)/(2*sigma^2)),2,n-3);

if plotflag == 0
    return;
end;


% Display results
disp('Parameters:'); disp('---------------');
fprintf(1,'Mesor = %g \nAmplitude = %g \nAcrophase = %g \n\n',M,Amp,phi);

%Plot orginal data and cosine fit
f = M + Amp*cos(w.*t+phi);

figure('name','Cosinor Analysis: Original data and fitted function');
plot(t,y); hold on;
xlabel('x-axis');
ylabel('y-axis');
plot(t,f,'r');
legend('Original', 'Cosinor');
xlim([min(t) max(t)]);


fprintf(1,'Zero Amplitude Test \n');
fprintf(1,'------------------------------------------------------\n');
fprintf(1,'Amplitude        0.95 Confidence Limits        P Value\n');
fprintf(1,'---------        ----------------------        -------\n');
fprintf(1,' %.2f               (%.2f to %.2f)             %g\n\n',Amp,CI_Amp_min,CI_Amp_max,p_3a);


function [CI_Amp_min, CI_Amp_max, CI_phi_min, CI_phi_max] = CIcalc(X,T,Z,beta,gamma,n,sigma,Amp,phi,alpha, varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Confidence Interval Calculations
%
% Description:
%   SubFuction of 'cosinor.m'. Finds individual confidence intervals for
%   Amplitude and Acrophase and plots a polar plot representation of the
%   cosinor fit.
%
%   Follows cosinor analysis of a time series as outlined by
%   Nelson et al. "Methods for Cosinor-Rhythmometry" Chronobiologica.
%   1979. Please consult reference.
%
% Parent Function:
%   'cosinor.m'
%
% Example: Run Parent Function
%   Define time series:
%       y = [102,96.8,97,92.5,95,93,99.4,99.8,105.5];
%       t = [97,130,167.5,187.5,218,247.5,285,315,337.5]/360;
%   Define cycle length and alpha:
%       w = 2*pi;
%       alpha = .05;
%   Run Code:
%       cosinor(t,y,w,alpha)
%
% Record of revisions:
%     Date           Programmmer        Description of change
%     =====          ===========        ======================
%     6/10/08        Casey Cox          Original Code
%     6/24/08        Casey Cox          Revisions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Find beta and gamma confidence region
if nargin <= 10
    plotflag = 1;
else
    plotflag = varargin{1};
end;

F_distr = finv(1-alpha,2,n-3);

A = X;
B = 2*T;
C = Z;
D = -2*X*beta - 2*T*gamma;
E = -2*T*beta - 2*Z*gamma;
F = X*beta^2 + 2*T*beta*gamma + Z*gamma^2 - (2/n)*sigma^2*F_distr;

g_max = -(2*A*E - D*B)/(4*A*C - B^2);

gamma_s = (g_max-Amp*2:Amp/1000:g_max+Amp*2);
beta_s1 = (-(B.*gamma_s + D) + sqrt((B.*gamma_s + D).^2 - 4*A*(C.*gamma_s.^2 + E.*gamma_s + F)))/(2*A);
beta_s2 = (-(B.*gamma_s + D) - sqrt((B.*gamma_s + D).^2 - 4*A*(C.*gamma_s.^2 + E.*gamma_s + F)))/(2*A);

%Isolate ellipse region
IND = find(real(beta_s1) ~= real(beta_s2));
gamma_s = gamma_s(IND); beta_s1 = beta_s1(IND); beta_s2 = beta_s2(IND);

%Determine if confidence region overlaps the pole.
if (range(gamma_s) >= max(gamma_s)) && ((range(beta_s1) >= max(beta_s1)) || (range(beta_s2) >= max(beta_s2)))
    % disp('!! Confidence region overlaps the pole. Confidence limits for Amplitude and Acrophase cannot be determined !!');disp(' ');
    
    CI_Amp_max = 0;
    CI_Amp_min = 0;
    CI_phi_max = 0;
    CI_phi_min = 0;
else
    %Confidence Intervals for Amplitude
    CI_Amp_max = max(max([sqrt(beta_s1.^2 + gamma_s.^2); sqrt(beta_s2.^2 + gamma_s.^2)],[],2));
    CI_Amp_min = min(min([sqrt(beta_s1.^2 + gamma_s.^2); sqrt(beta_s2.^2 + gamma_s.^2)],[],2));
    
    %Confidence Intervals for Acrophase
    theta = cat(2,atan(abs(gamma_s./beta_s1)), atan(abs(gamma_s./beta_s2)));
    a = sign(cat(2,beta_s1,beta_s2));
    b = sign(cat(2,gamma_s,gamma_s))*3;
    c = a + b;
    CIphi = zeros(length(c), 1);
    for ii = 1:length(c);
        if (c(ii) == 4 || c(ii) == 3)
            CIphi(ii) = -theta(ii);
            c(ii) = 1;
        elseif (c(ii) == 2 || c(ii) == -1)
            CIphi(ii) = -pi + theta(ii);
            c(ii) = 2;
        elseif (c(ii) == -4 || c(ii) == -3)
            CIphi(ii) = -pi - theta(ii);
            c(ii) = 3;
        elseif (c(ii) == -2 || c(ii) == 1)
            CIphi(ii) = -2*pi + theta(ii);
            c(ii) = 4;
        end
    end
    if max(c) - min(c) == 3
        CI_phi_max = min(CIphi(c == 1));
        CI_phi_min = max(CIphi(c == 4));
    else
        CI_phi_max = max(CIphi);
        CI_phi_min = min(CIphi);
    end
end

if ~plotflag
    return;
end;

%Polar Representation of Cosinor analysis
figure('name','Rhythm Parameter Estimates with Joint Confidence Region', 'position', [245 357 643 600]);
plot(gamma_s,beta_s1,'linewidth', 2); hold on;
plot(gamma_s,beta_s2,'linewidth', 2)
line([0 gamma], [0 beta], 'color','k','linewidth', 2);
line([gamma -2*Amp*sin(phi)], [beta 2*Amp*cos(phi)], 'color','k','linewidth', 2, 'linestyle',':');
xlabel('\gamma');
ylabel('\beta');
ylim([-Amp*2.5 Amp*2.5]);
xlim([-Amp*2.5 Amp*2.5]);
line([0 0], [-Amp*2 Amp*2], 'color','k','linestyle', '--');
line([-Amp*2 Amp*2], [0 0], 'color','k','linestyle', '--');

%Clock and Labels
theta_clock = (0:pi/60:2*pi)';
clock = ([Amp*2*cos(theta_clock) Amp*2*sin(theta_clock)]);
plot(clock(:,1),clock(:,2), 'k', 'linewidth', 1.5);
plot(clock(:,1)*1.2,clock(:,2)*1.2, 'k', 'linewidth', 1.5);
theta_clock = (0:pi/4:2*pi-pi/4)';
clock_labels_xy = ([Amp*2.2*cos(theta_clock) Amp*2.2*sin(theta_clock)]);
clock_labels = {'6:00'; '3:00'; '0:00'; '21:00'; '18:00'; '15:00'; '12:00'; '9:00'};
for jj=1:length(clock_labels)
    text(clock_labels_xy(jj,1), clock_labels_xy(jj,2), clock_labels(jj), 'horizontalalignment', 'center');
end
