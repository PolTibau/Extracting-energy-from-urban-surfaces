% HiperParameters
rho = 2300; %concrete density
k = 1.34; %thermal conductivity
cc = 1000; %specific heat of concrete
D = sqrt(rho*cc*k); 
h = 500; %Heat transfer coefficient
A = 0.2; %Albedo of concrete
d = k/(rho*cc);
Qsun = 800; %maximum solar irradiation at 14pm
Q = (1-A) .* Qsun; %solar irradiation with albedo applied
M = 293; %Inicial temperature of concrete at time 0

%Polynomial simulating sun going up and down during 16 hours
%a,b,c coefficients of second degree polynomial
K = 5;
a = (K-Q)/829440000;
b = (Q-K)/14400;
c = K;

%Optimable parameters
%These parameters can be modified to study how can we extract more energy
u = 2;               % Convection speed
z = 0.1;             % depth parameter
R = 0.01;             % R parameter

%PDE discretazion
L = 100;              % Length of the domain
Nx = 1000;           % Number of spatial points
Nt = 960;           % Number of time points (change as needed)
dx = L / (Nx-1);     % Spatial step size
dt = 0.1;        % Time step size
% Initialize solution matrix
T = 293 * ones(Nx, Nt);

% Time integration
for n = 1:Nt-1
    src = Tc(z, n*60, a, b, c, D, k, M);
    %disp(src);
    for i = 2:Nx-1
        Txx = (T(i+1, n) - 2*T(i, n) + T(i-1, n)) / dx^2;
        Tx = (T(i+1, n) - T(i-1, n)) / (2*dx);
        T(i, n+1) = T(i, n) + dt * (d * Txx - u * Tx - (2*d*h)/(R*k) * T(i, n) + (2*d*h)/(R*k) * src);
        %disp(T(i,n+1))
    end
    % Neumann Boundary Condition at x = L
    T(Nx, n+1) = T(Nx-1, n+1);
end



% Visualization
x = linspace(0, L, Nx);
times = [0,100,150,200,250,300,350,400,450,500,550,600,650,700,750,800]; % Example time points for visualization
for t = times
    plot(x, (T(:,t+1)));
    hold on;
end

xlabel('x');
ylabel('T(x,t)');
title('Temperature distribution over time');
legend(string(times));
h1 = annotation('textbox', [0.9, 0.8, 0.1, 0.1], 'String', "u = " + u);
set(h1, 'FontSize', 9, 'EdgeColor', 'none');
h2 = annotation('textbox', [0.9, 0.75, 0.1, 0.1], 'String', "z = " + z);
set(h2, 'FontSize', 9, 'EdgeColor', 'none');
h3 = annotation('textbox', [0.9, 0.7, 0.1, 0.1], 'String', "R = " + R);
set(h3, 'FontSize', 9, 'EdgeColor', 'none');


%Concrete temperature function
function src = Tc(z,t,a,b,c,D,k,M)
    src = (a.*(D.^4).*(z.^5).*erf(D.*z/(2.*k.*sqrt(t)))/(60.*(k.^5)))-(a.*(D.^4).*(z.^5)/(60.*(k.^5)))+(a.*(D.^3).*sqrt(t).*(z.^4).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t))/(30.*sqrt(pi).*(k.^4)))+(a.*(D.^2).*t.*(z.^3).*erf(D.*z/(2.*k.*sqrt(t)))/(3.*(k.^3)))-(a.*(D.^2).*t.*(z.^3)/(3.*(k.^3)))+(3.*a.*D.*(t.^(3/2)).*(z.^2).*exp(-(D.^2).*(z.^2)/(4.*t.*(k.^2)))/(5.*sqrt(pi).*(k.^2)))+(a.*(t.^2).*z.*erf(D.*z/(2.*k.*sqrt(t)))/k)-(a.*(t.^2).*z/k)+(16.*a.*(t.^(5/2)).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t))/(15.*sqrt(pi).*D))+(b.*(D.^2).*(z.^3).*erf(D.*z/(2.*k.*sqrt(t)))/(6.*(k.^3)))-(b.*(D.^2).*(z.^3)/(6.*(k.^3)))+(b.*D.*sqrt(t).*(z.^2).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t))/(3.*sqrt(pi).*(k.^2)))+(b.*t.*z.*erf(D.*z/(2.*k.*sqrt(t)))/k)-(b.*t.*z/k)+(4.*b.*(t.^(3/2)).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t))/(3.*sqrt(pi).*D))+(c.*z.*erf(D.*z/(2.*k.*sqrt(t)))/k)-(c.*z/k)+(2.*c.*sqrt(t).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t))/(sqrt(pi).*D)) + M;
end


