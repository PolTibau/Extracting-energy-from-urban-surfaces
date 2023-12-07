%Our solar irradiation starts from 0 at 6am and finishes to 0 at 22 pm. 
%Maximum solar irradiation is got at 14 pm with Qsun
%Concrete Albedo = 0.2-0.45
%Asphalt Albedo = 0.05 - 0.2
K = 5

L=293
A = 0.05
Qsun = 800
Q = (1-A) .* Qsun;
p = 2300;
cc = 1000;
k = 1.34;
D = sqrt(k .* cc .* p)
a = (K-Q)/829440000
b = (Q-K)/14400
c = K

% Define los valores de t que deseas graficar
t_values = [0, 5760, 11520, 17280, 23040, 28800, 34560, 40320, 46080, 51840, 57600];

% Define un rango de valores para z
z = linspace(0, 0.7, 1000);  % Ampliado el rango de z

% Crea una figura para el gráfico
figure;

% Calcula y grafica T(z, t) para cada valor de t en el mismo gráfico
hold on;

for i = 1:length(t_values)
    t_fixed = t_values(i);
    T_fixed = (a.*(D.^4).*(z.^5).*erf(D.*z/(2.*k.*sqrt(t_fixed)))/(60.*(k.^5)))-(a.*(D.^4).*(z.^5)/(60.*(k.^5)))+(a.*(D.^3).*sqrt(t_fixed).*(z.^4).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t_fixed))/(30.*sqrt(pi).*(k.^4)))+(a.*(D.^2).*t_fixed.*(z.^3).*erf(D.*z/(2.*k.*sqrt(t_fixed)))/(3.*(k.^3)))-(a.*(D.^2).*t_fixed.*(z.^3)/(3.*(k.^3)))+(3.*a.*D.*(t_fixed.^(3/2)).*(z.^2).*exp(-(D.^2).*(z.^2)/(4.*t_fixed.*(k.^2)))/(5.*sqrt(pi).*(k.^2)))+(a.*(t_fixed.^2).*z.*erf(D.*z/(2.*k.*sqrt(t_fixed)))/k)-(a.*(t_fixed.^2).*z/k)+(16.*a.*(t_fixed.^(5/2)).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t_fixed))/(15.*sqrt(pi).*D))+(b.*(D.^2).*(z.^3).*erf(D.*z/(2.*k.*sqrt(t_fixed)))/(6.*(k.^3)))-(b.*(D.^2).*(z.^3)/(6.*(k.^3)))+(b.*D.*sqrt(t_fixed).*(z.^2).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t_fixed))/(3.*sqrt(pi).*(k.^2)))+(b.*t_fixed.*z.*erf(D.*z/(2.*k.*sqrt(t_fixed)))/k)-(b.*t_fixed.*z/k)+(4.*b.*(t_fixed.^(3/2)).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t_fixed))/(3.*sqrt(pi).*D))+(c.*z.*erf(D.*z/(2.*k.*sqrt(t_fixed)))/k)-(c.*z/k)+(2.*c.*sqrt(t_fixed).*exp(-(D.^2).*(z.^2)/(4.*(k.^2).*t_fixed))/(sqrt(pi).*D))+L;
    plot(z, T_fixed, 'DisplayName', ['t = ', num2str(t_fixed)]);
end


xlabel('z');
ylabel('T(z, t)');
title('Graphic of T(z, t) for some values of t');
legend;
hold off;