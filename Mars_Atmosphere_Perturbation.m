%Initial Conditions
wm = 7.0882e-5;         %angular velocity of mars (rad/s)
R = 3396.2;             %equatorial radius of mars (km)
mu = 0.042828e6;        %gravitational parameter of Mars (km^3/s^2)

I = 55*pi/180;          %inclination of satellites
CdAm = (1/50);          %(m^2/kg)
CdAm = CdAm*(1/1000)^2; %drag coefficient*area/mass (km^2/kg)
e = 0;                  %eccentricity of circular orbit
a = R + 6855;           %initial semi-major axis = altitude + R_Mars


%Calculations
% z = 0;
% I0 = besseli(0, z)
% I1 = besseli(1, z)
% I2 = besseli(2, z)
rev = 0;

Tp = 2*pi*sqrt(a^3/mu);
atmdata = csvread('exponential_atmosphere_Wertz_1978.csv'); %read in data

n = sqrt(mu/a^3);       %mean motion (s^-1)
Q = 1-cos(I)*wm/n;
rho = exponential_atmosphere_density(atmdata, a-R);
rho = rho*(1000)^3;
rho = rho/100;          %approximating mars atmosphere density via using 1%
da_rev = -2*pi*CdAm*Q^2*a^2*rho
dI_rev = -pi*Q*CdAm*wm*a*rho*sin(I)/(2*n)
% while a-R > 6755 %while the altitude is > 200 continue calculating delta_a steps
%     n = sqrt(mu/a^3);       %mean motion (s^-1)
%     Q = 1-cos(I)*we/n;
%     rho = exponential_atmosphere_density(atmdata, a-R);
%     rho = rho*(1000)^3; %converting into km^3
%     da_rev = -2*pi*CdAm*Q^2*a^2*rho;
%     a = a+da_rev;
%     rev = rev+1;
% end
% display(rev)
% display(a)


function rho = exponential_atmosphere_density(atmdata, h)

%Calculate atmospheric density via an exponential atmosphere model
%
%Inputs:
%   h - scalar or array of heights in km
%
%Outputs:
%   rho - Atmospheric densities in kg/m^3 (same dimension as input)

h = h(:); %ensure input is a column array

%calculate density
R = zeros(length(h), 1);
L = length(atmdata);
for j = 1:length(h)
    for i = 2:L
        %ensuring that eq uses values for closest approximation of h in the
        %atmosphere data table ( next lowest <= h < next heighest)
        if (h(j) < atmdata(i,1)) && (atmdata(i-1,1) <= h(j))
        H = atmdata(i-1,3); %setting values for equation from standard atm
        rho1 = atmdata(i-1,2);
        h0 = atmdata(i-1,1);
        R(j) = rho1*exp((h0-h(j))/H); %precise rho equation for standard atm
        break
        end
    end
    if h(j) >= atmdata(L,1)
        H = atmdata(L,3); %setting values for equation from standard atm
        rho1 = atmdata(L,2);
        h0 = atmdata(L,1);
        R(j) = rho1*exp((h0-h(j))/H);
    end
end
rho = R;

%ensure output is column vector
rho = rho(:);
end