function [at,ellt,v1t,dv1t] = Lambert_Earth_Mars(r1,r2,v1,mu,dt)

%INPUTS
%   r1 - 3x1 position vector of trajectory start
%   r2 - 3x1 position vector of trajectory end (both r1 and r2 are with
%            respect to the same central body)
%   v1 - 3x1 velocity vector at trajectory start
%   mu - scalar gravitational parameter of central body (consistent
%               distance units with r1 and r2)
%   dt - scalar duration of transfer trajectory (must be positive;
%                consistent time units with mu)
%
%OUTPUTS
%   at   - scalar semi-major axis of min dv1 transfer
%   ellt - scalar semi-parameter of min dv1 transfer
%   v1t  - 3x1    velocity (with respect to central body at r1 of min dv1 
%                 transfer orbit)
%   dv1t - scalar magnitude of v1 - v1t for min dv1 transfer orbit

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CONSTANTS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

r1m = norm(r1);
r2m = norm(r2);
c = norm(r1 - r2);               %magnitude of r_1/2        
s = (1/2)*(r1m + r2m + c);       %constant s
dnu = acos(dot(r1/r1m,r2/r2m));  %change in true anomaly between r1 & r2
%(possibly edit dnu for sign change?)

%%%%%%%%%%%%%%%%%%%%%%% POSSIBLE TRAJECTORIES FOR DT %%%%%%%%%%%%%%%%%%%%%%
%find solutions
opts = optimoptions('fsolve','Display','off');
[a1,~,eflag1] = fsolve(@E1tfun,1,opts);
[a2,~,eflag2] = fsolve(@E2tfun,1,opts);
[a3,~,eflag3] = fsolve(@E3tfun,1,opts);
[a4,~,eflag4] = fsolve(@E4tfun,1,opts);
[a5,~,eflag5] = fsolve(@H1tfun,1,opts);
[a6,~,eflag6] = fsolve(@H2tfun,1,opts);
%we only care about solutions where eflag > 0

%%%%%%%%%%%%%%%%%%%%%% SEMIPARAMETERS/SEMI-MAJOR AXES %%%%%%%%%%%%%%%%%%%%%
k = 1;                  %counter for trajectories where eflag > 0
Ln = zeros(1000,1);     %semiparameter matrix (usually will only find 1-10 flagged trajectories)
An = zeros(1000,1);     %semi-major axis matrix
for i = 1:length(eflag1)        %possible trajectories for E1 orbits
   if eflag1(i) > 0             %see exit flag documentation for why must be positive
       al = 2*asin(sqrt(s/(2*a1(i))));         %constant angle
       bt = 2*asin(sqrt((s-c)/(2*a1(i))));     %constant angle
       
       Ln(k) = 4*a1(i)*(s-r1m)*(s-r2m)*(sin((al+bt)/2))^2/(c^2);
       An(k) = a1(i);
       k = k+1;
   end
end
for i = 1:length(eflag2)        %possible trajectories for E2 orbits
   if eflag2(i) > 0 
       al = 2*asin(sqrt(s/(2*a2(i))));         %constant angle
       bt = 2*asin(sqrt((s-c)/(2*a2(i))));     %constant angle
       
       Ln(k) = 4*a2(i)*(s-r1m)*(s-r2m)*(sin((al-bt)/2))^2/(c^2);
       An(k) = a2(i);
       k = k+1;
   end
end
for i = 1:length(eflag3)        %possible trajectories for E3 orbits
   if eflag3(i) > 0 
       al = 2*asin(sqrt(s/(2*a3(i))));         %constant angle
       bt = 2*asin(sqrt((s-c)/(2*a3(i))));     %constant angle
       
       Ln(k) = 4*a3(i)*(s-r1m)*(s-r2m)*(sin((al+bt)/2))^2/(c^2);
       An(k) = a3(i);
       k = k+1;
   end
end
for i = 1:length(eflag4)        %possible trajectories for E4 orbits
   if eflag4(i) > 0 
       al = 2*asin(sqrt(s/(2*a4(i))));         %constant angle
       bt = 2*asin(sqrt((s-c)/(2*a4(i))));     %constant angle
       
       Ln(k) = 4*a4(i)*(s-r1m)*(s-r2m)*(sin((al-bt)/2))^2/(c^2);
       An(k) = a4(i);
       k = k+1;
   end
end
for i = 1:length(eflag5)        %possible trajectories for H1 orbits
   if eflag5(i) > 0 
       gam = 2*asinh(sqrt(s/(-2*a5(i))));       %constant angle
       del = 2*asinh(sqrt((s-c)/(-2*a5(i))));   %constant angle
       
       Ln(k) = -4*a5(i)*(s-r1m)*(s-r2m)/c^2*(sinh((gam+del)/2))^2;
       An(k) = a5(i);
       k = k+1;
   end
end
for i = 1:length(eflag6)        %possible trajectories for H2 orbits
   if eflag6(i) > 0 
       gam = 2*asinh(sqrt(s/(-2*a6(i))));       %constant angle
       del = 2*asinh(sqrt((s-c)/(-2*a6(i))));   %constant angle
       
       Ln(k) = -4*a6(i)*(s-r1m)*(s-r2m)/c^2*(sinh((gam-del)/2))^2;
       An(k) = a6(i);
       k = k+1;
   end
end

%limit semiparameter and semimajor axis matrices to be correctly sized
%instead of size 1000x1
L = zeros(1,k);
A = zeros(1,k);
for i = 1:k
    L(i) = Ln(i);
    A(i) = An(i);
end

%%%%%%%%%%%%%%%%%%%%%%%%% DELTA V CALCULATIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%

v1s = zeros(3, k);  %vector velocity required to enter transfer orbit
dv = zeros(k,1);    %delta V array

%filling v1s with each possible trajectory's required velocity
for i = 1:k
v1s(:,i) = sqrt(mu*L(i))*(r2 - (1-r2m/L(i)*(1-cos(dnu)))*r1)/(r1m*r2m*sin(dnu));
end

%delta V requirement is the required velocity - initial velocity, for each
%possible trajectory
for i = 1:k
   dv(i) = norm(v1s(:,i) - v1);
end

%%%%%%%%%%%%%%%%%%%%%%%% MINIMUM DELTA V RESULTS %%%%%%%%%%%%%%%%%%%%%%%%%%
check = 0;
loc = 0;
dvmin = min(dv);
for i = 1:k
    if dv(i) == dvmin
        loc = i;        %location in all our matrices where dv_min occurs
        check = 1;      %minimum delta v was able to be found
        break
    end
end

at = 0;         %ensuring that Lambert Solver gives nonsense values
ellt = 0;       %if there was no possible delta Vs for the given
v1t = 100000;   %inputs to the function
dv1t = 100000;  % << possibly edit to have better solution >>

if(check == 1)
    %results
    at = A(loc);
    ellt = L(loc);
    v1t = v1s(:,loc);
    dv1t = dvmin;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%% EVENT FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    function val = E1tfun(a)
        Tp = 2*pi*sqrt(a^3/mu);             %orbital period for the given a
        al = 2*asin(sqrt(s/(2*a)));         %constant angle
        bt = 2*asin(sqrt((s-c)/(2*a)));     %constant angle
        
        %transfer time for an E1 orbit for the given a
        Te1 = Tp/(2*pi)*(al - sin(al) - (bt - sin(bt)));
        val = abs(Te1 - dt);    %if val=0 then E1 meets our dt requirement
    end

    function val = E2tfun(a)
        Tp = 2*pi*sqrt(a^3/mu);             %orbital period for given a
        al = 2*asin(sqrt(s/(2*a)));         %constant angle
        bt = 2*asin(sqrt((s-c)/(2*a)));     %constant angle
        
        %transfer time for an E2 orbit for the given a
        Te2 = Tp - Tp/(2*pi)*(al - sin(al) + (bt - sin(bt)));
        val = abs(Te2 - dt);    %if val=0 then E2 meets our dt requirement
    end

    function val = E3tfun(a)
        Tp = 2*pi*sqrt(a^3/mu);             %orbital period for given a
        al = 2*asin(sqrt(s/(2*a)));         %constant angle
        bt = 2*asin(sqrt((s-c)/(2*a)));     %constant angle
        
        %transfer time for an E3 orbit for the given a
        Te3 = Tp - Tp/(2*pi)*(al - sin(al) - (bt - sin(bt)));
        val = abs(Te3 - dt);    %if val=0 then E3 meets our dt requirement
    end

    function val = E4tfun(a)
        Tp = 2*pi*sqrt(a^3/mu);             %orbital period for given a
        al = 2*asin(sqrt(s/(2*a)));         %constant angle
        bt = 2*asin(sqrt((s-c)/(2*a)));     %constant angle
        
        %transfer time for an E4 orbit for the given a
        Te4 = Tp/(2*pi)*(al - sin(al) + (bt - sin(bt)));
        val = abs(Te4 - dt);    %if val=0 then E4 meets our dt requirement
    end

    function val = H1tfun(a)
        gam = 2*asinh(sqrt(s/(-2*a)));      %constant angle
        del = 2*asinh(sqrt((s-c)/(-2*a)));  %cosntant angle
        
        %transfer time for an E4 orbit for the given a
        Th1 = sqrt(-a^3/mu)*(sinh(gam)-gam-(sinh(del)-del));
        val = abs(Th1 - dt);    %if val=0 then E4 meets our dt requirement
    end

    function val = H2tfun(a)
        gam = 2*asinh(sqrt(s/(-2*a)));      %constant angle
        del = 2*asinh(sqrt((s-c)/(-2*a)));  %cosntant angle
        
        %transfer time for an E4 orbit for the given a
        Th2 = sqrt(-a^3/mu)*(sinh(gam)-gam+(sinh(del)-del));
        val = abs(Th2 - dt);    %if val=0 then E4 meets our dt requirement
    end
end