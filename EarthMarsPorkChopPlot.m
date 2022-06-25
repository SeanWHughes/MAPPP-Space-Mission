
%%%%%%%%%%%%%%%%%%%%%%%%%%%%% POSSIBLE TRAJECTORIES %%%%%%%%%%%%%%%%%%%%%%%

mus = 2.9591220828559093e-4;         %gravitational parameter for the sun (central body) (AU^3*day^-2)

%reading JPL Horizons ephemeris data for Mars & Earth
marsdata = readHorizons('D:\MATLABSCRIPTS\Spaceflight Mechanics\Group Project\Mars_Horizons.txt',true);
earthdata = readHorizons('D:\MATLABSCRIPTS\Spaceflight Mechanics\Group Project\Earth_Horizons.txt',true);

%converting raw cell array data into normal arrays
Jt = cell2mat(marsdata(:,1));        %juliandate time array
r_m = cell2mat(marsdata(:,3:5));     %mars position array (AU)
v_m = cell2mat(marsdata(:,6:8));     %mars velocity array (AU/day) 
r_e = cell2mat(earthdata(:,3:5));    %earth position array (AU)
v_e = cell2mat(earthdata(:,6:8));    %earth velocity array (AU/day)
%all trajectory arrays are now (124x3) matrices

%setting all trajectory arrays to (3x124) matrices
r_m = r_m.';
v_m = v_m.';
r_e = r_e.';
v_e = v_e.';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELTA T MATRIX %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%delta t matrices setup, mission duration is required to be one month
t01 = juliandate(datetime('2035-07-01'));   %earliest possible launch date
t02 = juliandate(datetime('2035-10-01'));   %latest possible launch date 
tf1 = juliandate(datetime('2035-08-01'));   %one month after earliest possible launch date
tf2 = juliandate(datetime('2035-11-01'));   %one month after latest possible launch date

tdep = t01:1:t02;                           %all possible departure days
tarr = tf1:1:tf2;                           %all possible arrival days

[d,r] = meshgrid(tdep,tarr);    %creating array of (tdep,tarr) coordinates
%r = 93x93 array where each row is tdep(:)
%d = 93x93 array where each column is tarr(:)
dta = r-d;         %gives all possible constant mission durations (93x93)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MIN DV SEARCH %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

m = length(tarr);
n = length(tdep);
dvt = zeros(m,n);  %initializing delta v matrix as same size as dta
v0t = [0;0;0];
r0t = [0;0;0];
for i = 1:m
    for j = 1:n
        if(dta(i,j) > 10)   %requiring that mission duration input is > 10 days (will require 30 days later)
        %[r_e,r_v,v_e,v_v] = EarthVenusPosVel(tdep(j),dta(i,j));
        [~,~,v1t,dv1t] = Lambert_Earth_Mars(r_e(:,j),r_m(:,j+dta(i,j)),v_e(:,j),mus,dta(i,j));
            if(dv1t ~= 100000)                          %dv1t = 100000 if no min delta V could be found
                dvt(i,j) = dv1t*149597870.7/(86400);    %conversions from AU/day --> km/s
%                 if(dvt(i,j) > 20 && dvt(i,j) < 21)
%                     r0t = r_e(:,j)
%                     v0t = v1t
%                 end
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%% PLOTTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Calender Date axis conversion from Savransky's github:
nmonthsdep = calmonths(caldiff([datetime(tdep(1),'ConvertFrom','JD'),datetime(tdep(end),'ConvertFrom','JD')]));
xticks = [juliandate(datetime(tdep(1),'ConvertFrom','JD')+caldays(5))];
for j = 1:nmonthsdep-1
    xticks(j+1) = juliandate(datetime(tdep(1),'ConvertFrom','JD')+caldays(5)+calmonths(j));
end
xticklabels = {};
for j = 1:length(xticks)
    xticklabels{j} = char(datetime(xticks(j),'ConvertFrom','JD','Format','dd-MMM'));
end

nmonthsarr = calmonths(caldiff([datetime(tarr(1),'ConvertFrom','JD'),datetime(tarr(end),'ConvertFrom','JD')]));
yticks = [juliandate(datetime(tarr(1),'ConvertFrom','JD')+caldays(5))];
for j = 1:nmonthsarr-1
    yticks(j+1) = juliandate(datetime(tarr(1),'ConvertFrom','JD')+caldays(5)+calmonths(j));
end
yticklabels = {};
for j = 1:length(yticks)
    yticklabels{j} = char(datetime(yticks(j),'ConvertFrom','JD','Format','dd-MMM'));
end
%finished conversion

%displays dashed line contours with x = departure time, y = arrival time, z = duration
contour(tdep,tarr,dta,20,'k--','ShowText','on'); 
hold on

%displays contours with x = departure time, y = arrival time, z = min dv
contour(tdep,tarr,dvt,50,'ShowText','on');
set(gca,'FontName','Times','FontSize',16,'XTick',xticks,...
    'XTickLabel',xticklabels,'YTick',yticks,'YTickLabel',yticklabels);
xlabel('Time of Departure (2035)');
ylabel('Time of Arrival (2035)');
title('Delta-V (Earth-Mars Transfer)');
hold on