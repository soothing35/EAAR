function [Width,Length,Diff_angle, Energy, Sum_Energy] = AngleEnergy(s,e,a)

close all;
clear all;

% % s = 2.4;
% % e = 10.4;
% % a = 1.2;
% 
%space 4.0 nm; innermost 30 nm; outermost 150 nm;
s = 3.4;
e = 12.715;
a = 1.7;
% % 
% % % %space 8.0 nm; innermost 30 nm; outermost 150 nm;
% % s = 1.7;
% % e = 8.827;
% % a = 0.85;


% 
% %space 11 nm; innermost 30 nm; outermost 150 nm;
% s = 1.24;
% e = 7.5;
% a = 0.62;

% %space 11 nm; innermost 30 nm; outermost 150 nm; One spiral has the same length with double spirals 
% s = 1.24;
% e = 10.536;
% a = 0.62;



bb=2.16;
%% input parameters
j = 1;
step_size = 3.074; %2.9*1.06 switch the pixel size to 1.

k0 = 0.00244; %6.8+-0.79
x0 = 8.75; % 8.92+-1.63


%% %%%%%%DONOT CHANGE ANY THING BELOW%%%%%%
%% judeg c value
if rem(floor(2*s),2) == 0
    c = 2*s-floor(2*s);
else
    c = 2*s-floor(2*s)+1;
end

%% Diameter Analysis
    j=1;

    for i=s+1:0.06:e
        
        tt = linspace(s*pi,i*pi,1000);
        x = bb*tt.*cos(2*tt-c*pi)/a;
        y = bb*tt.*sin(2*tt-c*pi)/a;
        
        L(j) = arclength(x,y);
        
        x1= bb*i*pi*cos(2*i*pi-c*pi)/a;
        y1 = bb*i*pi*sin(2*i*pi-c*pi)/a;
        
        Width_1 = sqrt(x1^2+y1^2);
        
        x2= bb*(i-0.5)*pi*cos(2*(i-0.5)*pi-c*pi)/a;
        y2 = bb*(i-0.5)*pi*sin(2*(i-0.5)*pi-c*pi)/a;
        Width_2 = sqrt(x2^2+y2^2);
        
        Width(j) = (Width_1 + Width_2);
        
        j=j+1;
    end
    
    Length=L'/step_size;
    Width=Width';
    
    figure(1);
    plot(Length,Width,'k');hold on;
    
%% 
    tt = linspace(s*pi,e*pi,1000);
    x = bb*tt.*cos(2*tt-c*pi)/a;
    y = bb*tt.*sin(2*tt-c*pi)/a;
    z = zeros(1,1000);
    L = arclength(x,y);
 
figure(2);
plot(x,y,'r');axis equal;hold on;
%plot3(x,y,z,'r');axis equal;hold on;
%% determine the step size in spacing data

    k(:,1) = x;
    k(:,2) = y;

    [cp1,cp2] = findControlPoints(k); 
    [Inter_points] = displayBezier(k,cp1,cp2);
    
    %% calculate the arclength of the curve
    px = Inter_points(:,1);
    py = Inter_points(:,2);
    [arclen,seglen] = arclength(px,py);


    %% calculate the spacing data

    point_num =  ceil(arclen/step_size);

    [linspace_x,linspace_y] = linspacearc(Inter_points(:,1),Inter_points(:,2),point_num);
    linspace_points = [linspace_x' linspace_y'];

    
    %% calculate the slope for spacing data( tangent line)
    %approximate the derivative (slope) at each point with the following second order differencing
    t = linspace_points(:,1)'; p = linspace_points(:,2)';
    n = size(linspace_points,1);
    td = [t(3),t(1:n-1)]; tu = [t(2:n),t(n-2)];
    pd = [p(3),p(1:n-1)]; pu = [p(2:n),p(n-2)];
    dpdt =(pu-pd)./(tu-td);
    slope = dpdt';
    
    Sum_Energy = [];
    Sum_Energy(1) = 0;


    %% 
    for i=1 : n-1

        if linspace_points(i,1) >= linspace_points((i+1),1)
            ang(i) = 90 -  (-atand(slope(i)));
        else
            ang(i) = 270 - (-atand(slope(i)));
        end
        
        %% determine the possible wrong calculation at angles around 0 or 180
        if i>= 2 && abs(ang(i)-ang(i-1))<= 270 && abs(ang(i)-ang(i-1))>= 90
            ang(i) = ang(i) -180;
            if ang(i) <=0
                ang(i)=ang(i)+360;
            end
        end
        
        ang(i);
        if i >=3
            diff(i-1) = min (abs(ang(i)-ang(i-1)),abs(abs(ang(i)-ang(i-1))-360));
            Energy(i-1) =  k0*(diff(i-1)-x0)^2;
            Sum_Energy(i-1) = Sum_Energy(i-2) + Energy(i-1);
        end
        
    end
Diff_angle = diff(2:end)';
Energy = Energy(2:end)';
Sum_Energy = Sum_Energy(2:end)';
% 
% 
figure(3);
plot(Diff_angle(:,1),'k');hold on;

figure(4);
plot(Energy,'k');hold on;

figure(5);
plot(Sum_Energy,'k');hold on;






