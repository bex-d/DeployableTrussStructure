clc;clf;clear
% profile on
%% Variable definition

N = 8;  %Number of sections. must be even
H = 2;  %Height of vertical members, m
R = 10; %Radius of outer ring, m
R_dash = 0.6; %initial (undeployed) radius, m

mass(1:4)= 1; mass(5:9)=  1;     %mass per length, kg/m.  truss types 1-4 are fixed length (steel or carbon fibre), 5-8 are variable (elastic)
E(1:4)=100;  E(5:9)=10;        %Youngs Modulus, Pa
A(1:4)= 1;   A(5:9)= 1;       %CSA, m2

C = 0.7;        %damping kg/s   

dt = 0.01;      %time step, must be <0.5
tf = 1000;      %final time
t=1:dt:tf+1 ;   %time

%%

[X0,~,l0]=structure(N,H,R,R);        %Gets the lengths of the trusses when the structure is fully deployed
[X,iX,~]=structure(N,H,R,R_dash);   %gets the node coordinates and connectivity matrix for the inital, undeployed model

M=massvector(X,iX,mass);           %Creates mass vector from the mass of each truss type

x0 = X(:,1); %define initial coordinates
y0 = X(:,2);
z0 = X(:,3);
initlength = 0;

nTrusses = length(iX);
nNodes  = length(X);
nTime  = length(t);
L0 = zeros(length(E),1);

for i=1:nTrusses
    L0(iX(i,3)) = l0(i);                                                %assign final length for truss types
end

x = zeros(nNodes,nTime); y = x; z = x;                                  %preallocated empty matrix

x(1:end,1:2) = [x0,x0];
y(1:end,1:2) = [y0,y0];
z(1:end,1:2) = [z0,z0];

myPlot(X,iX);
xlabel('X axis');ylabel('Y axis');zlabel('Z Axis');
title('Undeployed State');

waitforbuttonpress

title('Deploying'); 

for i=2:nTime-1              % 1:tf/dt+1            %Results Matrix
    
    t(i)-1 %prints the time
    F_x = zeros(nTrusses,1); F_y = F_x;  F_z = F_x;
    Fx = zeros(nNodes,1);    Fy = Fx;    Fz = Fx;

    for j=1:nTrusses
        
        n1  = iX(j,1);              %node 1
        n2  = iX(j,2);              %node 2
        typ = iX(j,3);              %truss type

        delta_x = x(n2,i)-x(n1,i);  %calc new changes in length
        delta_y = y(n2,i)-y(n1,i);
        delta_z = z(n2,i)-z(n1,i);
        

        L(i,j) = sqrt(delta_x^2 + delta_y^2 + delta_z^2);
        

        %        
        %         if  typ == 5                %type 5 truss is a rope with applied force.
        %             
        %             F(i,typ) =0.03;
        %             if i ==2
        %                 initlength = initlength + L(i,j);
        %             end
        %           
        %         elseif typ ==6
        %             
        %             F(i,typ) = E(typ)*A(typ)*log(L(i,j)/L0(typ)); 
        %             if F(i,typ)<0
        %                 F(i,typ)=0;
        %             end
        %        
        %         else   

        F(i,typ) = E(typ)*A(typ)*log(L(i,j)/L0(typ));  

        %         end

        F_x = F(i,typ) * delta_x/L(i,j);
        F_y = F(i,typ) * delta_y/L(i,j);
        F_z = F(i,typ) * delta_z/L(i,j);

        Fx(n1) = Fx(n1) - F_x;   %allocate forces from trusses to nodes
        Fy(n1) = Fy(n1) - F_y;
        Fz(n1) = Fz(n1) - F_z;

        Fx(n2) = Fx(n2) + F_x;
        Fy(n2) = Fy(n2) + F_y;
        Fz(n2) = Fz(n2) + F_z;
     
    end

    x(:,i+1) = 1./(2.*M + C*dt)  .*  ( C*dt*x(:,i-1) + 2*M.* ((2*x(:,i))-x(:,i-1))  - 2*dt^2*Fx );      %dynamic response (using central difference scheme)
    y(:,i+1) = 1./(2.*M + C*dt)  .*  ( C*dt*y(:,i-1) + 2*M.* ((2*y(:,i))-y(:,i-1))  - 2*dt^2*Fy );      %WITH DAMPING
    z(:,i+1) = 1./(2.*M + C*dt)  .*  ( C*dt*z(:,i-1) + 2*M.* ((2*z(:,i))-z(:,i-1))  - 2*dt^2*Fz );

    X = [x(:,i+1) y(:,i+1) z(:,i+1)];
    
    pause(0.00001); %10ms
    clf;
    myPlot(X,iX)
    xlim([-R*1.5 R*1.5]);ylim([-R*1.5 R*1.5]);zlim([-6*H 6*H]);
        
%     if  round(x(1,i+1),4) == round(x(1,i),4) ==  round(x(1,i-1),4)
%         pause(0.00001); %10ms
%         clf;
%         myPlot(X,iX)
%         xlim([-R*1.5 R*1.5]);ylim([-R*1.5 R*1.5]);zlim([-6*H 6*H]);
%        
%       %  E = F(i,5)*initlength
%     
%         filename = 'C:\Users\Owner\Downloads\diss\results\24.04\FORCE dt0.001 c0.7 F0.03.mat';    save( filename, 'F' );
%         filename = 'C:\Users\Owner\Downloads\diss\results\24.04\LENGTH dt0.001 c0.7 F0.03.mat';    save( filename, 'L' );
%         filename = 'C:\Users\Owner\Downloads\diss\results\24.04\COORD dt0.001 c0.7 F0.03.mat';    save( filename, 'x' );
%        % MaxF = maxforces(F)
%         
%         return
%         
%     else
%         if t(i)-1 == 999.9
% %         filename = 'C:\Users\Owner\Downloads\diss\results\24.04\COORD dt0.001 c0.7 F0.03.mat';    save( filename, 'x' );
% %         filename = 'C:\Users\Owner\Downloads\diss\results\24.04\FORCE dt0.001 c0.7 F0.03.mat';    save( filename, 'F' );
% %         filename = 'C:\Users\Owner\Downloads\diss\results\24.04\LENGTH dt0.001 c0.7 F0.03.mat';    save( filename, 'L' );
% % %         fintime = finishTime(x)
% %        % MaxF = maxforces(F)
%        
%     end

end
% profile off
% profile viewer