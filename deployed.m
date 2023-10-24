clc;clf;clear

N = 40;  %Number of sections. Even,
H = 2;  %Height of vertical members, m
R = 20; %Radius of outer ring, m
R_dash = 20;

mass(1:4)= 1; mass(5:9)=  1;     %mass per length, kg/m.  truss types 1-4 are fixed length, 5-8 are variable
E(1:4)=100;  E(5:9)=10;        %Youngs Modulus, Pa
A(1:4)= 1;   A(5:9)= 1;       %CSA, m2
   
dt = 0.1;       %time step <0.5
tf = 1000;         %final time
t=1:dt:tf+1 ;   %time

C = 0.7;          %no damping Ns/m kg/s

[X,iX,l0]=structure(N,H,R,R_dash);        %Gets the lengths of the trusses when the structure is fully deployed

M=massvector(X,iX,mass);           %Creates mass vector from the mass of each truss type

x0 = X(:,1);
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
