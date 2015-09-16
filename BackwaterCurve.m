% Matlab Program for the computation of a backwater curve after a weir
clear all

%User input
k=35;       %Strickler coefficient [m^(1/3)/s] (1/n where n is Manning n)
m=2;        %Side slope [-]
b=21;        %Channel width b [m]
Q=80;        %Initial flow rate [m^3/s]
S0=1e-3;    %Bottom slope [-]
beta=1;     %Coeff (?) [-]
y0=1.29;   %Lower boundary condition: water depth given by weir [m]

%Physical constants
g=9.81;     %Gravitational acceleration [m/s^2]

%Estimate of backwater length
length=1000;     %Length of cannal
nx=100;           %Number of discretization intervals
dx=length/nx;    %Spatial discretization lengt

% Direction of computation = direction of information flow: upstream (as flow in backwater is subcritical)

% Intialization of water surface and channel bottom
y(1:nx+1)=y0;   %Water depth
x(1:nx+1)=0;    %Distance from weir 
z(1:nx+1)=0;    %Cannal bottom height relative to downstream height
v(1:nx+1)=0;    %Water velocity
A(1:nx+1)=0;    %Cross-sectional area
dAdy(1:nx+1)=0; %Cross-sectional area
P(1:nx+1)=0;    %Wet perimeter
Sf(1:nx+1)=0;   %Friction slope
Fr2(1:nx+1)=0;  %Froude number (squared)

% Depth uniform flow (S0 = Sf):

% Computation of backwater curve
for i=1:nx+1
    % position along the channel
	x(i)=dx*(i-1);
   
	%Computation of water surface and channel bottom:
	% The h and z at position x+1 are h_old and zold,
	% the h and z at time x are hnew and znew.
     
	if i==1 % Position at the weir:
        y(i)=y0; % constant water depth at the weir
    else % Rest of the channel:
        y(i)=y(i-1)-dx*(S0-Sf(i-1))/(1-Fr2(i-1));
    end
    
    v(i)=Q/(b*y(i));
    A(i)=(b+m*y(i))*y(i);
    dAdy(i)=b+2*m*y(i);
    P(i)=b+2*y(i)*(1+m^2)^0.5;
    Sf(i)=Q*abs(Q)*P(i)^(4/3)/(k^2*A(i)^(10/3));
    Fr2(i)=beta*Q^2*dAdy(i)/(g*A(i)^3);
   
    % z(i)=S0*dx*i;   

end

wl=y+z;             % water level
y(end)

figure
   hold on
   plot(x,wl,'b','LineWidth',1);
   plot(x,z,'r','LineWidth',1);
   xlabel('Distance from weir')
   ylabel('Water elevation (m)')
   title('Backwater curve')
   legend('Water level','Channel Bottom','Location','NorthWest');
   hold off
 