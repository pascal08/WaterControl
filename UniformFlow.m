function d = UniformFlow(y_guess,Q,S0,b,m,k)

    % UniformFlow finds the equilibrium between bed- and friction slope.
    % Input variables:
    %   y        %Initial guess of uniform water depth  
    %   Q        %Initial flow rate [m^3/s]
    %   S0       %Bed slope of canal [-]
    %   b        %Channel width b [m]
    %   m        %Side slope [-]
    %   k        %Strickler coefficient [m^(1/3)/s] (1/n where n is Manning n)
    % Output variable:
    %   d        %Uniform water depth

    A = @(y) (b+m*y)*y;            % Cross-sectional area of water in canal
    P = @(y) b+2*y*(1+m^2)^0.5;    % Wet perimeter of water in canal
    Sf = @(y) Q*abs(Q)*P(y)^(4/3)/(k^2*A(y)^(10/3));  % Friction slope
    F = @(y) S0/Sf(y)-1;         % S0/Sf(y)-1
    
    dAdy = @(y) b+2*m*y;
    dPdy = 2*(1+m^2)^0.5;
    dFdy = @(y) 2/3*(F(y)+1)*(5/A(y)*dAdy(y)-2/P(y)*dPdy);
    
    d = NewtonRoot(F,dFdy,y_guess,1e-5,100);

end
