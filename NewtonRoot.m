function x_sol = NewtonRoot(fun,fun_der,x_guess,eps,imax)

    % NewtonRoot finds the root of fun = 0 near the point x_guess using Newton's method.
    % Input variables:
    %   fun       Name (string) of a function file that calculates 'fun' for a given x.
    %   fun_der   Name (string) of a function file that calculates the derivative of 'fun' for given x
    %   x_guess   Initial estimate of the solution
    %   eps       Maxiumum error
    %   imax      Maxiumum number of iterations
    % Output variable:
    %   x_sol        Solution

    x(1)=x_guess;
    
    for i=2:imax+1
        x(i) = x(i-1) - feval(fun,x(i-1))/feval(fun_der,x(i-1));
        if abs((x(i)-x(i-1))/x(i-1)) < eps
            x_sol = x(i);
            break
        end
    end
    
    if i==imax+1
        error('Solution was not obtained in %i iterations.\n',imax);
    end
    
end