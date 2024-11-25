% Author: Reagan Gakstatter / reg0052@auburn.edu
% Date: 2024-12-02
% Assignment Name: hw06


classdef hw06
    methods (Static)
        
       function ret = p1(func, a, b, n, option)
            % Implements composite quadrature rules for integration.
            %
            % :param func: The function to integrate, provided as a function handle.
            % :param a: The lower bound of the integration interval.
            % :param b: The upper bound of the integration interval.
            % :param n: The number of subintervals to use.
            % :param option: The quadrature rule (1: midpoint, 2: trapezoidal, 3: Simpson).
            % :return: Approximation of the integral over [a, b].
        
            h = (b - a) / n; % Step size
        
            if option == 1
                % Midpoint Rule
                midpoints = a + h/2 : h : b - h/2; % Midpoints of each subinterval
                ret = h * sum(func(midpoints));
            
            elseif option == 2
                % Trapezoidal Rule
                x = linspace(a, b, n + 1); % Subinterval endpoints
                y = func(x);
                ret = h * (0.5 * y(1) + sum(y(2:end-1)) + 0.5 * y(end));
            
            elseif option == 3
                % Simpson's Rule (n must be even)
                if mod(n, 2) ~= 0
                    error('Simpson''s rule requires an even number of subintervals (n).');
                end
                x = linspace(a, b, n + 1); % Subinterval endpoints
                y = func(x);
                ret = h / 3 * (y(1) + 4 * sum(y(2:2:end-1)) + 2 * sum(y(3:2:end-2)) + y(end));
            
            else
                error('Invalid option. Must be 1 (midpoint), 2 (trapezoidal), or 3 (Simpson).');
            end
        end

        
  
        function p2()
            % run with the following command: hw06.p2(). Do not edit this function.
            %
            % It checks the convergence of the composite quadrature rules implemented in p1.
            %
            % Here we use some examples, 
            % f_1(x) = exp(x) 
            % f_2(x) = (1 - x^2)^3, this function's first 2 derivatives at the endpoints are zero.
            % f_3(x) = (1 - x^2)^5, this function's first 4 derivatives at the endpoints are zero.
            % f_4(x) = (1 - x^2)^7, this function's first 6 derivatives at the endpoints are zero.

            % Run this function will plot the figures for the convergence of the composite quadrature rules.
            % Make comments about the results obtained from the plots. 
            %
            % > For instance, you can comment on the convergence rates of the quadrature rules, and how they compare to the theoretical rates.
            % > Here are a few example questions you can answer in your comments:
            % > Does the Runge phenomenon of f1 (Runge's function) lower the convergence rate?
            % > Does Simpson's rule have a faster convergence rate than the other two for these examples?
            % > Based on your observations, what kind of functions can have a faster convergence rate with the given composite quadrature rules?

            % Write your comments here.
            %
            % The midpoint and trapezoidal rules converged at the
            % second-order for all test functions. Simpson's rule converged
            % at the fourth-order and demonstrated better accuracy. The
            % exponential function f1x = exp(x) converged slower due to
            % Runge phenomenon. The polynomial functions had increasing
            % smoothness at endpoints to allow for smaller error and faster
            % convergence. The functions that experienced higher smoothness
            % and derivative continuity at endpoints are preferred.            
            %

            f = {  @(x)exp(x),  @(x) (1 - x.^2 ).^3, @(x)(1 - x.^2).^5,  @(x) (1 - x.^2).^7} ;  % Define the integrand
            exact = [exp(1) - exp(-1), 32/35, 512/693 , 4096/6435];  % Define the exact integral
            n = 2.^(1:8);  % Define the number of subintervals
            for k = 1 : length(f)

                error = zeros(3, length(n));  % Initialize the error matrix with zeros

                % Calculate the approximate integral and the error for each quadrature rule and number of subintervals
                for i = 1 : length(n)
                    error(1, i) = abs(hw06.p1(f{k},-1, 1, n(i), 1) - exact(k));
                    error(2, i) = abs(hw06.p1(f{k},-1, 1, n(i), 2) - exact(k));
                    error(3, i) = abs(hw06.p1(f{k},-1, 1, n(i), 3) - exact(k));
                end

                % Plot the error against the number of subintervals using a log-log scale
                figure(k);
    
                loglog(n, error(1, :), 'r-+', 'LineWidth', 2);
                hold on;
                loglog(n, error(2, :), 'g-d', 'LineWidth', 2);
                loglog(n, error(3, :), 'b-x', 'LineWidth', 2);

                loglog(n, 1./ n.^2, 'm--', 'LineWidth', 1);
                loglog(n, 1./ n.^4, 'k-.', 'LineWidth', 1);
                loglog(n, 1./ n.^6, 'm--d', 'LineWidth', 1);
                loglog(n, 1./ n.^8, 'k--o', 'LineWidth', 1);

                xlabel('Number of subintervals');
                ylabel('Absolute error');
                title(sprintf('Convergence of composite quadrature rules for %s', functions(f{k}).function));
                legend('Midpoint rule', 'Trapezoidal rule', 'Simpson''s rule', '2nd order convergence', '4th order convergence', '6th order convergence', '8th order convergence', 'Location', 'best');
                grid on;
                hold off;
            end

        end

        
        function ret = p3(func, a, b, max_refinements, option)
            n = 1; % Start with 1 interval
            tol = 1e-8; % Numerical tolerance
            max_iters = 10; % Prevent infinite loop
            
            for iter = 1:max_iters
                % Initialize Romberg table
                R = zeros(max_refinements, max_refinements);
                
                % Fill the first column using p1
                for i = 1:max_refinements
                    R(i, 1) = hw06.p1(func, a, b, n, option);
                    n = n * 2; % Refine by doubling intervals
                end
                
                % Fill the Romberg table
                for j = 2:max_refinements
                    for i = j:max_refinements
                        R(i, j) = (4^(j-1) * R(i, j-1) - R(i-1, j-1)) / (4^(j-1) - 1);
                    end
                end
                
                % Check for convergence
                if abs(R(max_refinements, max_refinements) - R(max_refinements-1, max_refinements-1)) < tol
                    break;
                end
            end
            
            ret = R(max_refinements, max_refinements);
        end



        function ret = p4()
            % Gauss quadrature rule with 6 nodes
            P6 = @(x) (1/16) * (231 * x.^6 - 315 * x.^4 + 105 * x.^2 - 5); % Legendre polynomial P6
            dP6 = @(x) (1/16) * (1386 * x.^5 - 1260 * x.^3 + 210 * x); % Derivative of P6
            
            roots = zeros(6, 1);
            weights = zeros(6, 1);
            intervals = [-1, -3/4; -3/4, -1/4; -1/4, 0; 0, 1/4; 1/4, 3/4; 3/4, 1];
            
            for i = 1:6
                % Find roots using bisection
                a = intervals(i, 1);
                b = intervals(i, 2);
               while abs(b - a) > 1e-14
                        c = (a + b) / 2;
                        if P6(c) == 0
                            break;
                        elseif P6(a) * P6(c) < 0
                            b = c;
                        else
                            a = c;
                        end
                    end
                    roots(i) = (a + b) / 2;
                end 
                % Compute weights
                for i = 1:6
                    weights(i) = 2 / ((1 - roots(i)^2) * dP6(roots(i))^2);
                end
                
                ret = [roots, weights];
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %                                                                                                             %
        % Helper functions below. Do not modify. You can create your own helper functions if needed.                  %
        %                                                                                                             %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % Helper functions for p4. The following function is used to evaluate the Legendre polynomial of degree 6.
        function val = legendre_poly_6(x)
            % Compute the Legendre polynomial of degree 6 at the point x.
            %
            % :param x: The point at which to evaluate the Legendre polynomial.
            % :return: The value of the Legendre polynomial of degree 6 at the point x.

            val = (231 * x^6 - 315 * x^4 + 105 * x^2 - 5) / 16;
        end

        % Helper functions for p5. The following function is used to evaluate the Legendre polynomial of degree n.
        function val = legendre_poly(n, x)
            % Compute the nth Legendre polynomial P_n at the point x.
            %
            % :param n: The degree of the Legendre polynomial.
            % :param x: The point at which to evaluate the Legendre polynomial.
            % :return: The value of the nth Legendre polynomial at the point x.

            if (n == 0)
                val = 1;
            elseif (n == 1)
                val = x;
            else
                val = hw06.legendre_poly(n-1, x) * x * (2 * n - 1)/n - (n - 1) * hw06.legendre_poly(n - 2, x) / n;
            end
        end

        function val = deriv_lengendre_poly(n, x)
            % Compute the derivative of the nth Legendre polynomial P_n at the point x.
            %   
            % :param n: The degree of the Legendre polynomial.
            % :param x: The point at which to evaluate the derivative of the Legendre polynomial.
            % :return: The value of the derivative of the nth Legendre polynomial at the point x.
            val = n / (x^2 - 1) * (x * hw06.legendre_poly(n, x) - hw06.legendre_poly(n - 1, x));
        end
    end
end
