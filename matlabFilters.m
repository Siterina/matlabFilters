classdef matlabFilters
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        %% Kalman Filter
        
        % Functions-coefficients of the system
        
        function [ k ] = k(t) 
            global kb;
            global tStroke;
            k = kb*sign(t - tStroke);
        end
        
        function [ a ] = a(t, x, condition)
            if(condition == 4)
                global a1cond4;
                a = a1cond4 * x;
            end
            if(condition == 3)
                global c3;
                global betaC3;
                global g;
                global R;
                global kb;
                global tStroke;
                k = kb*sign(t - tStroke);
                a(1, 1) = -c3*x(1)*x(1)*exp(-betaC3*x(3)) - g*sin(x(2));
                a(2, 1) = c3*k*x(1)*exp(-betaC3*x(3)) + cos(x(2))*(x(1)/(R + x(3))- g/x(1));
                a(3, 1) = x(1)*sin(x(2));
            end          
            if(condition == 2)
                global omegacond2;
                global alphacond2;
                global betacond2;
                a(1, 1) = x(2);
                a(2, 1) = -omegacond2*omegacond2*x(1) + alphacond2*x(2) - alphacond2*betacond2*x(1)*x(1)*x(2);
            end           
            if(condition == 1)
                global a0cond1;
                global a1cond1;
                global a3cond1;
                a = a0cond1 + a1cond1 * x + a3cond1 * x * x * x; 
            end
        end
        
        function [ a_hut ] = a_hut(t, x, p, condition)    
            if(condition == 4)
                global a1cond4;
                a_hut = a1cond4 * x;
            end
            if(condition == 2)
                global omegacond2;
                global alphacond2;
                global betacond2;
                a_hut(1, 1) = x(2);
                %a_hut(2, 1) = -omegacond2*omegacond2*x(1) + alphacond2*x(2) - alphacond2*betacond2*x(2)*(x(1)*x(1) + p(1, 1)); % из диплома
                a_hut(2, 1) = -omegacond2*omegacond2*x(1) + alphacond2*x(2) - alphacond2*betacond2*(x(1)*x(1)*x(2) + 2*x(1)*p(1, 2)+ x(2)*p(1, 1)); % посчитано
            end           
            if(condition == 1)
                global a0cond1;
                global a1cond1;
                global a3cond1;
                a = a0cond1 + a1cond1 * x + a3cond1 * x * x * x; 
                a_hut = a + 3 * a3cond1 * x * p;
            end
        end
        
        function [ c ] = c(t, x, condition)
            if(condition == 4)
                global c1cond4;
                global c2cond4;
                global omegacond4;
                c = c1cond4 * x + c2cond4 * sin(x + omegacond4*t);
            end
            if(condition == 3)
                global c3;
                global betaC3;
                global kb;
                global tStroke;
                k = kb*sign(t - tStroke);
                c(1, 1) = c3*x(1)*x(1)*exp(-betaC3*x(3))*(cos(x(2)) - k*sin(x(2)));
                c(2, 1) = c3*x(1)*x(1)*exp(-betaC3*x(3))*(sin(x(2)) + k*cos(x(2)));
            end     
            if(condition == 2)
                c(1, 1) = x(1);
                c(2, 1) = x(2);
            end            
            if(condition == 1)
                global c1cond1;
                global c2cond1;
                c = c1cond1 * x + c2cond1 * x * x;
            end
        end
        
        function [ c_hut ] = c_hut(t, x, p, condition) 
            if(condition == 4)
                global c1cond4;
                global c2cond4;
                global omegacond4;
                c_hut = c1cond4*x + c2cond4*exp(-p/2)*sin(x + omegacond4*t);
            end
            if(condition == 2)
                c_hut(1, 1) = x(1);
                c_hut(2, 1) = x(2);
            end            
            if(condition == 1)
                global c1cond1;
                global c2cond1;
                c = c1cond1 * x + c2cond1 * x * x;
                c_hut = c + c2cond1 * p;
            end
        end
        
        
        function [ B ] = B(t, x, condition)
            global gamma;
            if(condition == 4)
                global b1cond4;
                B = b1cond4;
            end
            if(condition == 3)
                B = zeros(3, 2);
            end            
            if(condition == 2)
                B(1, 1) = 0;
                B(1, 2) = 0;
                B(2, 1) = 0;
                B(2, 2) = gamma*x(1);
            end            
            if(condition == 1)
                global b0cond1;
                global b1cond1;
                B = (b0cond1 + b1cond1 * x)*gamma;
            end
        end
        
        function [ Bx ] = Bx(t, x, condition) % B'x
            global gamma;
            if(condition == 1)
                global b1cond1;
                Bx = b1cond1*gamma;
            end
        end
        
        function [ D ] = D(t, x, condition)
            global gamma;
            if(condition == 4)
                global d1cond4;
                D = d1cond4;
            end
            if(condition == 3)
                global c3;
                global betaC3;
                global sigmaM;
                global sigmaA;
                global kb;
                global tStroke;
                k = kb*sign(t - tStroke);
                D(1, 1) = sigmaA;
                D(1, 2) = 0;
                D(2, 1) = 0;
                D(2, 2) = sigmaA;
                
%                 D(1, 1) = sigmaA;
%                 D(1, 2) = 0;
%                 D(1, 3) = c3*x(1)*x(1)*exp(-betaC3*x(3))*(cos(x(2)) - k*sin(x(2)))*sigmaM;
%                 D(1, 4) = 0;
%                 D(2, 1) = 0;
%                 D(2, 2) = sigmaA;
%                 D(2, 3) = 0;
%                 D(2, 4) = c3*x(1)*x(1)*exp(-betaC3*x(3))*(sin(x(2)) + k*cos(x(2)))*sigmaM;
            end
            if(condition == 2)
                global d11cond2;
                global d22cond2;
                D(1, 1) = d11cond2;
                D(1, 2) = 0;
                D(2, 1) = 0;
                D(2, 2) = d22cond2;   
                D = D*gamma;
            end           
            if(condition == 1)
                global d0cond1;
                global d1cond1;
                D = (d0cond1 + d1cond1 * x)*gamma;
            end
        end
        
        function [ Dx ] = Dx(t, x, condition) %D'x
            global gamma;          
            if(condition == 1)
                global d1cond1;
                Dx = d1cond1*gamma;
            end
        end
        
        function [ R ] = R(t, x, condition)
            D = matlabFilters.D(t, x, condition);
            R = D * D.';
        end
        
        function [ Q ] = Q(t, x, condition)
            B = matlabFilters.B(t, x, condition);
            Q = B * B.';
        end
        
        function [ Q_hut ] = Q_hut(t, x, p, condition)
            if(condition == 4)
                global  b1cond4;
                Q_hut = b1cond4 * b1cond4;
            end
            if(condition == 2)
                global gamma;
                Q_hut(1, 1) = 0;
                Q_hut(1, 2) = 0;
                Q_hut(2, 1) = 0;
                Q_hut(2, 2) = gamma*gamma*(x(1)*x(1) + p(1, 1));
            end
            if(condition == 1)
                global  b1cond1;
                B = matlabFilters.B(t, x, condition);
                Q = B * B.';
                Q_hut = Q +  b1cond1 * b1cond1 * p;
            end
        end
        
        
        % Support functions for Kalman filter
        
        function [ A ] = A(t, x, condition)                               % матрица Якоби первых производных a(t, x) по x
            if(condition == 4)
                global a1cond4;
                A = a1cond4;
            end
            if(condition == 3)
                global c3;
                global betaC3;
                global g;
                global R;
                global kb;
                global tStroke;
                k = kb*sign(t - tStroke);
                A(1, 1) = -2*c3*x(1)*exp(-betaC3*x(3));
                A(1, 2) = -g*cos(x(2));
                A(1, 3) = c3*betaC3*x(1)*x(1)*exp(-betaC3*x(3));
                A(2, 1) = c3*k*exp(-betaC3*x(3)) + cos(x(2))*(1/(R + x(3))+ g/(x(1)*x(1)));
                A(2, 2) = -sin(x(2))*(x(1)/(R + x(3))- g/x(1));
                A(2, 3) = -c3*betaC3*k*x(1)*exp(-betaC3*x(3)) - cos(x(2))*x(1)/((R + x(3))*(R+x(3)));
                A(3, 1) = sin(x(2));
                A(3, 2) = x(1)*cos(x(2));
                A(3, 3) = 0;
            end           
            if(condition == 2)
                global omegacond2;
                global alphacond2;
                global betacond2;
                A(1, 1) = 0;
                A(1, 2) = 1;
                A(2, 1) = -omegacond2*omegacond2 - 2*alphacond2*betacond2*x(1)*x(2);
                A(2, 2) = alphacond2*(1 - betacond2*x(1)*x(1));
            end           
            if(condition == 1)
                global a0cond1;
                global a1cond1;
                global a3cond1;
                A = a1cond1 + 3 * a3cond1 * x * x;
            end
        end
        
        function [ A_hut ] = A_hut(t, x, p, condition)                               % матрица Якоби первых производных a_hut(t, x) по x     
            if(condition == 4)
                global a1cond4;
                A_hut = a1cond4;
            end
            if(condition == 2)
                global omegacond2;
                global alphacond2;
                global betacond2;
                A_hut(1, 1) = 0;
                A_hut(1, 2) = 1;
%                 A_hut(2, 1) = -omegacond2*omegacond2 - 2*alphacond2*betacond2*x(1)*x(2); % pdf
%                 A_hut(2, 2) = alphacond2 - alphacond2*betacond2*(x(1)*x(1) + p(1, 1));
                A_hut(2, 1) = -omegacond2*omegacond2 - 2*alphacond2*betacond2*(x(1)*x(2) + p(1, 2));
                A_hut(2, 2) = alphacond2 - alphacond2*betacond2*(x(1)*x(1) + p(1, 1));
            end           
            if(condition == 1)
                global a0cond1;
                global a1cond1;
                global a3cond1;
                A = a1cond1 + 3 * a3cond1 * x * x;
                A_hut = A + 3 * a3cond1 * p;
            end
        end
        
        function [ G ] = G(t, x, condition)                               % матрица Якоби первых производных с(t, x) по x
            if(condition == 4)
                global c1cond4;
                global c2cond4;
                global omegacond4;
                G = c1cond4  + c2cond4 * cos(x + omegacond4*t);
            end
            if(condition == 3)
                global c3;
                global betaC3;
                global kb;
                global tStroke;
                k = kb*sign(t - tStroke);
                G(1, 1) = 2*c3*x(1)*exp(-betaC3*x(3))*(cos(x(2)) - k*sin(x(2)));
                G(1, 2) = -c3*x(1)*x(1)*exp(-betaC3*x(3))*(sin(x(2)) + k*cos(x(2)));
                G(1, 3) = -betaC3*c3*x(1)*x(1)*exp(-betaC3*x(3))*(cos(x(2)) - k*sin(x(2)));
                G(2, 1) = 2*c3*x(1)*exp(-betaC3*x(3))*(sin(x(2)) + k*cos(x(2)));
                G(2, 2) = c3*x(1)*x(1)*exp(-betaC3*x(3))*(cos(x(2)) - k*sin(x(2)));
                G(2, 3) = -betaC3 * c3*x(1)*x(1)*exp(-betaC3*x(3))*(sin(x(2)) + k*cos(x(2)));
            end
            if(condition == 2)
                G(1, 1) = 1;
                G(1, 2) = 0;
                G(2, 1) = 0;
                G(2, 2) = 1;
           end            
           if(condition == 1)
                global c1cond1;
                global c2cond1;
                G = c1cond1 + 2 * c2cond1 * x;
           end
        end
        
        function [ G_hut ] = G_hut(t, x, p, condition)                               % матрица Якоби первых производных с_hut(t, x) по x
            if(condition == 4)
                global c1cond4;
                global c2cond4;
                global omegacond4;
                G_hut = c1cond4 + c2cond4*exp(-p/2)*cos(x + omegacond4*t);
            end
            if(condition == 2)
                G_hut(1, 1) = 1;
                G_hut(1, 2) = 0;
                G_hut(2, 1) = 0;
                G_hut(2, 2) = 1;
           end         
           if(condition == 1)
                global c1cond1;
                global c2cond1;
                G = c1cond1 + 2 * c2cond1 * x;
                G_hut = G;
           end
        end
        
        function [ Theta ] = Theta(t, x, p, condition)
            if(condition == 4)
                global c2cond4;
                global omegacond4;
                Theta = -p*p* c2cond4*exp(-p/2)*sin(x + omegacond4*t);
            end
            if(condition == 2)
                Theta(1, 1) = 0;
                Theta(1, 2) = 0;
           end 
            if(condition == 1)
                global c2cond1;
                Theta = 2 * c2cond1 * p * p;
            end
        end
        
        
        % structure functions
        
        function [ K ] = K(t, x, p, condition) 
            K = p * matlabFilters.G(t, x, condition).' /(matlabFilters.R(t, x, condition));
        end
        
        function [ K_hut ] = K_hut(t, x, p, condition) 
            K_hut = p * matlabFilters.G_hut(t, x, p, condition).' /(matlabFilters.R(t, x, condition));
        end
        
        function [ Ksi ] = Ksi(t, x, p, condition)                        
            Ksi = matlabFilters.A(t, x, condition) * p + p * matlabFilters.A(t, x, condition).' + ...
                matlabFilters.Q(t, x, condition) - matlabFilters.K(t, x, p, condition) * ... 
                matlabFilters.R(t, x, condition) * matlabFilters.K(t, x, p, condition).';
        end
        
        function [ Ksi_hut ] = Ksi_hut(t, x, p, condition)
            
            Ksi_hut = matlabFilters.A_hut(t, x, p, condition) * p + p * matlabFilters.A_hut(t, x, p, condition).' + ...
                matlabFilters.Q_hut(t, x, p, condition) - matlabFilters.K_hut(t, x, p, condition) * ... 
                matlabFilters.R(t, x, condition) * matlabFilters.K_hut(t, x, p, condition).';
        end
        
        function [ Xi ] = Xi(t, x, p, condition) 
            Xi = matlabFilters.K(t, x, p, condition);
        end
        
        function [ Xiz ] = Xiz(t, x, p, condition) 
            global d0cond1;
            global d1cond1;
            global c1cond1;
            global c2cond1;
            Xiz = p * (-2)*(-c2cond1*d0cond1+c2cond1*d1cond1*x+d1cond1*c1cond1)/(d0cond1+d1cond1*x)^3;
        end
        
        function [ Zeta ] = Zeta(t, x, p, condition) 
            Zeta = matlabFilters.a(t, x, condition) - matlabFilters.Xi(t, x, p, condition)*matlabFilters.c(t, x, condition);
        end
        
        function [ Pi ] = Pi(t, x, z, p, condition) 
            Pi = matlabFilters.Xi(t, z, p, condition)*matlabFilters.Dx(t, x, condition)*matlabFilters.B(t, x, condition) + ...
                matlabFilters.D(t, x, condition)*matlabFilters.Xiz(t, z, p, condition)* ...
                matlabFilters.Xi(t, z, p, condition)*matlabFilters.D(t, x, condition);
        end
        
        % Support functions
        
        function [ M ] = changeMx(i, x, Mx)
            M = Mx + (x - Mx)/i;
        end
        
        function [ D ] = changeDx(i, x, Mx, Dx)
            D = Dx + ((x - Mx) * (x - Mx) * (i - 1)/i - Dx)/i;
        end
        
        function [ M ] = mathExpectation(Temp, NX)
            global N;          
            M = zeros(NX, 1);
            for i = 1:N
                M = M + Temp(:, i);
            end
            M = M * 1/N;
        end
        
        
        function [ D ] = dispersion(Temp, M, NX)
            global N;
            if(NX == 1)
                D = 0;
                for i = 1:N
                    D = D + (Temp(i) - M)*(Temp(i) - M);
                end
                D = D * 1/(N-1);
            end
            if(NX > 1)
                D = zeros(NX, NX);
                for i = 1:N
                    for j = 1:NX
                        for k = 1:NX
                            D(j, k) = D(j, k) + (Temp(j, i) - M(j))*(Temp(k, i) - M(k));
                        end
                    end
                end
                for j = 1:NX
                    for k = 1:NX
                        D(j, k) = D(j, k)/(N-1);
                    end
                end
            end   
        end
        
        
        
        
        function [ ] = initialConditions(a0cond1G, a1cond1G, a3cond1G, b0cond1G, b1cond1G, c1cond1G, c2cond1G, d0cond1G, d1cond1G, MoXcond1G, SoXcond1G, ...
                                        omegacond2G, alphacond2G, betacond2G, d11cond2G, d22cond2G, SoXcond2G, SoYcond2G, MoXcond2G, MoYcond2G, gammaG, ...
                                        a1cond4G, b1cond4G, c1cond4G, c2cond4G, d1cond4G, omegacond4G, MoXcond4G, SoXcond4G)
                        
            global deltaTime;
            global N;
            global K;
            
            %1 condition
            global a0cond1;
            global a1cond1;
            global a3cond1;
            global c1cond1;
            global c2cond1;
            global b0cond1;
            global b1cond1;
            global d0cond1;
            global d1cond1;
            
            global MoXcond1;
            global MoYcond1;
            global SoXcond1;
            global SoYcond1;
            global SoZcond1;
            
            %2 condition
            global omegacond2;
            global alphacond2;
            global betacond2;
            global d11cond2;
            global d22cond2;
            global SoXcond2;
            global SoYcond2;
            global gamma;
            global MoXcond2;
            global MoYcond2;
            
            %3 condition
            global c3;
            global betaC3;
            global g;
            global R;
            global tStroke;
            global kb;
            global MoXv;
            global MoXteta;
            global MoXh;
            global SoXv;
            global SoXteta;
            global SoXh;
            global sigmaM;
            global sigmaA;
            
            %4 condition
            global a1cond4;
            global b1cond4;
            global c1cond4;
            global c2cond4;
            global d1cond4;
            global omegacond4
            global MoXcond4;
            global SoXcond4;
            
            
            %1 condition
            a0cond1 = a0cond1G;                                                        % a, g - параметры объекта
            a1cond1 = a1cond1G;
            a3cond1 = a3cond1G;
            c1cond1 = c1cond1G;                                                    % alpha, beta - параметры измерителя
            c2cond1 = c2cond1G;
            b0cond1 = b0cond1G;
            b1cond1 = b1cond1G;
            d0cond1 = d0cond1G;
            d1cond1 = d1cond1G;
                        
            MoXcond1 = MoXcond1G;                                                     % Характеристики нач. усл. объекта, измерителя и фильтра
            MoYcond1 = 0.0;
            SoXcond1 = SoXcond1G;
            SoYcond1 = 0.0;
            SoZcond1 = 0.0;

            % easy scalar
%             a0cond1 = 0;                                                      
%             a1cond1 = 0;
%             a3cond1 = -1;
%             c1cond1 = 1;                                                   
%             c2cond1 = 0;
%             b0cond1 = 0;
%             b1cond1 = 1;
%             d0cond1 = 1;
%             d1cond1 = 0;           
% 
%             MoXcond1 = 0.1;                                                     
%             MoYcond1 = 0.0;
%             SoXcond1 = 0.707;
%             SoYcond1 = 0.0;
%             SoZcond1 = 0.0; 

            %2 condition
            omegacond2 = omegacond2G;
            alphacond2 = alphacond2G;
            betacond2 = betacond2G;
            d11cond2 = d11cond2G;
            d22cond2 = d22cond2G;
            SoXcond2 = SoXcond2G;
            SoYcond2 = SoYcond2G;
            MoXcond2 = MoXcond2G;
            MoYcond2 = MoYcond2G;
            gamma = gammaG;
             
            
            %cond 3 - meters
            c3 = 0.0043333333;
            betaC3 = 0.00009;
            g = 3.72;
            R = 3400000;
            tStroke = 45;
            kb = 0.3;
            MoXv = 6000;
            MoXteta = -18*pi/180;
            MoXh = 100000;
            SoXv = 15;
            SoXteta = 1*pi/180;
            SoXh = 7000;
            sigmaM = 0;
            sigmaA = 200;
            
            
            %4 condition
            a1cond4 = a1cond4G;
            b1cond4 = b1cond4G;
            c1cond4 = c1cond4G;
            c2cond4 = c2cond4G;
            d1cond4 = d1cond4G;
            omegacond4 = omegacond4G;
            MoXcond4 = MoXcond4G;
            SoXcond4 = SoXcond4G;
            
           
        end
        
        
        function[ ] = main(condition, N1, K1, deltaTime1, buildTrajectory, trajectoryNumber, buildLAOF, buildLFOS, buildGAOF, buildGFOS, buildMx, useEuler, useHun)
            
            % Initial conditions
            global MoXcond1;
            global MoYcond1;
            global SoXcond1;
            global SoYcond1;
            global SoZcond1;
            global MoXcond4;
            global SoXcond4;
            
            global SoXcond2;
            global SoYcond2;
            global MoXcond2;
            global MoYcond2;
            
            global MoXv;
            global MoXteta;
            global MoXh;
            global SoXv;
            global SoXteta;
            global SoXh;
            global N;
            global K;
            global deltaTime;
            N = N1;
            K = ceil(K1 / deltaTime1);
            deltaTime = deltaTime1;

%             matlabFilters.initialConditions();
            if(condition == 4)
                NX = 1;
                NY = 1;
                NW = 1;
           end 
            if(condition == 3)
                NX = 3;
                NY = 2;
                NW = 2;
            end
            if(condition == 2)
                NX = 2;
                NY = 2;
                NW = 2;
            end
           if(condition == 1)
                NX = 1;
                NY = 1;
                NW = 1;
           end            
                        
            % arrays for final table            
            tTable = ones(1, 1, K);
            MxTable = ones(NX, 1, K);
            SxTable = ones(NX, 1, K);
            Sx_RTable = ones(NX, 1, K);
            Sx_HTable = ones(NX, 1, K);
            Meps_aofLTable = ones(NX, 1, K);
            Seps_aofLTable = ones(NX, 1, K);
            Crit_aofLTable = ones(NX, 1, K);
            ICrit_aofLTable = ones(NX, 1, K);
            Meps_fosLTable = ones(NX, 1, K);
            Seps_fosLTable = ones(NX, 1, K);
            Crit_fosLTable = ones(NX, 1, K);
            ICrit_fosLTable = ones(NX, 1, K);
            Meps_aofGTable = ones(NX, 1, K);
            Seps_aofGTable = ones(NX, 1, K);
            Crit_aofGTable = ones(NX, 1, K);
            ICrit_aofGTable = ones(NX, 1, K);
            Meps_fosGTable = ones(NX, 1, K);
            Seps_fosGTable = ones(NX, 1, K);
            Crit_fosGTable = ones(NX, 1, K);
            ICrit_fosGTable = ones(NX, 1, K);
            
            t = deltaTime;
            X = ones(NX, N);
            X_R = ones(NX, N);
            X_H = ones(NX, N);
            X_tilde = ones(NX, N);
            deltaY = ones(NY, N);
            deltaY_H = ones(NY, N);
            Z_aofL = ones(NX, N);
            Z_fosL = ones(NX, N);
            Z_aofL_H = ones(NX, N);
            Z_fosL_H = ones(NX, N);
            Z_fosL_H_tilde = ones(NX, N);
            Z_aofL_H_tilde = ones(NX, N);
            Z_aofG = ones(NX, N);
            Z_fosG = ones(NX, N);
            P_L = ones(NX, NX, N);
            P_L_H = ones(NX, NX, N);
            P_L_H_tilde = ones(NX, NX, N);
            P_G = ones(NX, NX, N);
            J_L = ones(NX, NX, N);
            J_L_tilde = ones(NX, NX, N);
            J_G = ones(NX, NX, N);
                      
            
%             t = zeros(1,100);
%             for n = 1:100
%                 A = rand(n,n);
%                 B = rand(n,n);
%                 tic;
%                 C = A*B;
%                 t(n) = toc;
%             end
%             plot(t)

            
            rng(10);
            if(condition == 4)
                for i = 1:N
                    X(:, i) = normrnd(MoXcond4, SoXcond4);
                    Z_aofL(:, i) = MoXcond4;
                    Z_fosL(:, i) = Z_aofL(:, i);
                    Z_aofG(:, i) = Z_aofL(:, i);
                    Z_fosG(:, i) = Z_aofL(:, i);
                    P_L(:, :, i) = SoXcond4 * SoXcond4;
                    P_G(:, :, i) = P_L(:, :, i);
                    J_L(:, :, i) = P_L(:, :, i);
                    J_G(:, :, i) = P_L(:, :, i);
                end
            end
            
            if(condition == 3)
                 for i = 1:N;
                    P_L(:, :, i) = [SoXv*SoXv, 0, 0; 0, SoXteta*SoXteta, 0; 0, 0, SoXh*SoXh];
                    P_G(:, :, i) = P_L(:, :, i);
                    J_L(:, :, i) = P_L(:, :, i);
                    J_G(:, :, i) = P_L(:, :, i);
                    Z_aofL(:, i) = [MoXv, MoXteta, MoXh];
                    Z_fosL(:, i) = Z_aofL(:, i);
                    Z_aofG(:, i) = Z_aofL(:, i);
                    Z_fosG(:, i) = Z_aofL(:, i);

                    X(1, i) = normrnd(MoXv, SoXv);
                    X(2, i) = normrnd(MoXteta, SoXteta);
                    X(3, i) = normrnd(MoXh, SoXh);
                 end
            end
            
            if(condition == 2)
                 for i = 1:N;
                    P_L(:, :, i) = [SoXcond2*SoXcond2, 0; 0, SoYcond2*SoYcond2];
                    P_G(:, :, i) = P_L(:, :, i);
                    J_L(:, :, i) = P_L(:, :, i);
                    J_G(:, :, i) = P_L(:, :, i);
                    Z_aofL(:, i) = [MoXcond2, MoYcond2];
                    Z_fosL(:, i) = Z_aofL(:, i);
                    Z_aofG(:, i) = Z_aofL(:, i);
                    Z_fosG(:, i) = Z_aofL(:, i);

                    X(1, i) = normrnd(MoXcond2, SoXcond2);
                    X(2, i) = normrnd(MoYcond2, SoYcond2);
                 end
            end
            
            if(condition == 1)
                for i = 1:N
                    X(:, i) = normrnd(MoXcond1, SoXcond1);
                    X_R(:, i) = X(:, i);
                    X_H(:, i) = X(:, i);
                    Z_aofL(:, i) = normrnd(MoXcond1, SoZcond1);
                    Z_fosL(:, i) = Z_aofL(:, i);
                    Z_aofL_H(:, i) = Z_aofL(:, i);
                    Z_fosL_H(:, i) = Z_aofL(:, i);
                    Z_aofG(:, i) = Z_aofL(:, i);
                    Z_fosG(:, i) = Z_aofL(:, i);
                    P_L(:, :, i) = SoXcond1 * SoXcond1;
                    P_G(:, :, i) = P_L(:, :, i);
                    P_L_H(:, :, i) = P_L(:, :, i);
                    P_L_H_tilde(:, :, i) = P_L(:, :, i);
                    J_L(:, :, i) = P_L(:, :, i);
                    J_L_tilde(:, :, i) = P_L(:, :, i);
                    J_G(:, :, i) = P_L(:, :, i);
                end
            end
            
            % step #0
            Mx = matlabFilters.mathExpectation(X, NX);
            Mx_R = Mx;
            Mx_H = Mx;
            Meps_aofL = matlabFilters.mathExpectation(X - Z_aofL, NX);
            Deps_aofL = matlabFilters.dispersion(X - Z_aofL, Meps_aofL, NX);
            Meps_fosL = matlabFilters.mathExpectation(X - Z_fosL, NX);
            Deps_fosL = matlabFilters.dispersion(X - Z_fosL, Meps_fosL, NX);
            Meps_aofL_H = matlabFilters.mathExpectation(X - Z_aofL_H, NX);
            Deps_aofL_H = matlabFilters.dispersion(X - Z_aofL_H, Meps_aofL_H, NX);
            Meps_fosL_H = matlabFilters.mathExpectation(X - Z_fosL_H, NX);
            Deps_fosL_H = matlabFilters.dispersion(X - Z_fosL_H, Meps_fosL_H, NX);
            Meps_aofG = matlabFilters.mathExpectation(X - Z_aofG, NX);
            Deps_aofG = matlabFilters.dispersion(X - Z_aofG, Meps_aofG, NX);
            Meps_fosG = matlabFilters.mathExpectation(X - Z_fosG, NX);
            Deps_fosG = matlabFilters.dispersion(X - Z_fosG, Meps_fosG, NX);
            
            Dx = matlabFilters.dispersion(X, Mx, NX);
            Dx_R = Dx;
            Dx_H = Dx;
            Dx_H_tilde = Dx;
            Mz_fosL = matlabFilters.mathExpectation(Z_fosL, NX);
            Dz_fosL = matlabFilters.dispersion(Z_fosL, Mz_fosL, NX);
            Mz_fosL_H = matlabFilters.mathExpectation(Z_fosL_H, NX);
            Dz_fosL_H = matlabFilters.dispersion(Z_fosL_H, Mz_fosL_H, NX);
            Dz_fosL_H_tilde = matlabFilters.dispersion(Z_fosL_H, Mz_fosL_H, NX);
            Mz_fosG = matlabFilters.mathExpectation(Z_fosG, NX);
            Dz_fosG = matlabFilters.dispersion(Z_fosG, Mz_fosG, NX);
            

            Sx = ones(NX, 1);
            Sx_R = ones(NX, 1);
            Sx_H = ones(NX, 1);
            Seps_aofL = ones(NX, 1);
            Crit_aofL = ones(NX, 1);
            Seps_aofL_H = ones(NX, 1);
            Crit_aofL_H = ones(NX, 1);
            %ICrit_aofL = ones(NX, 1);
            Seps_fosL = ones(NX, 1);
            Crit_fosL = ones(NX, 1);
            Seps_fosL_H = ones(NX, 1);
            Crit_fosL_H = ones(NX, 1);
           % ICrit_fosL = ones(NX, 1);
           Seps_aofG = ones(NX, 1);
            Crit_aofG = ones(NX, 1);
            Seps_fosG = ones(NX, 1);
            Crit_fosG = ones(NX, 1);
            for i = 1:NX
                Sx(i) = sqrt(Dx(i, i));
                Sx_R(i) = Sx(i);
                Sx_H(i) = Sx(i);
                Seps_aofL(i) = sqrt(Deps_aofL(i, i));
                Crit_aofL(i) = Meps_aofL(i)*Meps_aofL(i) + Deps_aofL(i, i);
                Seps_fosL(i) = sqrt(Deps_fosL(i, i));
                Crit_fosL(i) = Meps_fosL(i)*Meps_fosL(i) + Deps_fosL(i, i);
                Seps_aofL_H(i) = sqrt(Deps_aofL_H(i, i));
                Crit_aofL_H(i) = Meps_aofL_H(i)*Meps_aofL_H(i) + Deps_aofL_H(i, i);
                Seps_fosL_H(i) = sqrt(Deps_fosL_H(i, i));
                Crit_fosL_H(i) = Meps_fosL_H(i)*Meps_fosL_H(i) + Deps_fosL_H(i, i);
                Seps_aofG(i) = sqrt(Deps_aofG(i, i));
                Crit_aofG(i) = Meps_aofG(i)*Meps_aofG(i) + Deps_aofG(i, i);
                Seps_fosG(i) = sqrt(Deps_fosG(i, i));
                Crit_fosG(i) = Meps_fosG(i)*Meps_fosG(i) + Deps_fosG(i, i);
            end
            ICrit_aofL = Crit_aofL;
            ICrit_fosL = Crit_fosL;
            ICrit_aofL_H = Crit_aofL_H;
            ICrit_fosL_H = Crit_fosL_H;
            ICrit_aofG = Crit_aofG;
            ICrit_fosG = Crit_fosG;

            MxTable(:, :, 1) = Mx;
            SxTable(:, :, 1) = Sx;
            Sx_RTable(:, :, 1) = Sx;
            Sx_HTable(:, :, 1) = Sx;
            Meps_aofLTable(:, :, 1) = Meps_aofL;
            Seps_aofLTable(:, :, 1) = Seps_aofL;              
            Crit_aofLTable(:, :, 1) = Crit_aofL;
            ICrit_aofLTable(:, :, 1) = ICrit_aofL;   
            Meps_fosLTable(:, :, 1) = Meps_fosL;
            Seps_fosLTable(:, :, 1) = Seps_fosL;              
            Crit_fosLTable(:, :, 1) = Crit_fosL;
            ICrit_fosLTable(:, :, 1) = ICrit_fosL;
            Meps_aofL_HTable(:, :, 1) = Meps_aofL_H;
            Seps_aofL_HTable(:, :, 1) = Seps_aofL_H;              
            Crit_aofL_HTable(:, :, 1) = Crit_aofL_H;
            ICrit_aofL_HTable(:, :, 1) = ICrit_aofL_H;   
            Meps_fosL_HTable(:, :, 1) = Meps_fosL_H;
            Seps_fosL_HTable(:, :, 1) = Seps_fosL_H;              
            Crit_fosL_HTable(:, :, 1) = Crit_fosL_H;
            ICrit_fosL_HTable(:, :, 1) = ICrit_fosL_H;
            Meps_aofGTable(:, :, 1) = Meps_aofG;
            Seps_aofGTable(:, :, 1) = Seps_aofG;              
            Crit_aofGTable(:, :, 1) = Crit_aofG;
            ICrit_aofGTable(:, :, 1) = ICrit_aofG; 
            Meps_fosGTable(:, :, 1) = Meps_fosG;
            Seps_fosGTable(:, :, 1) = Seps_fosG;              
            Crit_fosGTable(:, :, 1) = Crit_fosG;
            ICrit_fosGTable(:, :, 1) = ICrit_fosG;
            tTable(:, 1) = 0;
                      
            sDeltaTime = sqrt(deltaTime);
             
            
            stepForTrajectory = trajectoryNumber; %from N
            XForTrajectoryPlot(:, 1) = X(:, stepForTrajectory);
            X_RForTrajectoryPlot(:, 1) = X_R(:, stepForTrajectory);
            X_HForTrajectoryPlot(:, 1) = X_H(:, stepForTrajectory);
            Z_aofLForTrajectoryPlot(:, 1) = Z_aofL(:, stepForTrajectory);
            Z_fosLForTrajectoryPlot(:, 1) = Z_fosL(:, stepForTrajectory);
            Z_aofGForTrajectoryPlot(:, 1) = Z_aofG(:, stepForTrajectory);
            Z_fosGForTrajectoryPlot(:, 1) = Z_fosG(:, stepForTrajectory);
            
            %progress bar
            
            progressBarCurrent = 0;
            progressBarStep = 1 / K;
            progressBarFunction = waitbar(progressBarCurrent, 'Выполняется фильтрация...');
            
            
            % time loop
            for k = 2:K                
                % realization loop
                
                progressBarCurrent = progressBarCurrent + progressBarStep;
                waitbar(progressBarCurrent, progressBarFunction, 'Выполняется фильтрация...');
                
                deltaV = normrnd(0, sDeltaTime, NW, 2, N);
                q = zeros(1,N);
                tic;

                for i = 1:N
                    
%                     deltaV1 = normrnd(0, sDeltaTime, NX, 1);
%                     deltaV2 = normrnd(0, sDeltaTime, NX, 1); 
                    
                    X(:, i) = X(:, i) + matlabFilters.a(t, X(:, i), condition)*deltaTime + ...
                        matlabFilters.B(t, X(:, i), condition)*deltaV(:, 1, i);

                    deltaY(:, i) = matlabFilters.c(t, X(:, i), condition)*deltaTime + ...
                        matlabFilters.D(t, X(:, i), condition)*deltaV(:, 2, i); 
                    
                    % Хьюн
                    if(condition == 1 && useHun)
                        fi = matlabFilters.a(t, X(:, i), condition) * deltaTime + ...
                            matlabFilters.B(t, X(:, i), condition)*deltaV(:, 1, i);
                        
                        fi_y = matlabFilters.c(t, X(:, i), condition) * deltaTime + ...
                            matlabFilters.D(t, X(:, i), condition)*deltaV(:, 2, i);
                        
                        t_tilde = t + deltaTime;
                        
                        X_tilde(:, i) = X(:, i) + fi;
                                                
                        psi = (matlabFilters.a(t_tilde, X_tilde(:, i), condition) - ...
                            1/2*matlabFilters.Bx(t_tilde, X_tilde(:, i), condition) * matlabFilters.B(t_tilde, X_tilde(:, i), condition)) * deltaTime + ...
                            matlabFilters.B(t_tilde, X_tilde(:, i), condition)*deltaV(:, 1, i);
                        
                        psi_y = (matlabFilters.c(t_tilde, X_tilde(:, i), condition) - ...
                            1/2*matlabFilters.Dx(t_tilde, X_tilde(:, i), condition) * matlabFilters.B(t_tilde, X_tilde(:, i), condition)) * deltaTime + ...
                            matlabFilters.D(t_tilde, X_tilde(:, i), condition)*deltaV(:, 1, i);
                        
                        X_H(:, i) = X_H(:, i) + 1/2*(fi + psi - 1/2 * matlabFilters.Bx(t, X(:, i), condition) * matlabFilters.B(t, X(:, i), condition)*deltaTime);
                        
                        deltaY_H(:, i) = 1/2*(fi_y + psi_y - 1/2 * matlabFilters.Dx(t, X(:, i), condition) * matlabFilters.B(t, X(:, i), condition)*deltaTime);
                        
                        
                        if(buildLAOF)
                            fi_z = matlabFilters.Zeta(t, Z_aofL_H(:, i), P_L_H(:, :, i), condition)*deltaTime + ...
                                matlabFilters.Xi(t, Z_aofL_H(:, i), P_L_H(:, :, i), condition) * fi_y;
                            
                            Z_aofL_H_tilde(:, i) = Z_aofL_H(:, i) + fi_z;
                            
                            fi_p = matlabFilters.Ksi(t, Z_aofL_H(:, i), P_L_H(:, :, i), condition);
                            
                            P_L_H_tilde(:, :, i) = P_L_H(:, :, i) + fi_p;
                            
                            psi_z = (matlabFilters.Zeta(t_tilde, Z_aofL_H_tilde(:, i), P_L_H_tilde(:, :, i), condition) - ...
                                1/2*matlabFilters.Pi(t_tilde, X_tilde(:, i), Z_aofL_H_tilde(:, i), P_L_H_tilde(:, :, i), condition))*deltaTime + ...
                                matlabFilters.Xi(t_tilde, Z_aofL_H_tilde(:, i), P_L_H_tilde(:, :, i), condition) * ...
                                (deltaY_H(:, i) - 1/2 * matlabFilters.Dx(t_tilde, X_tilde(:, i), condition) * ...
                                matlabFilters.B(t_tilde, X_tilde(:, i), condition)*deltaTime);
                            
                            Z_aofL_H(:, i) = Z_aofL_H(:, i) + 1/2*(fi_z + psi_z - 1/2 * ...
                                matlabFilters.Pi(t, X(:, i), Z_aofL_H(:, i), P_L_H(:, :, i),condition)*deltaTime);
                            
                            psi_p = matlabFilters.Ksi(t_tilde, Z_aofL_H_tilde(:, i), P_L_H_tilde(:, :, i), condition);
                            
                            P_L_H(:, :, i) = P_L_H(:, :, i) + 1/2*(fi_p + psi_p);
                            
                        end
                        
                        if(buildLFOS)
                            J_L(:, :, i) = Dx_H - Dz_fosL_H;
                            J_L_tilde(:, :, i) = Dx_H_tilde - Dz_fosL_H_tilde;
                            
                            fi_z = matlabFilters.Zeta(t, Z_fosL_H(:, i), J_L(:, :, i), condition)*deltaTime + ...
                                matlabFilters.Xi(t, Z_fosL_H(:, i), J_L(:, :, i), condition) * fi_y;
                            
                            Z_fosL_H_tilde(:, i) = Z_fosL_H(:, i) + fi_z;
                            
                            psi_z = (matlabFilters.Zeta(t_tilde, Z_fosL_H_tilde(:, i), J_L_tilde(:, :, i), condition) - ...
                                1/2*matlabFilters.Pi(t_tilde, X_tilde(:, i), Z_fosL_H_tilde(:, i), J_L_tilde(:, :, i), condition))*deltaTime + ...
                                matlabFilters.Xi(t_tilde, Z_fosL_H_tilde(:, i), J_L_tilde(:, :, i), condition) * ...
                                (deltaY_H(:, i) - 1/2 * matlabFilters.Dx(t_tilde, X_tilde(:, i), condition) * ...
                                matlabFilters.B(t_tilde, X_tilde(:, i), condition)*deltaTime);
                            
                            Z_fosL_H(:, i) = Z_fosL_H(:, i) + 1/2*(fi_z + psi_z - 1/2 * ...
                                matlabFilters.Pi(t, X(:, i), Z_fosL_H(:, i), J_L(:, :, i),condition)*deltaTime);
                        end
                        
                    end
                    
                    if(useEuler)
                        if(buildLAOF)
                            Z_aofL(:, i) = Z_aofL(:, i) + matlabFilters.a(t, Z_aofL(:, i), condition)*deltaTime + ...
                                matlabFilters.K(t, Z_aofL(:, i), P_L(:, :, i), condition)*(deltaY(:, i) - ...
                                matlabFilters.c(t, Z_aofL(:, i), condition)*deltaTime);
                            P_L(:, :, i) = P_L(:, :, i) +  matlabFilters.Ksi(t, Z_aofL(:, i), P_L(:, :, i), condition)*deltaTime;
                        end


                        if(buildLFOS)
                            Z_fosL(:, i) = Z_fosL(:, i) + matlabFilters.a(t, Z_fosL(:, i), condition)*deltaTime + ...
                                matlabFilters.K(t, Z_fosL(:, i), J_L(:, :, i), condition)*(deltaY(:, i) - ...
                                matlabFilters.c(t, Z_fosL(:, i), condition)*deltaTime);

                             J_L(:, :, i) = Dx - Dz_fosL;
                        end

                       if(buildGAOF)
                            Z_aofG(:, i) = Z_aofG(:, i) + matlabFilters.a_hut(t, Z_aofG(:, i), P_G(:, :, i), condition)*deltaTime + ...
                                matlabFilters.K_hut(t, Z_aofG(:, i), P_G(:, :, i), condition)*(deltaY(:, i) - ...
                                matlabFilters.c_hut(t, Z_aofG(:, i), P_G(:, :, i), condition)*deltaTime);

                            P_G(:, :, i) = P_G(:, :, i) +  matlabFilters.Ksi_hut(t, Z_aofG(:, i), P_G(:, :, i), condition)*deltaTime + ...
                                matlabFilters.R(t, Z_aofG(:, i), condition)* (deltaY(:, i) - matlabFilters.c_hut(t, Z_aofG(:, i), P_G(:, :, i), condition) * deltaTime) ...
                                * matlabFilters.Theta(t, Z_aofG(:, i), P_G(:, :, i), condition);
                       end

                       if(buildGFOS)
                            Z_fosG(:, i) = Z_fosG(:, i) + matlabFilters.a_hut(t, Z_fosG(:, i), J_G(:, :, i), condition)*deltaTime + ...
                                matlabFilters.K_hut(t, Z_fosG(:, i), J_G(:, :, i), condition)*(deltaY(:, i) - ...
                                matlabFilters.c_hut(t, Z_fosG(:, i), J_G(:, :, i), condition)*deltaTime);

                             J_G(:, :, i) = Dx - Dz_fosG;
                       end
                    end
                         
                end
                

                q(5) = toc;
                
                % end of realization loop

                % calculation of characteristics
                Mx = matlabFilters.mathExpectation(X, NX);
                Dx = matlabFilters.dispersion(X, Mx, NX);
                
                Mx_R = matlabFilters.mathExpectation(X_R, NX);
                Dx_R = matlabFilters.dispersion(X_R, Mx_R, NX);
                
                Mx_H = matlabFilters.mathExpectation(X_H, NX);
                Dx_H = matlabFilters.dispersion(X_H, Mx_H, NX);
                
                Mx_H_tilde = matlabFilters.mathExpectation(X_tilde, NX);
                Dx_H_tilde = matlabFilters.dispersion(X_tilde, Mx_H_tilde, NX);
    
                if(buildLAOF)
                    Meps_aofL = matlabFilters.mathExpectation(X - Z_aofL, NX);
                    Deps_aofL = matlabFilters.dispersion(X - Z_aofL, Meps_aofL, NX);
                    Meps_aofL_H = matlabFilters.mathExpectation(X - Z_aofL_H, NX);
                    Deps_aofL_H = matlabFilters.dispersion(X - Z_aofL_H, Meps_aofL_H, NX);
                end
                if(buildLFOS)
                    Meps_fosL = matlabFilters.mathExpectation(X - Z_fosL, NX);
                    Deps_fosL = matlabFilters.dispersion(X - Z_fosL, Meps_fosL, NX);
                    Mz_fosL = matlabFilters.mathExpectation(Z_fosL, NX);
                    Dz_fosL = matlabFilters.dispersion(Z_fosL, Mz_fosL, NX);
                    Meps_fosL_H = matlabFilters.mathExpectation(X - Z_fosL_H, NX);
                    Deps_fosL_H = matlabFilters.dispersion(X - Z_fosL_H, Meps_fosL_H, NX);
                    Mz_fosL_H = matlabFilters.mathExpectation(Z_fosL_H, NX);
                    Dz_fosL_H = matlabFilters.dispersion(Z_fosL_H, Mz_fosL_H, NX);
                    Mz_fosL_H_tilde = matlabFilters.mathExpectation(Z_fosL_H_tilde, NX);
                    Dz_fosL_H_tilde = matlabFilters.dispersion(Z_fosL_H_tilde, Mz_fosL_H_tilde, NX);
                end
                if(buildGAOF)
                    Meps_aofG = matlabFilters.mathExpectation(X - Z_aofG, NX);
                    Deps_aofG = matlabFilters.dispersion(X - Z_aofG, Meps_aofG, NX);
                end
                if(buildGFOS)
                    Meps_fosG = matlabFilters.mathExpectation(X - Z_fosG, NX);
                    Deps_fosG = matlabFilters.dispersion(X - Z_fosG, Meps_fosG, NX);
                    Mz_fosG = matlabFilters.mathExpectation(Z_fosG, NX);
                    Dz_fosG = matlabFilters.dispersion(Z_fosG, Mz_fosG, NX);
                end
                
                
                Sx = ones(NX, 1);
                Sx_R = ones(NX, 1);
                Sx_H = ones(NX, 1);
                Seps_aofL = ones(NX, 1);
                Crit_aofL = ones(NX, 1);
                ICrit_aofL = ones(NX, 1);
                Seps_fosL = ones(NX, 1);
                Crit_fosL = ones(NX, 1);
                ICrit_fosL = ones(NX, 1);
                Seps_aofG = ones(NX, 1);
                Crit_aofG = ones(NX, 1);
                ICrit_aofG = ones(NX, 1);
                Seps_fosG = ones(NX, 1);
                Crit_fosG = ones(NX, 1);
                ICrit_fosG = ones(NX, 1);
                for i = 1:NX
                    Sx(i) = sqrt(Dx(i, i));
                    Sx_R(i) = sqrt(Dx_R(i, i));
                    Sx_H(i) = sqrt(Dx_H(i, i));
                    if(buildLAOF)
                        Seps_aofL(i) = sqrt(Deps_aofL(i, i));
                        Crit_aofL(i) = Meps_aofL(i)*Meps_aofL(i) + Deps_aofL(i,i);
                        Seps_aofL_H(i) = sqrt(Deps_aofL_H(i, i));
                        Crit_aofL_H(i) = Meps_aofL_H(i)*Meps_aofL_H(i) + Deps_aofL_H(i,i);
                    end
                    if(buildLFOS)
                        Seps_fosL(i) = sqrt(Deps_fosL(i, i));
                        Crit_fosL(i) = Meps_fosL(i)*Meps_fosL(i) + Deps_fosL(i,i);
                        Seps_fosL_H(i) = sqrt(Deps_fosL_H(i, i));
                        Crit_fosL_H(i) = Meps_fosL_H(i)*Meps_fosL_H(i) + Deps_fosL_H(i,i);
                    end
                    if(buildGAOF)
                        Seps_aofG(i) = sqrt(Deps_aofG(i, i));
                        Crit_aofG(i) = Meps_aofG(i)*Meps_aofG(i) + Deps_aofG(i,i);
                    end
                    if(buildGFOS)
                        Seps_fosG(i) = sqrt(Deps_fosG(i, i));
                        Crit_fosG(i) = Meps_fosG(i)*Meps_fosG(i) + Deps_fosG(i,i);
                    end
                end   
                if(buildLAOF)
                    ICrit_aofL = ICrit_aofLTable(:, :, k-1) + (Crit_aofL - ICrit_aofLTable(:, :, k-1)) / k;
                    ICrit_aofL_H = ICrit_aofL_HTable(:, :, k-1) + (Crit_aofL_H - ICrit_aofL_HTable(:, :, k-1)) / k;
                end
                if(buildLFOS)
                    ICrit_fosL = ICrit_fosLTable(:, :, k-1) + (Crit_fosL - ICrit_fosLTable(:, :, k-1)) / k;
                    ICrit_fosL_H = ICrit_fosL_HTable(:, :, k-1) + (Crit_fosL_H - ICrit_fosL_HTable(:, :, k-1)) / k;
                end
                if(buildGAOF)
                    ICrit_aofG = ICrit_aofGTable(:, :, k-1) + (Crit_aofG - ICrit_aofGTable(:, :, k-1)) / k;
                end
                if(buildGFOS)
                    ICrit_fosG = ICrit_fosGTable(:, :, k-1) + (Crit_fosG - ICrit_fosGTable(:, :, k-1)) / k;
                end
                
                MxTable(:, :, k) = Mx;
                SxTable(:, :, k) = Sx;
                Sx_RTable(:, :, k) = Sx_R;
                Sx_HTable(:, :, k) = Sx_H;
                Meps_aofLTable(:, :, k) = Meps_aofL;
                Seps_aofLTable(:, :, k) = Seps_aofL;                
                Crit_aofLTable(:, :, k) = Crit_aofL;
                ICrit_aofLTable(:, :, k) = ICrit_aofL;
                Meps_fosLTable(:, :, k) = Meps_fosL;
                Seps_fosLTable(:, :, k) = Seps_fosL;                
                Crit_fosLTable(:, :, k) = Crit_fosL;
                ICrit_fosLTable(:, :, k) = ICrit_fosL;
                Meps_aofL_HTable(:, :, k) = Meps_aofL_H;
                Seps_aofL_HTable(:, :, k) = Seps_aofL_H;                
                Crit_aofL_HTable(:, :, k) = Crit_aofL_H;
                ICrit_aofL_HTable(:, :, k) = ICrit_aofL_H;
                Meps_fosL_HTable(:, :, k) = Meps_fosL_H;
                Seps_fosL_HTable(:, :, k) = Seps_fosL_H;                
                Crit_fosL_HTable(:, :, k) = Crit_fosL_H;
                ICrit_fosL_HTable(:, :, k) = ICrit_fosL_H;   
                Meps_aofGTable(:, :, k) = Meps_aofG;
                Seps_aofGTable(:, :, k) = Seps_aofG;                
                Crit_aofGTable(:, :, k) = Crit_aofG;
                ICrit_aofGTable(:, :, k) = ICrit_aofG;
                Meps_fosGTable(:, :, k) = Meps_fosG;
                Seps_fosGTable(:, :, k) = Seps_fosG;                
                Crit_fosGTable(:, :, k) = Crit_fosG;
                ICrit_fosGTable(:, :, k) = ICrit_fosG;
                tTable(:, k) = t;
                
                t = (k) * deltaTime;
                
                XForTrajectoryPlot(:, k) = X(:, stepForTrajectory);
                X_RForTrajectoryPlot(:, k) = X_R(:, stepForTrajectory);
                X_HForTrajectoryPlot(:, k) = X_H(:, stepForTrajectory);
                if(buildLAOF)
                    Z_aofLForTrajectoryPlot(:, k) = Z_aofL(:, stepForTrajectory); 
                end
                if(buildLFOS)
                    Z_fosLForTrajectoryPlot(:, k) = Z_fosL(:, stepForTrajectory);       
                end
                if(buildGAOF)
                    Z_aofGForTrajectoryPlot(:, k) = Z_aofG(:, stepForTrajectory); 
                end
                if(buildGFOS)
                    Z_fosGForTrajectoryPlot(:, k) = Z_fosG(:, stepForTrajectory);       
                end
            end
            % end of time loop
            
            delete(progressBarFunction);

            % displaying the table
%             if(condition == 1)
%                 disp('     time       Mx     Meps_aofL    Sx    Seps_aofL   Crit_aofL  ICrit_aofL');
%                 answMatrix = [tTable, MxTable, Meps_aofLTable, SxTable, Seps_aofLTable, Crit_aofLTable, ICrit_aofLTable];
%                 for i = 1:K
%                     disp(answMatrix(:,:,i));
%                 end
%                 disp('     time       Mx     Meps_fosL    Sx    Seps_fosL   Crit_fosL  ICrit_fosL');
%                 answMatrix = [tTable, MxTable, Meps_fosLTable, SxTable, Seps_fosLTable, Crit_fosLTable, ICrit_fosLTable];
%                 for i = 1:K
%                     disp(answMatrix(:,:,i));
%                 end
%             end
%             if(condition == 2)
%                 if(buildLAOF)
%                     disp('     time       Mx1       Mx2   MepsAofL1  MepsAofL2    Sx1      Sx2    SepsAofL1 SepsAofL2  CritAofL1 CritAofL2 ICritAofL1 ICritAofL2');
%                     answMatrix = [tTable, MxTable(1, :, :), MxTable(2, :, :), Meps_aofLTable(1, :, :), Meps_aofLTable(2, :, :), ...
%                         SxTable(1, :, :), SxTable(2, :, :), Seps_aofLTable(1, :, :), Seps_aofLTable(2, :, :), Crit_aofLTable(1, :, :), Crit_aofLTable(2, :, :), ...
%                         ICrit_aofLTable(1, :, :), ICrit_aofLTable(2, :, :)];
%                     for i = 1:K
%                         disp(answMatrix(:,:,i));
%                     end
%                 end
%                 
%                 if(buildLFOS)
%                     disp('     time       Mx1       Mx2   MepsFosL1  MepsFosL2    Sx1      Sx2    SepsFosL1 SepsFosL2  CritFosL1 CritFosL2 ICritFosL1 ICritFosL2');
%                     answMatrix = [tTable, MxTable(1, :, :), MxTable(2, :, :), Meps_fosLTable(1, :, :), Meps_fosLTable(2, :, :), ...
%                         SxTable(1, :, :), SxTable(2, :, :), Seps_fosLTable(1, :, :), Seps_fosLTable(2, :, :), Crit_fosLTable(1, :, :), Crit_fosLTable(2, :, :), ...
%                         ICrit_fosLTable(1, :, :), ICrit_fosLTable(2, :, :)];
%                     for i = 1:K
%                         disp(answMatrix(:,:,i));
%                     end
%                 end
%                 
%                 
%                 if(buildGAOF)
%                     disp('     time       Mx1       Mx2   MepsAofG1  MepsAofG2    Sx1      Sx2    SepsAofG1 SepsAofG2  CritAofG1 CritAofG2 ICritAofG1 ICritAofG2');
%                     answMatrix = [tTable, MxTable(1, :, :), MxTable(2, :, :), Meps_aofGTable(1, :, :), Meps_aofGTable(2, :, :), ...
%                         SxTable(1, :, :), SxTable(2, :, :), Seps_aofGTable(1, :, :), Seps_aofGTable(2, :, :), Crit_aofGTable(1, :, :), Crit_aofGTable(2, :, :), ...
%                         ICrit_aofGTable(1, :, :), ICrit_aofGTable(2, :, :)];
%                     for i = 1:K
%                         disp(answMatrix(:,:,i));
%                     end
%                 end
%                 
%                 
%                 if(buildGFOS)
%                     disp('     time       Mx1       Mx2   MepsFosG1  MepsFosG2    Sx1      Sx2    SepsFosG1 SepsFosG2  CritFosG1 CritFosG2 ICritFosG1 ICritFosG2');
%                     answMatrix = [tTable, MxTable(1, :, :), MxTable(2, :, :), Meps_fosGTable(1, :, :), Meps_fosGTable(2, :, :), ...
%                         SxTable(1, :, :), SxTable(2, :, :), Seps_fosGTable(1, :, :), Seps_fosGTable(2, :, :), Crit_fosGTable(1, :, :), Crit_fosGTable(2, :, :), ...
%                         ICrit_fosGTable(1, :, :), ICrit_fosGTable(2, :, :)];
%                     for i = 1:K
%                         disp(answMatrix(:,:,i));
%                     end
%                 end
%             end
            if(condition == 4)
%                 disp('     time       Mx     Meps_aofL    Sx    Seps_aofL   Crit_aofL  ICrit_aofL');
%                 answMatrix = [tTable, MxTable, Meps_aofLTable, SxTable, Seps_aofLTable, Crit_aofLTable, ICrit_aofLTable];
%                 for i = 1:K
%                     disp(answMatrix(:,:,i));
%                 end
%                 disp('     time       Mx     Meps_fosL    Sx    Seps_fosL   Crit_fosL  ICrit_fosL');
%                 answMatrix = [tTable, MxTable, Meps_fosLTable, SxTable, Seps_fosLTable, Crit_fosLTable, ICrit_fosLTable];
%                 for i = 1:K
%                     disp(answMatrix(:,:,i));
%                 end
%                 disp('     time       Mx     Meps_aofG    Sx    Seps_aofG   Crit_aofG  ICrit_aofG');
%                 answMatrix = [tTable, MxTable, Meps_aofGTable, SxTable, Seps_aofGTable, Crit_aofGTable, ICrit_aofGTable];
%                 for i = 1:K
%                     disp(answMatrix(:,:,i));
%                 end
%                 disp('     time       Mx     Meps_fosG    Sx    Seps_fosG   Crit_fosG  ICrit_fosG');
%                 answMatrix = [tTable, MxTable, Meps_fosGTable, SxTable, Seps_fosGTable, Crit_fosGTable, ICrit_fosGTable];
%                 for i = 1:K
%                     disp(answMatrix(:,:,i));
%                 end
                disp('     time       Mx    Meps_aofL Meps_fosL Meps_aofG Meps_fosG    ');
                answMatrix = [tTable, MxTable, Meps_aofLTable, Meps_fosLTable, Meps_aofGTable, Meps_fosGTable];
                for i = 1:K
                    disp(answMatrix(:,:,i));
                end
                
                disp('     time       Sx    Seps_aofL Seps_fosL Seps_aofG Seps_fosG  ');
                answMatrix = [tTable, SxTable, Seps_aofLTable, Seps_fosLTable, Seps_aofGTable, Seps_fosGTable];
                for i = 1:K
                    disp(answMatrix(:,:,i));
                end
            end


            
%             plot Sx, Seps, Crit
            SxToPlot = zeros(1, NX);
            Sx_RToPlot = zeros(1, NX);
            Sx_HToPlot = zeros(1, NX);
            SeAofLToPlot = zeros(1, NX);
            SeFosLToPlot = zeros(1, NX);
            SeAofL_HToPlot = zeros(1, NX);
            SeFosL_HToPlot = zeros(1, NX);
            SeAofGToPlot = zeros(1, NX);
            SeFosGToPlot = zeros(1, NX);
%             sqrtCrit_aofLToPlot = zeros(1, NX);
%             sqrtCrit_fosLToPlot = zeros(1, NX);
%             ICrit_aofLToPlot = zeros(1, NX);
%             ICrit_fosLToPlot = zeros(1, NX);
            MxToPlot = zeros(1, NX);
            MeAofLToPlot = zeros(1, NX);
            MeFosLToPlot = zeros(1, NX);
            MeAofGToPlot = zeros(1, NX);
            MeFosGToPlot = zeros(1, NX);
            tToPlot = zeros(1, NX);
            
            
            for j = 1:NX
                for i = 1:K
                    SxToPlot(i) = SxTable(j, 1, i);
                    Sx_RToPlot(i) = Sx_RTable(j, 1, i);
                    Sx_HToPlot(i) = Sx_HTable(j, 1, i);
                    SeAofLToPlot(i) = Seps_aofLTable(j, 1, i);
                    SeFosLToPlot(i) = Seps_fosLTable(j, 1, i);
                    SeAofL_HToPlot(i) = Seps_aofL_HTable(j, 1, i);
                    SeFosL_HToPlot(i) = Seps_fosL_HTable(j, 1, i);
                    SeAofGToPlot(i) = Seps_aofGTable(j, 1, i);
                    SeFosGToPlot(i) = Seps_fosGTable(j, 1, i);
%                     sqrtCrit_aofLToPlot(i) = sqrt(Crit_aofLTable(j, 1, i));
%                     sqrtCrit_fosLToPlot(i) = sqrt(Crit_fosLTable(j, 1, i));
%                     ICrit_aofLToPlot(i) = ICrit_aofLTable(j, 1, i);
%                     ICrit_fosLToPlot(i) = ICrit_fosLTable(j, 1, i);

                    tToPlot(i) = tTable(1, 1, i);
                end               
                if(condition == 3) 
                    if(j == 2)
%                         SxToPlot = SxToPlot*180/pi;
%                         SeAofLToPlot = SeAofLToPlot*180/pi;
%                         SeFosLToPlot = SeFosLToPlot*180/pi;
                    else
                        SxToPlot = SxToPlot;
                        SeAofLToPlot = SeAofLToPlot;
                        SeFosLToPlot = SeFosLToPlot;
                    end
                end

%                 subplot(NX,1,j);
%                 plot(tToPlot.', SxToPlot.');
%                 hold on;
%                 if(buildLAOF)
%                     plot(tToPlot.', SeAofLToPlot.', 'o-', 'color', [240, 110, 50]/255);
%                     hold on;
%                 end
% 
%                 if(buildLFOS)
%                     plot(tToPlot.', SeFosLToPlot.', 'kx--');
%                     hold on;
%                 end
%                 
%                 if(buildGAOF)
%                     plot(tToPlot.', SeAofGToPlot.', '*-', 'color', [200, 180, 10]/255);
%                     hold on;
%                 end
%                 
%                 if(buildGFOS)
%                     plot(tToPlot.', SeFosGToPlot.', 's--', 'color', [0, 175, 250]/255);
%                     hold on;
%                 end

                if(useEuler)
                    figure('name', 'Sx, Se');
                    subplot(1,NX,j);
                    plot(tToPlot.', SxToPlot.');
                    hold on;
                    if(buildLAOF)
                        plot(tToPlot.', SeAofLToPlot.', 'o-', 'color', [240, 110, 50]/255);
                        hold on;
                    end

                    if(buildLFOS)
                        plot(tToPlot.', SeFosLToPlot.', 'kx--');
                        hold on;
                    end

                    if(buildGAOF)
                        plot(tToPlot.', SeAofGToPlot.', '*-', 'color', [200, 180, 10]/255);
                        hold on;
                    end

                    if(buildGFOS)
                        plot(tToPlot.', SeFosGToPlot.', 's--', 'color', [0, 175, 250]/255);
                        hold on;
                    end

    %                 plot(tToPlot.', sqrtCrit_aofLToPlot.', 'color', [0, 110, 50]/255);
    %                 plot(tToPlot.', sqrtCrit_fosLToPlot.', 'color', [200, 180, 10]/255);

                    % legend format
                    if(buildLAOF && buildLFOS && buildGAOF && buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'), strcat('Se', int2str(j), ' FOS-L'), strcat('Se', int2str(j), ' AOF-G'), strcat('Se', int2str(j), ' FOS-G'));
                    end
                    if(buildLAOF && buildLFOS && buildGAOF && ~buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'), strcat('Se', int2str(j), ' FOS-L'), strcat('Se', int2str(j), ' AOF-G'));
                    end
                    if(buildLAOF && buildLFOS && ~buildGAOF && buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'), strcat('Se', int2str(j), ' FOS-L'), strcat('Se', int2str(j), ' FOS-G'));
                    end
                    if(buildLAOF && buildLFOS && ~buildGAOF && ~buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'), strcat('Se', int2str(j), ' FOS-L'));
                    end
                    if(buildLAOF && ~buildLFOS && buildGAOF && buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'), strcat('Se', int2str(j), ' AOF-G'), strcat('Se', int2str(j), ' FOS-G'));
                    end
                    if(buildLAOF && ~buildLFOS && buildGAOF && ~buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'), strcat('Se', int2str(j), ' AOF-G'));
                    end
                    if(buildLAOF && ~buildLFOS && ~buildGAOF && buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'), strcat('Se', int2str(j), ' FOS-G'));
                    end
                    if(buildLAOF && ~buildLFOS && ~buildGAOF && ~buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'));
                    end
                    if(~buildLAOF && buildLFOS && buildGAOF && buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' FOS-L'), strcat('Se', int2str(j), ' AOF-G'), strcat('Se', int2str(j), ' FOS-G'));
                    end
                    if(~buildLAOF && buildLFOS && buildGAOF && ~buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' FOS-L'), strcat('Se', int2str(j), ' AOF-G'));
                    end
                    if(~buildLAOF && buildLFOS && ~buildGAOF && buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' FOS-L'), strcat('Se', int2str(j), ' FOS-G'));
                    end
                    if(~buildLAOF && buildLFOS && ~buildGAOF && ~buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' FOS-L'));
                    end
                    if(~buildLAOF && ~buildLFOS && buildGAOF && buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-G'), strcat('Se', int2str(j), ' FOS-G'));
                    end
                    if(~buildLAOF && ~buildLFOS && buildGAOF && ~buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-G'));
                    end
                    if(~buildLAOF && ~buildLFOS && ~buildGAOF && buildGFOS)
                        legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' FOS-G'));
                    end
                    if(~buildLAOF && ~buildLFOS && ~buildGAOF && ~buildGFOS)
                        legend(strcat('Sx', int2str(j)));
                    end

                    if(condition == 1)
                        str = 'Полиномиальный пример: ';
                    end
                    if(condition == 2)
                        str = 'Осциллятор Ван-дер-Поля: ';
                    end
                    if(condition == 3)
                        str = 'Спуск: ';
                    end
                    if(condition == 4)
                        str = 'Тригонометрический пример: ';
                    end
                    title(['\fontsize{20}', str, ' N=', int2str(N), ', deltaTime=', num2str(deltaTime), '. ']);
                    xlabel('t'); 
                end
                
                
                if(useHun)
                    figure('name', 'Sx_H, Se_H');
                    subplot(1,NX,j);
                    plot(tToPlot.', Sx_HToPlot.');
                    hold on;
                    if(buildLAOF)
                        plot(tToPlot.', SeAofL_HToPlot.', 'o-', 'color', [240, 110, 50]/255);
                        hold on;
                    end

                    if(buildLFOS)
                        plot(tToPlot.', SeFosL_HToPlot.', 'kx--');
                        hold on;
                    end

                    % legend format
                    if(buildLAOF && buildLFOS)
                        legend('Sx_H', 'Se_H AOF-L', 'Se_H FOS-L');
                    end
                    if(buildLAOF && ~buildLFOS)
                        legend('Sx_H', 'Se_H AOF-L');
                    end
                    if(~buildLAOF && buildLFOS)
                        legend('Sx_H', 'Se_H FOS-L');
                    end
                    if(~buildLAOF && ~buildLFOS)
                        legend('Sx_H');
                    end

                    str = 'Полиномиальный пример: ';

                    title(['\fontsize{20}', str, ' N=', int2str(N), ', deltaTime=', num2str(deltaTime), '. ']);
                    xlabel('t'); 
                end
                
            end     
               
            
            
            if(useEuler && useHun)
                figure('name', 'Sx-Sx_H, Se-Se_H');

                plot(tToPlot.', SxToPlot.', '*-');
                hold on;
                plot(tToPlot.', Sx_HToPlot.', 's--');
                hold on;
                if(buildLAOF)
                    plot(tToPlot.', SeAofLToPlot.');
                    hold on;
                    plot(tToPlot.', SeAofL_HToPlot.');
                    hold on;
                end
                if(buildLFOS)
                    plot(tToPlot.', SeFosLToPlot.', 'x-');
                    hold on;
                    plot(tToPlot.', SeFosL_HToPlot.', 'o--');
                end
                if(buildLAOF && buildLFOS)
                    legend('Sx', 'Sx_H', 'Se AOF-L', 'Se_H AOF-L', 'Se FOS-L', 'Se_H FOS-L');
                end
                if(~buildLAOF && buildLFOS)
                    legend('Sx', 'Sx_H', 'Se FOS-L', 'Se_H FOS-L');
                end
                if(buildLAOF && ~buildLFOS)
                    legend('Sx', 'Sx_H', 'Se AOF-L', 'Se_H AOF-L');
                end
                if(~buildLAOF && ~buildLFOS)
                    legend('Sx', 'Sx_H');
                end
                 str = 'Полиномиальный пример: ';

                 title(['\fontsize{20}', str, ' N=', int2str(N), ', deltaTime=', num2str(deltaTime), '. ']);
%                     
            end
            
            
            
            
            
            % plot trajectory 
            if(buildTrajectory)
                
                if(buildLAOF)
                    figure('name', 'trajectory: X(t), Z(t) AOF-L and confidence interval');
                    for j = 1:NX
                        subplot(NX,1,j);
                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).');
                        hold on;
                        plot(tToPlot.'.', Z_aofLForTrajectoryPlot(j, :).', 'o-', 'color', [240, 110, 50]/255);
                        hold on;
                        xlabel('t');

                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).' + 3*SeAofLToPlot.', '--', 'color', [150, 150, 150]/255);
                        hold on;
                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).' - 3*SeAofLToPlot.', '--', 'color', [150, 150, 150]/255);
                        hold on;

                        legend(strcat('X', int2str(j)), strcat('Z', int2str(j), ' AOF-L'), 'X + 3Se AOF-L', 'X - 3Se AOF-L');
                        title(['\fontsize{20}', 'Траектория j= ', int2str(stepForTrajectory), ' из ', int2str(N)]);
                    end
                end
                
                if(buildLFOS)
                    figure('name', 'trajectory: X(t), Z(t) FOS-L and confidence interval');
                    for j = 1:NX             
                        subplot(NX,1,j);
                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).');
                        hold on;
                        plot(tToPlot.'.', Z_fosLForTrajectoryPlot(j, :).', 'kx--');
                        xlabel('t');
                        hold on;

                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).' + 3*SeFosLToPlot.', '--', 'color', [150, 150, 150]/255);
                        hold on;
                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).' - 3*SeFosLToPlot.', '--', 'color', [150, 150, 150]/255);
                        hold on;

                        legend(strcat('X', int2str(j)), strcat('Z', int2str(j), ' FOS-L'), 'X + 3Se FOS-L', 'X - 3Se FOS-L');
                        title(['\fontsize{20}', 'Траектория j= ', int2str(stepForTrajectory), ' из ', int2str(N)]);
                    end
                end
                
                if(buildGAOF)
                    figure('name', 'trajectory: X(t), Z(t) AOF-G and confidence interval');
                    for j = 1:NX
                        subplot(NX,1,j);
                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).');
                        hold on;
                        plot(tToPlot.'.', Z_aofGForTrajectoryPlot(j, :).', '*-', 'color', [200, 180, 10]/255);
                        hold on;
                        xlabel('t');

                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).' + 3*SeAofGToPlot.', '--', 'color', [150, 150, 150]/255);
                        hold on;
                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).' - 3*SeAofGToPlot.', '--', 'color', [150, 150, 150]/255);
                        hold on;

                        legend(strcat('X', int2str(j)), strcat('Z', int2str(j), ' AOF-G'), 'X + 3Se AOF-G', 'X - 3Se AOF-G');
                        title(['\fontsize{20}', 'Траектория j= ', int2str(stepForTrajectory), ' из ', int2str(N)]);
                    end
                end
                
                if(buildGFOS)
                    figure('name', 'trajectory: X(t), Z(t) FOS-G and confidence interval');
                    for j = 1:NX             
                        subplot(NX,1,j);
                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).');
                        hold on;
                        plot(tToPlot.'.', Z_fosGForTrajectoryPlot(j, :).', 's--', 'color', [0, 175, 250]/255);
                        xlabel('t');
                        hold on;

                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).' + 3*SeFosGToPlot.', '--', 'color', [150, 150, 150]/255);
                        hold on;
                        plot(tToPlot.'.', XForTrajectoryPlot(j, :).' - 3*SeFosGToPlot.', '--', 'color', [150, 150, 150]/255);
                        hold on;

                        legend(strcat('X', int2str(j)), strcat('Z', int2str(j), ' FOS-G'), 'X + 3Se FOS-G', 'X - 3Se FOS-G');
                        title(['\fontsize{20}', 'Траектория j= ', int2str(stepForTrajectory), ' из ', int2str(N)]);
                    end
                end
            end
            
            
            % доделать гаоф и гфос далее
            
            
            %build Mx
            if(buildMx)
                figure('name', 'Mx, Me');
                
                maxMx = zeros(1, NX);
                maxMeAOFL = zeros(1, NX);
                maxMeFOSL = zeros(1, NX);
                for jj = 1:NX
                    for ii = 1:K
                        if(maxMx(1, jj) < abs(MxTable(jj, 1, ii)))
                            maxMx(1, jj) = abs(MxTable(jj, 1, ii));
                        end
                        if(maxMeAOFL(1, jj) < abs(Meps_aofLTable(jj, 1, ii)))
                            maxMeAOFL(1, jj) = abs(Meps_aofLTable(jj, 1, ii));
                        end
                        if(maxMeFOSL(1, jj) < abs(Meps_fosLTable(jj, 1, ii)))
                            maxMeFOSL(1, jj) = abs(Meps_fosLTable(jj, 1, ii));
                        end
                    end
                end
                
                for j = 1:NX
                    for i = 1:K
                        MxToPlot(i) = MxTable(j, 1, i);
                        MeAofLToPlot(i) = Meps_aofLTable(j, 1, i);
                        MeFosLToPlot(i) = Meps_fosLTable(j, 1, i);
                    end  
                    subplot(NX,1,j);
                    plot(tToPlot.', MxToPlot.');
                    hold on;
                    if(buildLAOF)
                        plot(tToPlot.', MeAofLToPlot.', 'o-', 'color', [240, 110, 50]/255);
                        hold on;
                    end

                    if(buildLFOS)
                        plot(tToPlot.', MeFosLToPlot.', 'kx--');
                    end
                    
                    
                  
                    % legend format
                    if(buildLAOF && buildLFOS)
                        legend(strcat('Mx', num2str(j)), strcat('Me', num2str(j), ' AOF-L'), strcat('Me', num2str(j), ' FOS-L'));
                    end
                    if(buildLAOF && ~buildLFOS)
                        legend(strcat('Mx', num2str(j)), strcat('Me', num2str(j), ' AOF-L'));
                    end
                    if(~buildLAOF && buildLFOS)
                        legend(strcat('Mx', num2str(j)), strcat('Me', num2str(j), ' FOS-L'));
                    end
                    if(~buildLAOF && ~buildLFOS)
                        legend(strcat('Mx', num2str(j)));
                    end
                    
                    title(['\fontsize{20}', 'max|Mx|=', num2str(maxMx(1, j), '%10.2e'), '; max|Me AOF-L|=', num2str(maxMeAOFL(1, j), '%10.2e'), '; max|Me FOS-L|=', num2str(maxMeFOSL(1, j), '%10.2e')]);
                end
            end
            
            % Только для скрина 3 примера
            if(buildMx)
                figure('name', 'Mx, Sx');
                for j = 1:NX
                    for i = 1:K
                        MxToPlot(i) = MxTable(j, 1, i);
                        SxToPlot(i) = SxTable(j, 1, i);
                    end  
                    if(condition == 3) 
                        if(j == 2)
                            %SxToPlot = SxToPlot*180/pi;
                        end
                    end
                    subplot(1,NX,j);
                    plot(tToPlot.', MxToPlot.');
                    hold on;
                    plot(tToPlot.', SxToPlot.');
                    legend(strcat('Mx', num2str(j)), strcat('Sx', num2str(j)));
                    if(condition == 1)
                        str = 'Полиномиальный пример: ';
                    end
                    if(condition == 2)
                        str = 'Осциллятор Ван-дер-Поля: ';
                    end
                    if(condition == 3)
                        str = 'Спуск: ';
                    end
                    title(['\fontsize{20}', str, ' N=', int2str(N), ', deltaTime=', num2str(deltaTime), '. ']);
                end
            
                
            end
            

     
        end        
        % end of Kalman filter function
        
        
       
       
       
       
    end
    
end

