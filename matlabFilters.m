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
        
        function [ c ] = c(t, x, condition)
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
        
        function [ B ] = B(t, x, condition)
            global gamma;
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
        
        function [ D ] = D(t, x, condition)
            global gamma;
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
        
        function [ R ] = R(t, x, condition)
            D = matlabFilters.D(t, x, condition);
            R = D * D.';
        end
        
        function [ Q ] = Q(t, x, condition)
            B = matlabFilters.B(t, x, condition);
            Q = B * B.';
        end
        
        
        % Support functions for Kalman filter
        
        function [ A ] = A(t, x, condition)                               % матрица Якоби первых производных a(t, x) по x
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
        
        function [ G ] = G(t, x, condition)                               % матрица Якоби первых производных с(t, x) по x
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
        
        
        % structure functions
        
        function [ K ] = K(t, x, p, condition) 
            K = p * matlabFilters.G(t, x, condition).' /(matlabFilters.R(t, x, condition));
%             p
%             matlabFilters.G(t, x, condition)
%             matlabFilters.G(t, x, condition).'
%             matlabFilters.R(t, x, condition)
        end
        
        function [ Ksi ] = Ksi(t, x, p, condition)                        
            Ksi = matlabFilters.A(t, x, condition) * p + p * matlabFilters.A(t, x, condition).' + ...
                matlabFilters.Q(t, x, condition) - matlabFilters.K(t, x, p, condition) * ... 
                matlabFilters.R(t, x, condition) * matlabFilters.K(t, x, p, condition).';
        end
        
        
        % Support functions
        
        function [ M ] = changeMx(i, x, Mx)
            M = Mx + (x - Mx)/i;
        end
        
        function [ D ] = changeDx(i, x, Mx, Dx)
            D = Dx + ((x - Mx) * (x - Mx) * (i - 1)/i - Dx)/i;
        end
        
        function [ M ] = mathExpectation(Temp, condition)
            global N;          
            M = zeros(condition, 1);
            for i = 1:N
                M = M + Temp(:, i);
            end
            M = M * 1/N;
        end
        
        
        function [ D ] = dispersion(Temp, M, condition)
            global N;
            if(condition == 1)
                D = 0;
                for i = 1:N
                    D = D + (Temp(i) - M)*(Temp(i) - M);
                end
                D = D * 1/(N-1);
            end
            if(condition > 1)
                D = zeros(condition, condition);
                for i = 1:N
                    for j = 1:condition
                        for k = 1:condition
                            D(j, k) = D(j, k) + (Temp(j, i) - M(j))*(Temp(k, i) - M(k));
                        end
                    end
                end
                for j = 1:condition
                    for k = 1:condition
                        D(j, k) = D(j, k)/(N-1);
                    end
                end
            end   
        end
        
        
        
        
        function [ ] = initialConditions(a0cond1G, a1cond1G, a3cond1G, b0cond1G, b1cond1G, c1cond1G, c2cond1G, d0cond1G, d1cond1G, MoXcond1G, SoXcond1G, ...
                                        omegacond2G, alphacond2G, betacond2G, d11cond2G, d22cond2G, SoXcond2G, SoYcond2G, MoXcond2G, MoYcond2G, gammaG)
                        
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
             
            
            %cond 3 - kilometers
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
            sigmaM = 0.01;
            sigmaA = 0.00002;
           
        end
        
        
        function[ ] = main(condition, N1, K1, deltaTime1, buildTrajectory, trajectoryNumber, buildLAOF, buildLFOS, buildMx)
            
            % Initial conditions
            global MoXcond1;
            global MoYcond1;
            global SoXcond1;
            global SoYcond1;
            global SoZcond1;
            
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
            Meps_aofLTable = ones(NX, 1, K);
            Seps_aofLTable = ones(NX, 1, K);
            Crit_aofLTable = ones(NX, 1, K);
            ICrit_aofLTable = ones(NX, 1, K);
            Meps_fosLTable = ones(NX, 1, K);
            Seps_fosLTable = ones(NX, 1, K);
            Crit_fosLTable = ones(NX, 1, K);
            ICrit_fosLTable = ones(NX, 1, K);
            
            t = deltaTime;
            X = ones(NX, N);
            deltaY = ones(NY, N);
            Z_aofL = ones(NX, N);
            Z_fosL = ones(NX, N);
            P = ones(NX, NX, N);
            J = ones(NX, NX, N);
                      
            
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
            if(condition == 3)
                 for i = 1:N;
                    P(:, :, i) = [SoXv*SoXv, 0, 0; 0, SoXteta*SoXteta, 0; 0, 0, SoXh*SoXh];
                    J(:, :, i) = [SoXv*SoXv, 0, 0; 0, SoXteta*SoXteta, 0; 0, 0, SoXh*SoXh];
                    Z_aofL(:, i) = [MoXv, MoXteta, MoXh];
                    Z_fosL(:, i) = Z_aofL(:, i);

                    X(1, i) = normrnd(MoXv, SoXv);
                    X(2, i) = normrnd(MoXteta, SoXteta);
                    X(3, i) = normrnd(MoXh, SoXh);
                 end
            end
            
            if(condition == 2)
                 for i = 1:N;
                    P(:, :, i) = [SoXcond2*SoXcond2, 0; 0, SoYcond2*SoYcond2];
                    J(:, :, i) = [SoXcond2*SoXcond2, 0; 0, SoYcond2*SoYcond2];
                    Z_aofL(:, i) = [MoXcond2, MoYcond2];
                    Z_fosL(:, i) = Z_aofL(:, i);

                    X(1, i) = normrnd(MoXcond2, SoXcond2);
                    X(2, i) = normrnd(MoYcond2, SoYcond2);
                 end
            end
            
            if(condition == 1)
                for i = 1:N
                    X(:, i) = normrnd(MoXcond1, SoXcond1);
                    Z_aofL(:, i) = normrnd(MoXcond1, SoZcond1);
                    Z_fosL(:, i) = Z_aofL(:, i);
                    P(:, :, i) = SoXcond1 * SoXcond1;
                end
            end
            
            % step #0
            Mx = matlabFilters.mathExpectation(X, condition);
            Meps_aofL = matlabFilters.mathExpectation(X - Z_aofL, condition);
            Meps_fosL = matlabFilters.mathExpectation(X - Z_fosL, condition);
            Dx = matlabFilters.dispersion(X, Mx, condition);
            Deps_aofL = matlabFilters.dispersion(X - Z_aofL, Meps_aofL, condition);
            Deps_fosL = matlabFilters.dispersion(X - Z_fosL, Meps_aofL, condition);
            Mz = matlabFilters.mathExpectation(Z_fosL, condition);
            Dz = matlabFilters.dispersion(Z_fosL, Mz, condition);
            

            Sx = ones(NX, 1);
            Seps_aofL = ones(NX, 1);
            Crit_aofL = ones(NX, 1);
            %ICrit_aofL = ones(NX, 1);
            Seps_fosL = ones(NX, 1);
            Crit_fosL = ones(NX, 1);
           % ICrit_fosL = ones(NX, 1);
            for i = 1:NX
                Sx(i) = sqrt(Dx(i, i));
                Seps_aofL(i) = sqrt(Deps_aofL(i, i));
                Crit_aofL(i) = Meps_aofL(i)*Meps_aofL(i) + Deps_aofL(i, i);
                Seps_fosL(i) = sqrt(Deps_fosL(i, i));
                Crit_fosL(i) = Meps_fosL(i)*Meps_fosL(i) + Deps_fosL(i, i);
            end
            ICrit_aofL = Crit_aofL;
            ICrit_fosL = Crit_fosL;

            MxTable(:, :, 1) = Mx;
            Meps_aofLTable(:, :, 1) = Meps_aofL;
            SxTable(:, :, 1) = Sx;
            Seps_aofLTable(:, :, 1) = Seps_aofL;              
            Crit_aofLTable(:, :, 1) = Crit_aofL;
            ICrit_aofLTable(:, :, 1) = ICrit_aofL;   
            Meps_fosLTable(:, :, 1) = Meps_fosL;
            Seps_fosLTable(:, :, 1) = Seps_fosL;              
            Crit_fosLTable(:, :, 1) = Crit_fosL;
            ICrit_fosLTable(:, :, 1) = ICrit_fosL;
            tTable(:, 1) = 0;
                      
            sDeltaTime = sqrt(deltaTime);
             
            
            stepForTrajectory = trajectoryNumber; %from N
            XForTrajectoryPlot(:, 1) = X(:, stepForTrajectory);
            Z_aofLForTrajectoryPlot(:, 1) = Z_aofL(:, stepForTrajectory);
            Z_fosLForTrajectoryPlot(:, 1) = Z_fosL(:, stepForTrajectory);
            
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
   
                    Z_aofL(:, i) = Z_aofL(:, i) + matlabFilters.a(t, Z_aofL(:, i), condition)*deltaTime + ...
                        matlabFilters.K(t, Z_aofL(:, i), P(:, :, i), condition)*(deltaY(:, i) - ...
                        matlabFilters.c(t, Z_aofL(:, i), condition)*deltaTime);
                    P(:, :, i) = P(:, :, i) +  matlabFilters.Ksi(t, Z_aofL(:, i), P(:, :, i), condition)*deltaTime;
                    
                    Z_fosL(:, i) = Z_fosL(:, i) + matlabFilters.a(t, Z_fosL(:, i), condition)*deltaTime + ...
                        matlabFilters.K(t, Z_fosL(:, i), J(:, :, i), condition)*(deltaY(:, i) - ...
                        matlabFilters.c(t, Z_fosL(:, i), condition)*deltaTime);

                     J(:, :, i) = Dx - Dz;
                         
                end
                q(5) = toc;
                
                % end of realization loop

                % calculation of characteristics
                Mx = matlabFilters.mathExpectation(X, condition);
                Meps_aofL = matlabFilters.mathExpectation(X - Z_aofL, condition);
                Meps_fosL = matlabFilters.mathExpectation(X - Z_fosL, condition);
                Dx = matlabFilters.dispersion(X, Mx, condition);
                Deps_aofL = matlabFilters.dispersion(X - Z_aofL, Meps_aofL, condition);
                Deps_fosL = matlabFilters.dispersion(X - Z_fosL, Meps_aofL, condition);
                Mz = matlabFilters.mathExpectation(Z_fosL, condition);
                Dz = matlabFilters.dispersion(Z_fosL, Mz, condition);
                
                Sx = ones(NX, 1);
                Seps_aofL = ones(NX, 1);
                Crit_aofL = ones(NX, 1);
                %ICrit_aofL = ones(NX, 1);
                Seps_fosL = ones(NX, 1);
                Crit_fosL = ones(NX, 1);
                %ICrit_fosL = ones(NX, 1);
                for i = 1:NX
                    Sx(i) = sqrt(Dx(i, i));
                    Seps_aofL(i) = sqrt(Deps_aofL(i, i));
                    Crit_aofL(i) = Meps_aofL(i)*Meps_aofL(i) + Deps_aofL(i,i);
                    Seps_fosL(i) = sqrt(Deps_fosL(i, i));
                    Crit_fosL(i) = Meps_fosL(i)*Meps_fosL(i) + Deps_fosL(i,i);
                end          
                ICrit_aofL = ICrit_aofLTable(:, :, k-1) + (Crit_aofL - ICrit_aofLTable(:, :, k-1)) / k;
                ICrit_fosL = ICrit_fosLTable(:, :, k-1) + (Crit_fosL - ICrit_fosLTable(:, :, k-1)) / k;
                
                MxTable(:, :, k) = Mx;
                Meps_aofLTable(:, :, k) = Meps_aofL;
                SxTable(:, :, k) = Sx;
                Seps_aofLTable(:, :, k) = Seps_aofL;                
                Crit_aofLTable(:, :, k) = Crit_aofL;
                ICrit_aofLTable(:, :, k) = ICrit_aofL;
                Meps_fosLTable(:, :, k) = Meps_fosL;
                Seps_fosLTable(:, :, k) = Seps_fosL;                
                Crit_fosLTable(:, :, k) = Crit_fosL;
                ICrit_fosLTable(:, :, k) = ICrit_fosL;   
                tTable(:, k) = t;
                
                t = (k) * deltaTime;
                
                XForTrajectoryPlot(:, k) = X(:, stepForTrajectory); 
                Z_aofLForTrajectoryPlot(:, k) = Z_aofL(:, stepForTrajectory); 
                Z_fosLForTrajectoryPlot(:, k) = Z_fosL(:, stepForTrajectory);                 
            end
            % end of time loop
            
            delete(progressBarFunction);

            % displaying the table
            if(condition == 1)
                disp('     time       Mx     Meps_aofL    Sx    Seps_aofL   Crit_aofL  ICrit_aofL');
                answMatrix = [tTable, MxTable, Meps_aofLTable, SxTable, Seps_aofLTable, Crit_aofLTable, ICrit_aofLTable];
                for i = 1:K
                    disp(answMatrix(:,:,i));
                end
                disp('     time       Mx     Meps_fosL    Sx    Seps_fosL   Crit_fosL  ICrit_fosL');
                answMatrix = [tTable, MxTable, Meps_fosLTable, SxTable, Seps_fosLTable, Crit_fosLTable, ICrit_fosLTable];
                for i = 1:K
                    disp(answMatrix(:,:,i));
                end
            end
            if(condition == 2)
                disp('     time       Mx1       Mx2   MepsAofL1  MepsAofL2    Sx1      Sx2    SepsAofL1 SepsAofL2  CritAofL1 CritAofL2 ICritAofL1 ICritAofL2');
                answMatrix = [tTable, MxTable(1, :, :), MxTable(2, :, :), Meps_aofLTable(1, :, :), Meps_aofLTable(2, :, :), ...
                    SxTable(1, :, :), SxTable(2, :, :), Seps_aofLTable(1, :, :), Seps_aofLTable(2, :, :), Crit_aofLTable(1, :, :), Crit_aofLTable(2, :, :), ...
                    ICrit_aofLTable(1, :, :), ICrit_aofLTable(2, :, :)];
                for i = 1:K
                    disp(answMatrix(:,:,i));
                end
                
                disp('     time       Mx1       Mx2   MepsFosL1  MepsFosL2    Sx1      Sx2    SepsFosL1 SepsFosL2  CritFosL1 CritFosL2 ICritFosL1 ICritFosL2');
                answMatrix = [tTable, MxTable(1, :, :), MxTable(2, :, :), Meps_fosLTable(1, :, :), Meps_fosLTable(2, :, :), ...
                    SxTable(1, :, :), SxTable(2, :, :), Seps_fosLTable(1, :, :), Seps_fosLTable(2, :, :), Crit_fosLTable(1, :, :), Crit_fosLTable(2, :, :), ...
                    ICrit_fosLTable(1, :, :), ICrit_fosLTable(2, :, :)];
                for i = 1:K
                    disp(answMatrix(:,:,i));
                end
            end
            
%             plot Sx, Seps, Crit
            SxToPlot = zeros(1, NX);
            SeAofLToPlot = zeros(1, NX);
            SeFosLToPlot = zeros(1, NX);
%             sqrtCrit_aofLToPlot = zeros(1, NX);
%             sqrtCrit_fosLToPlot = zeros(1, NX);
%             ICrit_aofLToPlot = zeros(1, NX);
%             ICrit_fosLToPlot = zeros(1, NX);
            MxToPlot = zeros(1, NX);
            MeAofLToPlot = zeros(1, NX);
            MeFosLToPlot = zeros(1, NX);
            tToPlot = zeros(1, NX);
            
            figure('name', 'Sx, Se');
            for j = 1:NX
                for i = 1:K
                    SxToPlot(i) = SxTable(j, 1, i);
                    SeAofLToPlot(i) = Seps_aofLTable(j, 1, i);
                    SeFosLToPlot(i) = Seps_fosLTable(j, 1, i);
%                     sqrtCrit_aofLToPlot(i) = sqrt(Crit_aofLTable(j, 1, i));
%                     sqrtCrit_fosLToPlot(i) = sqrt(Crit_fosLTable(j, 1, i));
%                     ICrit_aofLToPlot(i) = ICrit_aofLTable(j, 1, i);
%                     ICrit_fosLToPlot(i) = ICrit_fosLTable(j, 1, i);

                    tToPlot(i) = tTable(1, 1, i);
                end               
                if(condition == 3) 
                    if(j == 2)
                        SxToPlot = SxToPlot*180/pi;
                        SeAofLToPlot = SeAofLToPlot*180/pi;
                        SeFosLToPlot = SeFosLToPlot*180/pi;
                    else
                        SxToPlot = SxToPlot;
                        SeAofLToPlot = SeAofLToPlot;
                        SeFosLToPlot = SeFosLToPlot;
                    end
                end

                subplot(condition,1,j);
                plot(tToPlot.', SxToPlot.');
                hold on;
                if(buildLAOF)
                    plot(tToPlot.', SeAofLToPlot.', 'o-', 'color', [240, 110, 50]/255);
                    hold on;
                end

                if(buildLFOS)
                    plot(tToPlot.', SeFosLToPlot.', 'kx--');
                end
                
                % legend format
                if(buildLAOF && buildLFOS)
                    legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'), strcat('Se', int2str(j), ' FOS-L'));
                end
                if(buildLAOF && ~buildLFOS)
                    legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' AOF-L'));
                end
                if(~buildLAOF && buildLFOS)
                    legend(strcat('Sx', int2str(j)), strcat('Se', int2str(j), ' FOS-L'));
                end
                if(~buildLAOF && ~buildLFOS)
                    legend(strcat('Sx', int2str(j)));
                end

                if(condition == 1)
                    str = 'Скалярный пример: ';
                end
                if(condition == 2)
                    str = 'Осциллятор Ван-дер-Поля: ';
                end
                if(condition == 3)
                    str = 'Спуск ЛА на планету: ';
                end
                title(strcat(str, ' N=', int2str(N), ', deltaTime=', num2str(deltaTime), '. '));
                xlabel('t'); 
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
                        title(strcat('Траектория j= ', int2str(stepForTrajectory), ' из ', int2str(N)));
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
                        title(strcat('Траектория j= ', int2str(stepForTrajectory), ' из ', int2str(N)));
                    end
                end
            end
            
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
                    
                    title(strcat('max|Mx|=', num2str(maxMx(1, j), '%10.2e'), '; max|Me AOF-L|=', num2str(maxMeAOFL(1, j), '%10.2e'), '; max|Me FOS-L|=', num2str(maxMeFOSL(1, j), '%10.2e')));
                end
            end
            

     
        end        
        % end of Kalman filter function
        
        
       
       
       
       
    end
    
end

