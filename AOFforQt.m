classdef AOFforQt
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
    end
    
    methods(Static)
        %% Kalman Filter
        
        % Functions-coefficients of the system
        
        function [ a ] = a(t, x, condition)
            if(condition == 2)
                global omega;
                global alpha;
                global beta;
               
                a(1, 1) = x(2);
                a(2, 1) = -omega*omega*x(1) + alpha*x(2) - alpha*beta*x(1)*x(1)*x(2);
            end           
            if(condition == 1)
                global a0;
                global a1;
                global a3;
                a = a0 + a1 * x + a3 * x * x * x;
            end
        end
        
        function [ c ] = c(t, x, condition)
            if(condition == 2)
                c(1, 1) = x(1);
                c(2, 1) = x(2);
            end            
            if(condition == 1)
                global alpha1;
                global alpha2;
                c = alpha1 * x + alpha2 * x * x;
            end
        end
        
        function [ B ] = B(t, x, condition)
            if(condition == 2)
                B(1, 1) = 0;
                B(1, 2) = 0;
                B(2, 1) = 0;
                B(2, 2) = x(1);
            end            
            if(condition == 1)
                global g0;
                global g1;
                B = g0 + g1 * x;
            end
        end
        
        function [ D ] = D(t, x, condition)
            if(condition == 2)
                global d11;
                global d22;
                D(1, 1) = d11;
                D(1, 2) = 0;
                D(2, 1) = 0;
                D(2, 2) = d22;   
            end           
            if(condition == 1)
                global beta1;
                global beta2;
                D = beta1 + beta2 * x;
            end
        end
        
        function [ R ] = R(t, x, condition)
            D = AOFforQt.D(t, x, condition);
            R = D * D.';
        end
        
        function [ Q ] = Q(t, x, condition)
            B = AOFforQt.B(t, x, condition);
            Q = B * B.';
        end
        
        
        % Support functions for Kalman filter
        
        function [ A ] = A(t, x, condition)                               % матрица якоби первых производных a(t, x) по x
            if(condition == 2)
                global omega;
                global alpha;
                global beta;
                A(1, 1) = 0;
                A(1, 2) = 1;
                A(2, 1) = -omega*omega - 2*alpha*beta*x(1)*x(2);
                A(2, 2) = alpha*(1 - beta*x(1)*x(1));
            end           
            if(condition == 1)
                global a0;
                global a1;
                global a3;
                A = a1 + 3 * a3 * x * x;
            end
        end
        
        function [ G ] = G(t, x, condition)                               % матрица якоби первых производных с(t, x) по x
           if(condition == 2)
                G(1, 1) = 1;
                G(1, 2) = 0;
                G(2, 1) = 0;
                G(2, 2) = 1;
           end            
           if(condition == 1)
                global alpha1;
                global alpha2;
                G = alpha1 + 2 * alpha2 * x;
           end
        end
        
        
        % structure functions
        
        function [ K ] = K(t, x, p, condition) 
            K = p * AOFforQt.G(t, x, condition).' * inv(AOFforQt.R(t, x, condition));
        end
        
        function [ Ksi ] = Ksi(t, x, p, condition)                        
            Ksi = AOFforQt.A(t, x, condition) * p + p * AOFforQt.A(t, x, condition).' + ...
                AOFforQt.Q(t, x, condition) - AOFforQt.K(t, x, p, condition) * ... 
                AOFforQt.R(t, x, condition) * AOFforQt.K(t, x, p, condition).';
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
            if(condition == 2)
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
        
        
        
        
        function [ ] = initialConditions()
            global a0;
            global a1;
            global a3;
            global alpha1;
            global alpha2;
            global g0;
            global g1;
            global beta1;
            global beta2;
            
            global omega;
            global alpha;
            global beta;
            global d11;
            global d22;
            
            global deltaTime;
            global N;
            global K;
            
            global MoX;
            global MoY;
            global SoX;
            global SoY;
            global SoZ;

            
            a0 = 0;                                                        % a, g - параметры объекта
            a1 = -0.5;
            a3 = -1;
            alpha1 = 1;                                                    % alpha, beta - параметры измерител€
            alpha2 = 0.5;
            g0 = 1.5;
            g1 = 0;
            beta1 = 0.5;
            beta2 = 0.1;
            
            omega = 0.1*pi;
            alpha = 2;
            beta = 1;
            d11 = 0.1;
            d22 = 0.1;
            
            N = 1500;
            K = 7;
            deltaTime = 0.0005;
            
            MoX = 0.1;                                                     % ’арактеристики нач. усл. объекта, измерител€ и фильтра
            MoY = 0.0;
            SoX = 0.6;
            SoY = 0.0;
            SoZ = 0.0; 
        end
        
        
        function[ ] = main(condition)
            
            % Initial conditions
            global deltaTime;
            global N;
            global K;
            global MoX;
            global MoY;
            global SoX;
            global SoY;
            global SoZ;

            AOFforQt.initialConditions();
            if(condition == 2)
                NX = 2;
                NY = 2;
            end
           if(condition ==1)
                NX = 1;
                NY = 1;
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
            deltaY = ones(NX, N);
            Z_aofL = ones(NX, N);
            Z_fosL = ones(NX, N);
            P = ones(NX, NY, N);
            J = ones(NX, NY, N);
                      
            rng(10);
            if(condition == 2)
                for i = 1:N
                    P(:, :, i) = [5, 0; 0, 5];
                    J(:, :, i) = [5, 0; 0, 5];
                    Z_aofL(:, i) = [10, -3];
                    Z_fosL(:, i) = Z_aofL(:, i);
                    X(1, i) = normrnd(10, sqrt(5));
                    X(2, i) = normrnd(-3, sqrt(5));

                end
            end
            
            if(condition == 1)
                for i = 1:N
                    X(:, i) = normrnd(MoX, SoX);
                    Z_aofL(:, i) = normrnd(MoX, SoZ);
                    Z_fosL(:, i) = Z_aofL(:, i);
                    P(:, :, i) = SoX * SoX;
                end
            end
            
            % step #0
            Mx = AOFforQt.mathExpectation(X, condition);
            Meps_aofL = AOFforQt.mathExpectation(X - Z_aofL, condition);
            Meps_fosL = AOFforQt.mathExpectation(X - Z_fosL, condition);
            Dx = AOFforQt.dispersion(X, Mx, condition);
            Deps_aofL = AOFforQt.dispersion(X - Z_aofL, Meps_aofL, condition);
            Deps_fosL = AOFforQt.dispersion(X - Z_fosL, Meps_aofL, condition);

            Sx = ones(NX, 1);
            Seps_aofL = ones(NX, 1);
            Crit_aofL = ones(NX, 1);
            ICrit_aofL = ones(NX, 1);
            Seps_fosL = ones(NX, 1);
            Crit_fosL = ones(NX, 1);
            ICrit_fosL = ones(NX, 1);
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
             
            
            stepForTrajectory = 777; %from 1500
            XForTrajectoryPlot(:, 1) = X(:, stepForTrajectory);
            
            % time loop
            for k = 2:K                
                % realization loop
                for i = 1:N
                    
                    deltaV1 = normrnd(0, sDeltaTime, NX, 1);
                    deltaV2 = normrnd(0, sDeltaTime, NX, 1); 
                   
                    X(:, i) = X(:, i) + AOFforQt.a(t, X(:, i), condition)*deltaTime + ...
                        AOFforQt.B(t, X(:, i), condition)*deltaV1;

                    deltaY(:, i) = AOFforQt.c(t, X(:, i), condition)*deltaTime + ...
                        AOFforQt.D(t, X(:, i), condition)*deltaV2;                             
   
                    Z_aofL(:, i) = Z_aofL(:, i) + AOFforQt.a(t, Z_aofL(:, i), condition)*deltaTime + ...
                        AOFforQt.K(t, Z_aofL(:, i), P(:, :, i), condition)*(deltaY(:, i) - ...
                        AOFforQt.c(t, Z_aofL(:, i), condition)*deltaTime);
                    P(:, :, i) = P(:, :, i) +  AOFforQt.Ksi(t, Z_aofL(:, i), P(:, :, i), condition)*deltaTime;
                    
                    Z_fosL(:, i) = Z_fosL(:, i) + AOFforQt.a(t, Z_fosL(:, i), condition)*deltaTime + ...
                        AOFforQt.K(t, Z_fosL(:, i), J(:, :, i), condition)*(deltaY(:, i) - ...
                        AOFforQt.c(t, Z_fosL(:, i), condition)*deltaTime);
                    Mx = AOFforQt.mathExpectation(X, condition);
                    Mz = AOFforQt.mathExpectation(Z_fosL, condition);
                    Dx = AOFforQt.dispersion(X, Mx, condition);
                    Dz = AOFforQt.dispersion(Z_fosL, Mz, condition);
                    J(:, :, i) = Dx - Dz;
                         
                end
                % end of realization loop

                % calculation of characteristics
                Mx = AOFforQt.mathExpectation(X, condition);
                Meps_aofL = AOFforQt.mathExpectation(X - Z_aofL, condition);
                Meps_fosL = AOFforQt.mathExpectation(X - Z_fosL, condition);
                Dx = AOFforQt.dispersion(X, Mx, condition);
                Deps_aofL = AOFforQt.dispersion(X - Z_aofL, Meps_aofL, condition);
                Deps_fosL = AOFforQt.dispersion(X - Z_fosL, Meps_aofL, condition);
                
                Sx = ones(NX, 1);
                Seps_aofL = ones(NX, 1);
                Crit_aofL = ones(NX, 1);
                ICrit_aofL = ones(NX, 1);
                Seps_fosL = ones(NX, 1);
                Crit_fosL = ones(NX, 1);
                ICrit_fosL = ones(NX, 1);
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
                
            end
            % end of time loop

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
            
            % plot Sx, Seps, Crit
            figure('name', 'Sx; Seps; Crit');
            for j = 1:NX
                for i = 1:K
                    SxToPlot(i) = SxTable(j, 1, i);
                    SeAofLToPlot(i) = Seps_aofLTable(j, 1, i);
                    SeFosLToPlot(i) = Seps_fosLTable(j, 1, i);
                    sqrtCrit_aofLToPlot(i) = sqrt(Crit_aofLTable(j, 1, i));
                    sqrtCrit_fosLToPlot(i) = sqrt(Crit_fosLTable(j, 1, i));
                    tToPlot(i) = tTable(1, 1, i);
                end                

                subplot(1,2,j);
                plot(tToPlot.', SxToPlot.');
                hold on;
                plot(tToPlot.', SeAofLToPlot.', 'color', [240, 110, 50]/255);
                hold on;
                plot(tToPlot.', sqrtCrit_aofLToPlot.', 'color', [0, 110, 50]/255);
                hold on;
                plot(tToPlot.', SeFosLToPlot.', 'k');
                hold on;
                plot(tToPlot.', sqrtCrit_fosLToPlot.', 'color', [200, 180, 10]/255);
                legend(strcat('Sx', int2str(j)), strcat('SeAOF-L', int2str(j)), strcat('sqrtCrit-aofL', int2str(j)), strcat('SeFOS-L', int2str(j)), strcat('sqrtCrit-fosL', int2str(j)));
            end     
           
            if(condition == 2)
               % plot trajectory 
               figure('name', 'trajectory; X(t)');
               subplot(2,2,1);
               plot(XForTrajectoryPlot(1, :).', XForTrajectoryPlot(2, :).');
               xlabel('X1'); 
               ylabel('X2');
               subplot(2,2,3);
               plot(tToPlot.'.', XForTrajectoryPlot(1, :).');
               xlabel('t'); 
               ylabel('X1');
               subplot(2,2,4);
               plot(tToPlot.'.', XForTrajectoryPlot(2, :).');
               xlabel('t'); 
               ylabel('X2');
            end
               
        end        
        % end of Kalman filter function
        
        
       
       
       
       
    end
    
end

