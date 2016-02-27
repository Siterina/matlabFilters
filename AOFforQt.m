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
        
        function [ A ] = A(t, x, condition)                               % ������� ����� ������ ����������� a(t, x) �� x
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
        
        function [ G ] = G(t, x, condition)                               % ������� ����� ������ ����������� �(t, x) �� x
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

            
            a0 = 0;                                                        % a, g - ��������� �������
            a1 = -0.5;
            a3 = -1;
            alpha1 = 1;                                                    % alpha, beta - ��������� ����������
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
            
            MoX = 0.1;                                                     % �������������� ���. ���. �������, ���������� � �������
            MoY = 0.0;
            SoX = 0.6;
            SoY = 0.0;
            SoZ = 0.0; 
        end
        
        
        function[ ] = KalmanFilter(condition)
            
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
            MepsLTable = ones(NX, 1, K);
            SepsLTable = ones(NX, 1, K);
            CritLTable = ones(NX, 1, K);
            ICritLTable = ones(NX, 1, K);
            
            
            t = deltaTime;
            X = ones(NX, N);
            deltaY = ones(NX, N);
            Z = ones(NX, N);
            P = ones(NX, NY, N);
            
            
            rng(10);
            if(condition == 2)
                for i = 1:N
                    P(:, :, i) = [5, 0; 0, 5];
                    Z(:, i) = [10, -3];
                    X(1, i) = normrnd(10, sqrt(5));
                    X(2, i) = normrnd(-3, sqrt(5));

                end
            end
            
            if(condition == 1)
                for i = 1:N
                    X(:, i) = normrnd(MoX, SoX);
                    Z(:, i) = normrnd(MoX, SoZ);
                    P(:, :, i) = SoX * SoX;
                end
            end
            
            % step #0
                Mx = AOFforQt.mathExpectation(X, condition);
                MepsL = AOFforQt.mathExpectation(X - Z, condition);
                Dx = AOFforQt.dispersion(X, Mx, condition);
                DepsL = AOFforQt.dispersion(X - Z, MepsL, condition);
                Sx = ones(NX, 1);
                Sx(1) = sqrt(Dx(1, 1));
                Sx(NX) = sqrt(Dx(NX, NX));
                
                SepsL = ones(NX, 1);
                SepsL(1) = sqrt(DepsL(1, 1));
                SepsL(NX) = sqrt(DepsL(NX, NX));
                
                CritL = ones(NX, 1);
                CritL(1) = MepsL(1)*MepsL(1) + DepsL(1, 1);
                CritL(NX) = MepsL(NX)*MepsL(NX) + DepsL(NX, NX);

                ICritL = ones(NX, 1);
                ICritL = CritL;

                MxTable(:, :, 1) = Mx;
                MepsLTable(:, :, 1) = MepsL;
                SxTable(:, :, 1) = Sx;
                SepsLTable(:, :, 1) = SepsL;
                
                CritLTable(:, :, 1) = CritL;
                ICritLTable(:, :, 1) = ICritL;
                
                tTable(:, 1) = 0;
            
            
            sDeltaTime = sqrt(deltaTime);

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
   
                    Z(:, i) = Z(:, i) + AOFforQt.a(t, Z(:, i), condition)*deltaTime + ...
                        AOFforQt.K(t, Z(:, i), P(:, :, i), condition)*(deltaY(:, i) - ...
                        AOFforQt.c(t, Z(:, i), condition)*deltaTime);
                    P(:, :, i) = P(:, :, i) +  AOFforQt.Ksi(t, Z(:, i), P(:, :, i), condition)*deltaTime;
                         
                end
                % end of realization loop

                % calculation of characteristics
                Mx = AOFforQt.mathExpectation(X, condition);
                MepsL = AOFforQt.mathExpectation(X - Z, condition);
                Dx = AOFforQt.dispersion(X, Mx, condition);
                DepsL = AOFforQt.dispersion(X - Z, MepsL, condition);
                Sx = ones(NX, 1);
                Sx(1) = sqrt(Dx(1, 1));
                Sx(NX) = sqrt(Dx(NX, NX));
                
                SepsL = ones(NX, 1);
                SepsL(1) = sqrt(DepsL(1, 1));
                SepsL(NX) = sqrt(DepsL(NX, NX));
                
                CritL = ones(NX, 1);
                CritL(1) = MepsL(1)*MepsL(1) + DepsL(1,1);
                CritL(NX) = MepsL(NX)*MepsL(NX) + DepsL(NX,NX);

                ICritL = ones(NX, 1);
                ICritL = ICritLTable(:, :, k-1) + (CritL - ICritLTable(:, :, k-1)) / k;
                

                MxTable(:, :, k) = Mx;
                MepsLTable(:, :, k) = MepsL;
                SxTable(:, :, k) = Sx;
                SepsLTable(:, :, k) = SepsL;
                
                CritLTable(:, :, k) = CritL;
                ICritLTable(:, :, k) = ICritL;
                
                tTable(:, k) = t;
                
                t = (k) * deltaTime;
               
            end
            % end of time loop

            % displaying the table
            if(condition == 1)
                disp('     time       Mx       MepsL      Sx      SepsL      CritL     ICritL');
                answMatrix = [tTable, MxTable, MepsLTable, SxTable, SepsLTable, CritLTable, ICritLTable];
                for i = 1:K
                    disp(answMatrix(:,:,i));
                end
                for i = 1:K               
                    SxToPlot(i) = SxTable(1, 1, i);
                    SeToPlot(i) = SepsLTable(1, 1, i);
                    sqrtCritToPlot(i) = sqrt(CritLTable(1, 1, i));
                    tToPlot(i) = tTable(1, 1, i);
                end
            
                figure('name', 'result');
               
                plot(tToPlot.', SxToPlot.');
                hold on;
                plot(tToPlot.', SeToPlot.');
                hold on;
                plot(tToPlot.', sqrtCritToPlot.');
                legend('SxL', 'SeL', 'sqrtCritL');
            end
            if(condition == 2)
                disp('     time       Mx1       Mx2     MepsL1     MepsL2     Sx1       Sx2     SepsL1    SepsL2    CritL1    CritL2    ICritL1   ICritL2');
                answMatrix = [tTable, MxTable(1, :, :), MxTable(2, :, :), MepsLTable(1, :, :), MepsLTable(2, :, :), ...
                    SxTable(1, :, :), SxTable(2, :, :), SepsLTable(1, :, :), SepsLTable(2, :, :), CritLTable(1, :, :), CritLTable(2, :, :), ...
                    ICritLTable(1, :, :), ICritLTable(2, :, :)];
                for i = 1:K
                    disp(answMatrix(:,:,i));
                end
                for i = 1:K
                Sx1ToPlot(i) = SxTable(1, 1, i);
                Se1ToPlot(i) = SepsLTable(1, 1, i);
                sqrtCrit1ToPlot(i) = sqrt(CritLTable(1, 1, i));
                
                Sx2ToPlot(i) = SxTable(2, 1, i);
                Se2ToPlot(i) = SepsLTable(2, 1, i);
                sqrtCrit2ToPlot(i) = sqrt(CritLTable(2, 1, i));
                tToPlot(i) = tTable(1, 1, i);
            end
            
            figure('name', 'result');
            subplot(1,2,1);
            plot(tToPlot.', Sx1ToPlot.');
            hold on;
            plot(tToPlot.', Se1ToPlot.', 'color', [240, 110, 50]/255);
            hold on;
            plot(tToPlot.', sqrtCrit1ToPlot.', 'color', [0, 110, 50]/255);
            legend('SxL1', 'SeL1', 'sqrtCritL1');
            
            subplot(1,2,2);
            plot(tToPlot.', Sx2ToPlot.');
            hold on;
            plot(tToPlot.', Se2ToPlot.', 'color', [240, 110, 50]/255);
            hold on;
            plot(tToPlot.', sqrtCrit2ToPlot.', 'color', [0, 110, 50]/255);
            legend('SxL2', 'SeL2', 'sqrtCritL2');
                
            end

            
           
            
            
        end        
        % end of Kalman filter function
        
        
       
       
       
       
    end
    
end
