
%submitted by Pilli SAI SNEHITH (234103433)
clc
clear all;
format long
% Define the load matrix N
N = [100; 0; 0; 0; 0; 0];

% Material properties
E1 = 38.6;
E2 = 8.27;
nu12 = 0.28;
nu21 = (nu12*E2)/E1;
G12 = 4.14;
uts1=1062;
ucs1=610;
uts2=31;
ucs2=118;
uss=72;
alpha_1=8.6*1.0e-6;
alpha_2=22.1*1.0e-6;
alpha_12=0;
Tdiff=50;
C=0;
% Laminate properties
theta = [0,45,-45,90,90,-45,45,0]; % Fiber orientation angles in degrees for each layer
r=length(unique(theta));
thickness =0.125; % Thickness of each layer
num_layers = length(theta);T=thickness*0.5*num_layers;
thetaf=[];
% Z coordinates for each layer
Z = linspace(-T,T,num_layers+1);
%disp(Z)
% Calculate Q_bar matrix for the current layer
q11 = E1 / (1 - nu12 * nu21);
q22 = E2 / (1 - nu12 * nu21);
q12 = nu12 * E2 / (1 - nu12 * nu21);
q66 = G12;
q16=0;
q26=0;
q=[q11, q12, q16;q12, q22, q26;q16 q26, q66];

% Initialize zeros arrays to store Q_bar matrices, A matrices, B matrices, and D matrices for each layer
Q_bars = zeros(3, 3);
alpha1 = [alpha_1;alpha_2;alpha_12];

hx=0;hy=0;hxy=0;
hmx=0;hmy=0;hmxy=0;
NH = [hx;hy;hxy]; % Initialize thermal forces for each layer
MH = [hmx;hmy;hmxy]; % Initialize thermal moments for each layerr

for j=1:r
    fprintf("case%d\n",j);
    if isempty(setdiff(theta,thetaf))
        break
    else
        A = zeros(3, 3);
        B = zeros(3, 3);
        D = zeros(3, 3);
        nx=0;ny=0;nxy=0;
        mx=0;my=0;mxy=0;
        NT = [nx;ny;nxy]; % Initialize thermal forces for each layer
        MT = [mx;my;mxy]; % Initialize thermal moments for each layerr

        for i = 1:num_layers
            if any(thetaf==theta(i))

                continue
            else
                c = cosd(theta(i));
                s = sind(theta(i));

                T=[c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2 - s^2];
                alphaT= inv(T) * alpha1;
                x=alphaT(1);
                y=alphaT(2);
                z=2*alphaT(3);
                alpha=[x;y;z];
                Q11 = q11 * c ^ 4 + 2 * (q12 + 2 * q66) * s ^ 2 * c ^ 2 + q22 * s ^ 4;
                Q12 = (q11 + q22 - 4 * q66) * c ^ 2 * s ^ 2 + q12 * (s ^ 4 + c ^ 4);
                Q22 = q11 * s ^ 4 + 2 * (q12 + 2 * q66) * c ^ 2 * s ^ 2 + q22 * c ^ 4;
                Q16 = (q11 - q12 - 2 * q66) * s * c ^ 3 - (q22 - q12 - 2 * q66) * s ^ 3 * c;
                Q26 = (q11 - q12 - 2 * q66) * s ^ 3 * c - (q22 - q12 - 2 * q66) * s * c ^ 3;
                Q66 = (q11 + q22 - 2 * q12 - 2 * q66) * c ^ 2 * s ^ 2 + q66 * (s ^ 4 + c ^ 4);


                Q_bar = [Q11, Q12, Q16; Q12, Q22, Q26; Q16, Q26, Q66];
                % Store Q_bar matrix for the current layer
                % Q_bars{i} = Q_bar;
                Zk_minus_1 = Z(i); % Distance of upper surface of the current layer from mid-surface
                Zk = Z(i+1); % Distance of lower surface of the current layer from mid-surface

                % Calculate thermal forces for the current layer
                NT =NT+ Tdiff * Q_bar* alpha * (Zk - Zk_minus_1) ;
                % Calculate thermal moments for the current layer
                MT =MT+ 0.5 * Tdiff * Q_bar* alpha * ((Zk)^2 - (Zk_minus_1)^2) ;
                %claculate hygroscpic stress
                NH =NH+ C * Q_bar* alpha * (Zk - Zk_minus_1) ;
                % Calculate hygroscpic moments for the current layer
                MH =MH+ 0.5 *C * Q_bar* alpha * ((Zk)^2 - (Zk_minus_1)^2) ;
                % Calculate A matrix for the current layer
                Zk_minus_1 = Z(i); % Distance of upper surface of the current layer from mid-surface
                Zk = Z(i+1); % Distance of lower surface of the current layer from mid-surface
                A = A+Q_bar * (Zk - Zk_minus_1);

                % Store A matrix for the current layer
                %A_matrices{i} = A;

                % Calculate B matrix for the current layer
                B = B+0.5 * Q_bar * ((Zk)^2 - (Zk_minus_1)^2);

                % Store B matrix for the current layer
                %B_matrices{i} = B;

                % Calculate D matrix for the current layer
                D = D+(1/3) * Q_bar * ((Zk)^3 - (Zk_minus_1)^3);

                % Store D matrix for the current layer
                %D_matrices{i} = D;
            end
        end
        % Combine thermal forces and moments into a single vector
        NFT = [NT*1.e3; MT*1.e3];


        % Display thermal forces
        % disp('Thermal forces and moments:');
        % disp(NFT);

        % Form the matrix [A B; B D]
        matrix_ABBD = [A, B; B, D];

        % Display matrix ABBD
        % disp('Matrix [A B; B D]:');
        %disp(matrix_ABBD);


        % Calculate the strain matrix e
        e = inv(matrix_ABBD) * N;
        er= inv(matrix_ABBD) * NFT*1.0e-3;
        eh=inv(matrix_ABBD) * NFT*1.0e-3;

        % Display the resulting strain matrix e
        % disp('Resulting strain matrix [ex; ey; exy; Kx; Ky; Kxy]:');
        % disp(e);
        % fprintf('Resulting thermal strain matrix ');
        % disp(er);


        % Initialize zeros arrays to store strain matrices above and below each layer
        e_above = zeros(num_layers, 1);
        e_below = zeros(num_layers, 1);
        E_residual = zeros(num_layers, 1);
        et_above = zeros(num_layers, 1);
        et_below = zeros(num_layers, 1);

        % Extract E and K from strain matrix e
        E = e(1:3); % Reference strain
        K = e(4:6); % Curvature
        % Extract E and K from strain matrix e
        Er = er(1:3); % Reference strain
        Kr = er(4:6); % Curvature
        a=0;
        for i = 1:num_layers

            % Calculate strain above the current layer
            if Z(i)<0
                e_above_i = E + Z(i) * K;
                %e_above{i} = e_above_i;
            else
                % Calculate strain below the current layer
                e_below_i = E + Z(i+1) * K;
                %e_below{i} = e_below_i;
                %disp(e_below{i})
            end
        end


        for i = 1:num_layers
            % Extract the elements from e_below_i
            exb = e_below_i(1);
            exa= e_above_i(1);
            eyb = e_below_i(2);
            eya= e_above_i(2);
            exy_b =0.5* e_below_i(3); % Modify the third element
            exy_a =0.5* e_above_i(3);
            % Store the modified elements in E_below_i matrix
            E_below_i = [exb; eyb; exy_b];
            E_above_i = [exa; eya; exy_a];

            % Store E_below_i matrix for the current layer
            % E_below{i} = E_below_i;
            % E_above{i} = E_above_i;


        end
        d=0;
        n=0;
        S=zeros(3,num_layers);
        for i = 1:num_layers
            c = cosd(theta(i));
            s = sind(theta(i));
            T=[c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2 - s^2];
            %fprintf('principal strains in layer%d\n',i);
            if (Z(i)<0)
                e_i=T*e_above_i;
            else
                e_i=T*e_below_i;

            end
            for k=1:3
                e1_i(1)=e_i(1);
                e1_i(2)=e_i(2);
                e1_i(3)=2*e_i(3);
            end
            if any(thetaf==theta(i))
                S(:,i)=zeros(3,1);
                continue
            else

                S(:,i)=S(:,i)+q*transpose(e1_i);

                % fprintf('Layer %d:\n', i);
                % disp(e1_i)


                if S(1,i)>0
                    sr1 = S(1,i) / uts1;
                else
                    sr1 = abs(S(1,i)) / ucs1;
                end

                if S(2,i)>0
                    sr2 = S(2,i) / uts2;
                else
                    sr2 =abs( S(2,i)) / ucs2;
                end

                sr3 = abs(S(3,i)) / uss;

               % fprintf('Stress ratios: sr1=%.4f, sr2=%.4f, sr3=%.4f\n', sr1, sr2, sr3);

                % Check for failure
                max_sr = max([sr1, sr2, sr3]);


                l=max_sr;
                if (l>=d)
                    d=l;
                    a=i;
                    %disp(e1_i)

                    if sr1 == max_sr
                        %fprintf("failure due to longitudinal tension\n");
                        n=1;
                        %disp(n)
                    elseif sr2 == max_sr
                        % fprintf("failure due to transverse tension\n");
                        n=2;
                        % disp(n)
                    else
                        n=3;
                        %fprintf("failure due to shear\n");
                        %disp(n)
                    end
                    s1=S(n,i);

                end
            end
        end
        fprintf('principal stress in layers');
        disp(S)
         
        % fprintf('%d  lamina fails \n', a);
        Nx1=N/d;
        if j==1
            LPF=Nx1;
        end
        % disp(d)
        % disp(Nx1)
        % disp(s1)
        % disp(n)


        %%thermal
        b=0;

        Rs=zeros(3,num_layers);
        for i = 1:num_layers
            c = cosd(theta(i));
            s = sind(theta(i));
            R=[1,0,0;0,1,0;0,0,2];
            T=[c^2, s^2, 2*s*c; s^2, c^2, -2*s*c; -s*c, s*c, c^2 - s^2];
            alphaT= inv(T) * alpha1;
            x=alphaT(1);
            y=alphaT(2);
            z=2*alphaT(3);
            alpha=[x;y;z];
            ET_i=Tdiff*alpha;
            %fprintf("thermal strain in %dth lamina\n  ",i);
            %disp(ET_{i})
            % Calculate strain above the current layer

            if Z(i)<0
                et_above_i = Er + Z(i)*Kr-ET_i;
                %et_above{i} = et_above_i-ET_{i};
                %fprintf("residual strain in %dth are\n",i)
                %disp(et_above{i})
                Rs(:,i)=Rs(:,i)+q*R*T*inv(R)*(et_above_i).*1.0e3;
                if(i==a)
                    b=Rs(n,i);
                    % fprintf("residual stress in layer %d are\n",i)
                    % disp(b)
                end
            else
                % Calculate strain below the current layer
                et_below_i = Er + Z(i+1)*Kr -ET_i;
                %et_below{i} = et_below_i-ET_{i};
                %disp(et_below{i})
                Rs(:,i)=Rs(:,i)+q*R*T*inv(R)*(et_below_i).*1.0e3;
                if(i==a)
                    b=Rs(n,i);
                    %fprintf("residual stress in layer %d are\n",i)
                    %disp(b)
                end
            end
        end

        %disp(Rs)

        st=s1+b;
        % disp(st);
        if n==1
            if s1>0
               % fprintf("strength is:",uts1)
                x1=uts1-b;
                d1=s1/x1;
            else
                %fprintf("strength is:",ucs1)
                x1=ucs1-b;
                d1=s1/x1;
            end
        end
        if n==2
            if s1>0
               % fprintf("strength is:",uts2)
                x1=uts2-b;
                d1=s1/x1;
            else
                %fprintf("strength is:%f",ucs2)
                x1=ucs2-b;
                d1=s1/x1;
            end
        end
        if n==3
           %fprintf("strength is:%f",uss)
            x1=uss-s1;
            d1=s1/x1;
        end
        %fprintf('new load in failure lamina after considering thermal stress')
        %disp(x1)
        % fprintf("stress ratio thermal")
        % disp(d1)
        Nx1=N/d1;
        %fprintf('Failure  load considering thermal stress')
        %disp(Nx1)
        thetaf=[thetaf,theta(a)];
       % fprintf('%d angle lamina fails\n',thetaf(j))
        N=Nx1;

    end

end
fprintf("order of failure of theta")
disp(thetaf)
fprintf("FPF is");
disp(LPF);
fprintf("LPF is");
disp(N)


