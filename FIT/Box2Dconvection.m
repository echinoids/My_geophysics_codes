% Script modificado del libro Turcotte, 2014, Geodynamics, 3ra Ed,
% Cambridge University Press, pag 548)

%INPUTS:

% KAPPA = difusividad termica, m²/s
% ETA = viscosidad, Pa s
% RHO = Densidad, kg/m³
% ALPHA = coef de expansion termica, 1/K
% G = gravedad, m/s²
% TTOP = Temperatura tope, K
% TBOTTOM = Temperatura base, K 
% X = ancho de la caja, m
% Y = altura de la caja, m

%OUTPUTS:

% x_temp, y_temp : coordenadas nodos temperatura, Km
% x_vel, y_vel : coordenadas nodos velocidad, Km
% Temp = temperatura del fluido, K
% vx_out,vy_out = componentes x e y del campo de velocidad del fluido, m/s


% 2D steady state convection in a rectangular box
% for fluid with constant viscosity
% and constant thermal diffusivity
%
% Numerical solution is obtained
% based on finite-differences
% for pressure-velocity formulation
% with fully staggered grid
% (for numerical details see Gerya. T., 2010,
% Introduction to numerical geodynamic modelling,
% Cambridge University Press)

function [x_temp,y_temp,x_vel,y_vel,Temp,vx_out,vy_out] = Box2Dconvection(KAPPA,ETA,RHO,ALPHA,G,TTOP,TBOTTOM,X,Y);

NX=51; % model resolution in horizontal direction
NY=51; % model resolution in vertical direction

% Mechanical boundary conditions: -1=free slip, +1=no slip
BTOP=-1; % Boundary condition at the top
BBOTTOM=-1; % Boundary condition at the bottom
BLEFT=-1; % Boundary condition at the left wall
BRIGHT=-1; % Boundary condition at the right wall

% Convergence criterion for LSQ(Tnew-Told), K
DTMAX=0.1;

% Computing ============================================
% Grid steps, m
dx=X/(NX-1); % Horizontal
dy=Y/(NY-1); % vertical
% Number of unknowns
NTK=(NX+1)*(NY+1); % Temperature
NVP=(NX+1)*(NY+1); % Vx,Vy,Pr nodes
% Matrices for mechanical solution
L=sparse(NVP*3,NVP*3);
R=zeros(NVP*3,1);
% Matrices for thermal solution
LT=sparse(NTK,NTK);
RT=zeros(NTK,1);
% Initial temperature distribution
TK=ones((NY+1),(NX+1))*(TTOP+TBOTTOM)/2;
TK(1:3,:)=TTOP;
TK(1:fix(NY/2),1:3)=TTOP;
TK(NY-2:NY+1,:)=TBOTTOM;
TK(fix(NY/2)+1:NY+1,NX-2:NX+1)=TBOTTOM;
% Arrays for Vx,Vy,Pr
vx=zeros((NY+1),NX);
vy=zeros(NY,(NX+1));
pr=zeros((NY-1),(NX-1));
pscale=ETA/dx; % Pressure scaling

% Iteration begin
DT=DTMAX*100;
iter=0;
n_iter=500;
while(DT>DTMAX && iter<=n_iter)
    iter=iter+1;
    
    for j=1:NX+1
        for i=1:NY+1
            % Global indexes
            kx=((j-1)*(NY+1)+(i-1))*3+1;
            ky=kx+1;
            kp=kx+2;
            % Vx
            % Boundary Conditions
            if(j==1 || j==NX || j==NX+1 || i==1 || i==NY+1)
                % Ghost Vx nodes and left and right walls: Vx=0
                if(j==1 || j==NX || j==NX+1)
                L(kx,kx)=1;
                R(kx)=0;
                end
                % Top boundary
                if(i==1 && j>1 && j<NX)
                L(kx,kx)=1;
                L(kx,kx+3)=BTOP;
                R(kx)=0;
                end
                % Bottom boundary
                if(i==NY+1 && j>1 && j<NX)
                L(kx,kx)=1;
                L(kx,kx-3)=BBOTTOM;
                R(kx)=0;
                end
            else
                % Composing x-Stokes equation:
                % ETA*(d2Vx/dx^2+d2Vx/dy^2)-dP/dx=0
                L(kx,kx)=-2*ETA/dx/dx-2*ETA/dy/dy;
                L(kx,kx-(NY+1)*3)=ETA/dx/dx;
                L(kx,kx+(NY+1)*3)=ETA/dx/dx;
                L(kx,kx-3)=ETA/dy/dy;
                L(kx,kx+3)=ETA/dy/dy;
                L(kx,kp)=pscale/dx;
                L(kx,kp+(NY+1)*3)=-pscale/dx;
                R(kx)=0;
            end

            % Vy
            % Boundary Conditions
            if(j==1 || j==NX+1 || i==1 || i==NY || i==NY+1)
                % Ghost Vy nodes and top and Bottom: Vy=0
                if(i==1 || i==NY || i==NY+1)
                L(ky,ky)=1;
                R(ky)=0;
                end
                % Left boundary
                if(j==1 && i>1 && i<NY)
                L(ky,ky)=1;
                L(ky,ky+(NY+1)*3)=BLEFT;
                R(ky)=0;
                end
                % Right boundary
                if(j==NX+1 && i>1 && i<NY)
                L(ky,ky)=1;
                L(ky,ky-(NY+1)*3)=BRIGHT;
                R(ky)=0;
                end
            else
                % Composing x-Stokes equation:
                % ETA*(d2Vy/dx^2+d2Vy/dy^2)-dP/dx=-RHO*(1-ALPHA*(T-TTOP))*G
                L(ky,ky)=-2*ETA/dx/dx-2*ETA/dy/dy;
                L(ky,ky-(NY+1)*3)=ETA/dx/dx;
                L(ky,ky+(NY+1)*3)=ETA/dx/dx;
                L(ky,ky-3)=ETA/dy/dy;
                L(ky,ky+3)=ETA/dy/dy;
                L(ky,kp)=pscale/dy;
                L(ky,kp+3)=-pscale/dy;
                R(ky)=-RHO*(1-ALPHA*((TK(i,j)+TK(i+1,j))/2-TTOP))*G;
            end

            % Pr
            % Boundary Conditions
            if(j==1 || j==NX+1 || i==1 || i==NY+1 || (j==2 && i==2))
                L(kp,kp)=1;
                R(kp)=0;
                % Upper left corner pressure
                if(j==2 && i==2)
                    L(kp,kp)=1*pscale;
                    R(kp)=RHO*G*dy/2;
                end
           else
                % Composing continuity equation:
                % dVx/dx+dVy/dy=0
                L(kp,kx-(NY+1)*3)=-1/dx;
                L(kp,kx)=1/dx;
                L(kp,ky-3)=-1/dy;
                L(kp,ky)=1/dy;
                R(kp)=0;
            end
        end
    end
    
    % Solving matrix
    S=L\R;
    % Reloading solution for Vx
    for j=1:NX+1
        for i=1:NY+1
            % Global indexes
            kx=((j-1)*(NY+1)+(i-1))*3+1;
            ky=kx+1;
            kp=kx+2;
            if(j<NX+1)
                vx(i,j)=S(kx);
            end
            if(i<NY+1)
                vy(i,j)=S(ky);
            end
            if(j>1 && i>1 && j<NX+1 && i<NY+1)
                pr(i-1,j-1)=S(kp);
            end
        end
    end
    
    % Thermal solution
    % Composing matrix
    for j=1:NX+1
        for i=1:NY+1
            % Global index
            kt=(j-1)*(NY+1)+i;
            % Boundary Conditions
            if(j==1 || j==NX+1 || i==1 || i==NY+1)
                % Top boundary: T=TTOP
                if(i==1)
                    LT(kt,kt)=0.5;
                    LT(kt,kt+1)=0.5;
                    RT(kt)=TTOP;
                end
                % Bottom boundary
                if(i==NY+1)
                    LT(kt,kt)=0.5;
                    LT(kt,kt-1)=0.5;
                    RT(kt)=TBOTTOM;
                end
                % Left boundary
                if(j==1 && i>1 && i<NY+1)
                    LT(kt,kt)=1;
                    LT(kt,kt+(NY+1))=-1;
                    RT(kt)=0;
                end
                % Right boundary
                if(j==NX+1 && i>1 && i<NY+1)
                    LT(kt,kt)=1;
                    LT(kt,kt-(NY+1))=-1;
                    RT(kt)=0;
                end
            else
                % Composing Eulerian steady heat conservation equation:
                % KAPPA*(d2T/dx^2+d2Vx/dy^2)-vx*dT/dx-vy*dT/dy=0
                LT(kt,kt)=-2*KAPPA/dx/dx-2*KAPPA/dy/dy;
                LT(kt,kt-(NY+1))=KAPPA/dx/dx;
                LT(kt,kt+(NY+1))=KAPPA/dx/dx;
                LT(kt,kt-1)=KAPPA/dy/dy;
                LT(kt,kt+1)=KAPPA/dy/dy;
                % Vx Velocity in the temperature node
                vxcur=(vx(i,j)+vx(i,j-1))/2;
                LT(kt,kt-(NY+1))=LT(kt,kt-(NY+1))+vxcur/dx/2;
                LT(kt,kt+(NY+1))=LT(kt,kt+(NY+1))-vxcur/dx/2;
                % Vy Velocity in the temperature node
                vycur=(vy(i,j)+vy(i-1,j))/2;
                LT(kt,kt-1)=LT(kt,kt-1)+vycur/dy/2;
                LT(kt,kt+1)=LT(kt,kt+1)-vycur/dy/2;
                RT(kt)=0;
            end
        end
    end
    % Solving matrix
    ST=LT\RT;
    % Reloading solution for TK1
    TK1=TK;
    DT=0;
    for j=1:NX+1
        for i=1:NY+1
            % Global index
            kt=(j-1)*(NY+1)+i;
            TK1(i,j)=ST(kt);
            DT=DT+(TK1(i,j)-TK(i,j)).^2;
        end
    end
    
    % Computing LSQ temperature change
    DT=(DT/NTK)^0.5;
    % Computing ============================================
    % Nusselt number
    NUSSELT=0;
    for j=2:NX
        NUSSELT=NUSSELT+Y*(TK1(2,j)-TK1(1,j))/dy*dx/X/(TBOTTOM-TTOP);
    end
    disp(['(K) for iteration = ',num2str(iter),...
    '  LSQ(Tnew-Told) = ',num2str(DT)])   
    
    % T-Cells coordinates
    xcel=-dx/2:dx:X+dx/2;
    ycel=-dy/2:dy:Y+dy/2;
    
    TK=TK1;
   
end

x_temp = xcel/1000;
y_temp = ycel/1000;
x_vel = xcel(2:2:end-1)/1000;
y_vel = ycel(2:2:end-1)/1000;
Temp = TK1;
vx_out = vx(2:2:end-1,2:2:end);
vy_out = vy(2:2:end-1,2:2:end-1);

  