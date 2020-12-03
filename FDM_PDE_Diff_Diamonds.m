% Crank-Nicolson method for calculating diffusion equation for diffusion of
% He in diamonds together with build up of 4He
function [C,C_tot] = FDM_PDE_Diff_Diamonds(Inp)
%Inp is an array all the input parameters that includes boundary and initial 
%conditions and numerical parameters.

%Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=Inp(1); %radius of diamond in cm
time=Inp(2); %time in G years
D=Inp(3); %Diffusion in cm2/s
U8=Inp(4); %ppm
k=Inp(5); %Th/U ratio (mole?)

% Initial and boundary conditions %%%%%%%%%%%%%%%%

C0=Inp(6);  %%% initial He concentration in diamond
Cp=Inp(6); % present He concentrations for bw calc
Cm=Inp(7); %%% boundary conditions 0 He concentration in the mantle

% Numerical parameters for discretization;

td=Inp(8); %temporal discretization
N=Inp(9); %spatial discretization

production=Inp(10);

fw=Inp(11); % 1 for fw model 0 for bw model

D=D*3600*24*365.26*1e9; % consistancy in units D - cm2/Gy

l8=1.5513e-1; %238U decay constant 1/Gy
l5=9.8485e-1; %235U decay constant 1/Gy
l2=0.49475e-1; %232Th decay constant 1/Gy
a1=U8*8*l8;
a2=7/137.88*U8*l5;
a3=U8*k*6*l2;
%FINITE DIFFERENCE METHOD%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Spatial matrix%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dr=R/N; 
dt=time/td;
r=[dr:dr:(R-dr)];

AA=zeros(N-1,N-1);
rr=r(1);
rad=1/rr; % sphere coordinates term
%rad=0;
A1=1/2*(rad-1/dr);
A2=dr/(D*dt)+1/dr;
A3=-1/2*(rad+1/dr);
AA(1,1)=A1+A2;
AA(1,2)=A3;

BB=zeros(N-1,N-1);
B1=-A1;
B2=dr/(D*dt)-1/dr;
B3=-A3;
BB(1,1)=B1+B2;
BB(1,2)=B3;

for i=2:(N-2)
    rr=r(i);
    rad=1/rr; % sphere coordinates term
    A1=1/2*(rad-1/dr);
    A2=dr/(D*dt)+1/dr;
    A3=-1/2*(rad+1/dr);

    
    B1=-A1;
    B2=dr/(D*dt)-1/dr;
    B3=-A3;


    AA(i,i-1)=A1;
    AA(i,i)=A2;
    AA(i,i+1)=A3;
    
    BB(i,i-1)=B1;
    BB(i,i)=B2;
    BB(i,i+1)=B3;
    
end

rr=r(N-1);
rad=1/rr; % sphere coordinates term
A1=1/2*(rad-1/dr);
A2=dr/(D*dt)+1/dr;
A3=-1/2*(rad+1/dr);
AA(N-1,N-2)=A1;
AA(N-1,N-1)=A2;

B1=-A1;
B2=dr/(D*dt)-1/dr;
B3=-A3;
BB(N-1,N-2)=B1;
BB(N-1,N-1)=B2;

CN=zeros(N-1,1);
CN(N-1,1)=2*Cm*A3; % vector dealing with boundary conditions at r=R (check this)
C(1:(N-1),1)=C0;  % initial conditions

Cbw(1:(N-1),td+1)=Cp;
t=0:time/td:time;
f=0;
if fw==1 
    
    invAA=inv(AA);
    
    K1=invAA*BB;
    
    K2=invAA*ones(N-1,1);
    
    K3=invAA*CN;
else
        % Matrixes for backward calc
    invBB=inv(BB);
    K1b=invBB*AA;
    K2b=invBB*ones(N-1,1);
    K3b=invBB*CN;
end
    k=td+1;
for j=2:(td+1)
        k=k-1;
        if fw==1
            if production==1;
                f=(a1*exp(l8*(time-t(j)))+a2*exp(l5*(time-t(j)))+a3*exp(l2*(time-t(j))))*dr/D;
            end
            % forward model based on initial conditions T years ago
            C(:,j)=K1*C(:,j-1)+K2*f+K3;
          
            C_tot(j)=3*sum(C(:,j).*r'.^2)/R^3*dr+3*(C(end,j)+Cm)/2*dr/R;  %check this
        else
            % backward model based on present conditions
            if production==1;
                f=(a1*exp(l8*(time-t(k)))+a2*exp(l5*(time-t(k)))+a3*exp(l2*(time-t(k))))*dr/D;
               
            end
            Cbw(:,k)=K1b*Cbw(:,k+1)-K2b*f-K3b;
            Cbw_tot(k)=3*sum(Cbw(:,k).*r'.^2)/R^3*dr+3*Cm*dr/R;
        end
    end
    if fw==1
        C_tot(1)=C0;
        C=[C(1,:);C;ones(1,td+1).*Cm];
    else
        Cbw_tot(td+1)=Cp;
        C=[Cbw(1,:);Cbw;ones(1,td+1).*Cm];
        C_tot=Cbw_tot;
    end
    
    
        
