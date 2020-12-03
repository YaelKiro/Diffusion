% Crank-Nicolson method for calculating diffusion equation for diffusion of
% He in diamonds together with build up of 4He from radioactive decay of U and Th
function [C,C_tot] = FDM_PDE_Diff_Diamonds(Inp)
%Inp is an array of all the input parameters that includes boundary and initial
%conditions and numerical parameters.
%Output is the concentration of He as a function of distance from the
%diamond's center (r) and time (C - 2D matrix).
%C_tot is the total concentration in diamond vs. time (1D matrix)

%Parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% UNITS HERE ARE A SUGGESTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
R=Inp(1);               %radius of diamond in cm
time=Inp(2);            %time in G years
D=Inp(3);               %Diffusion coefficient in cm2/s
D=D*3600*24*365.26*1e9; % transfer to units D - cm2/Gy
U8=Inp(4);              %238U mol/g
k=Inp(5);               %Th/U molar ratio
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initial and boundary conditions %%%%%%%%%%%%%%%%%%%%%%%%%%%
C0=Inp(6); % initial He concentration in diamond
Cm=Inp(7); % boundary conditions, He concentration in mantle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Numerical parameters for discretization
td=Inp(8); %temporal discretization - number of time steps
N=Inp(9);  %spatial discretization - number of nodes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

production=Inp(10); % if He (4He) is produced by radioactive decay set - 1, for 3He set - 0
%%%%%%%%%% DECAY CONSTANTS %%%%%%%%%%%%%%%%%%
l8=1.5513e-1;  %238U decay constant 1/Gy
l5=9.8485e-1;  %235U decay constant 1/Gy
l2=0.49475e-1; %232Th decay constant 1/Gy
a1=U8*8*l8;
a2=7/137.88*U8*l5;
a3=U8*k*6*l2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%FINITE DIFFERENCE METHOD - CRANK-NICOLSON %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

dr=R/N;             % Size of element (cm)
dt=time/td;         % Size of time step (Gy)
r=[dr:dr:(R-dr)];   % spatial array

% Building spatial matrices
%%%% SEE MATRICES IN METHOD SECTION AT WEISS ET AL. 2021, NATURE COMMUNICATIONS
AA=zeros(N-1,N-1);

% For r=0 (r(1))%%%%%%%%%%
rr=r(1);
rad=1/rr;                % sphere coordinates term
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
%%%%%%%%%%%%%%%%%%%%

% Spatial loop
for i=2:(N-2)
    rr=r(i);
    rad=1/rr;            % sphere coordinates term
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

% For R=r
rr=r(N-1);
rad=1/rr;               % sphere coordinates term
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

% Dealing with boundary and initial conditions
CN=zeros(N-1,1);
CN(N-1,1)=2*Cm*A3; % vector dealing with boundary conditions at r=R, see detail in method section at Weiss et al., 2021, NCM

C(1:(N-1),1)=C0;   % initial conditions


t=0:time/td:time;
f=0; % production of He term

invAA=inv(AA);

K1=invAA*BB;

K2=invAA*ones(N-1,1);

K3=invAA*CN;

k=td+1;
% Temporal loop
for j=2:(td+1)
    k=k-1;
    if production==1;
        f=(a1*exp(l8*(time-t(j)))+a2*exp(l5*(time-t(j)))+a3*exp(l2*(time-t(j))))*dr/D;
    end
    
    C(:,j)=K1*C(:,j-1)+K2*f+K3; % Concentration of He with distance from diamond's center and time (2D matrix)
    
    C_tot(j)=3*sum(C(:,j).*r'.^2)/R^3*dr+3*(C(end,j)+Cm)/2*dr/R;  % Total concentration of helium in diamond with time (1D matrix)
    
end
% Adding initial conditions to final results    
C_tot(1)=C0;
C=[C(1,:);C;ones(1,td+1).*Cm];
