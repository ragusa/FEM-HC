function F=driver

% clear the console screen
clc; close all;clf
% load the data structure with info pertaining to the physical problem
global dat
dat.condu=@condu;
dat.esrc=@esrc;
dat.hgap=15764;
dat.hconv=20000;
dat.tcool=50;
dat.width=[0.003175 0.034823 0.036];
bc.left.type=2; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=400; % (that data is C in: kdu/dn=C // kdu/dn = hconv(C-Tcool) // u=C)
bc.rite.type=2;
bc.rite.C=80;
dat.bc=bc; clear bc;

verif_hc_eq;

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function verif_hc_eq
global dat

cd=dat.condu; src=dat.esrc; hgap=dat.hgap; hconv=dat.hconv;
tcool=dat.tcool; L=dat.width; 

% general form of the solution:
% Zone 1 : T1 = B1*x + E1
% dT1/dx= B1
% Zone 2 : T2 = -q/(2*k2)*(x.^2) + B2*x + E2 = Y2*(x.^2) + B2*x + E2
% dT2/dx= -q/k2*x + B2 = Y1*x + B2
% Zone 3 : T3 = B3*x + E3
% dT3/dx= B3

% definition of 2 coefficients
Y1=-src(L(2))/cd(L(1));
Y2=-src(L(2))/(2*cd(L(2)));

switch dat.bc.left.type
    case 0 % Neumann
        % k1*du1/dn=C on the left becomes: -k1*du1/dx=C 
        % <==> -k1*B1=C <==> B1=-C/k1
        mat(1,1:6) =[1,0,0,0,0,0];
        b(1) = -dat.bc.left.C / cd(L(1));
    case 1 % Robin
        % k1*du1/dn = hconv(C-Tcool) on the left becomes: -k1*du1/dx=hconv(C-Tcool) 
        % <==> -k1*B1=hconv(C-Tcool) 
        mat(1,1:6) =[-cd(L(1)),0,0,0,0,0]; %cd(L(1))=cd(L(0)) in the function definition
        b(1) = hconv*(dat.bc.left.C-tcool);
    case 2 % Dirichlet
        % u1=C <==> E1=C
        mat(1,1:6) =[0,1,0,0,0,0];
        b(1) = dat.bc.left.C;
end
switch dat.bc.rite.type
    case 0 % Neumann
        % k3*du3/dn=C on the right becomes: k3*du3/dx=C 
        % <==> k3*B3=C <==> B3=C/k3
        mat(6,1:6) =[0,0,0,0,1,0];
        b(6) = dat.bc.rite.C / cd(L(3));
    case 1 % Robin
        % k3*du3/dn = hconv(C-Tcool) on the right becomes: k3*du3/dx=hconv(C-Tcool) 
        % <==> k3*B3=hconv(C-Tcool) 
        mat(6,1:6) =[0,0,0,0,cd(L(3)),0];
        b(6) = hconv*(dat.bc.rite.C-tcool);
    case 2 % Dirichlet
        % u3=C <==> B3*L+E3=C
        mat(6,1:6) =[0,0,0,0,L(3),1];
        b(6) = dat.bc.rite.C;
end

% fixed conditions
% continuity of T and flux between zone 1 and zone 2 (interface L1)
% T1(L1)=T2(L1) <==> B1*L1+E1=Y2*(L1^2)+B2*L1+E2
% <==> B1*L1+E1-B2*L1-E2=Y2*(L1^2)
mat(2,1:6) =[L(1),1,-L(1),-1,0,0];
b(2) =Y2*L(1)*L(1);
% phi1(L1)=phi2(L1) <==> k1*B1=k2*(Y1*L1+B2)
% <==> (k1/k2)*B1-B2=Y1*L1
mat(3,1:6) =[cd(L(1))/cd(L(2)),0,-1,0,0,0];
b(3) =Y1*L(1);

% discontinuity of T between zone 2 and zone 3 (interface L2)
% T2(L2)=Y2*(L2^2)+B2*L2+E2
% T3(L2)=B3*L2+E3
% Tg=(T2(L2)+T3(L2))/2
% -k2*dT2/dx=hgap(T2(L2)-Tg)
% <==> -k2(-q/k2*L2+B2)=hgap(T2(L2)-T3(L2))/2
% <==> (2k2+L2*hgap)B2+hgap*E2-L2*hgap*B3-hgap*E3=-hgap*Y2*(L2^2)+2*q*L2
mat(4,1:6) =[0,0,2*cd(L(2))+L(2)*hgap,hgap,-L(2)*hgap,-hgap];
b(4) =-hgap*Y2*L(2)*L(2)+2*src(L(2))*L(2);
% -k3*dT3/dx=hgap(Tg-T3(L2))
% <==> -k3*B3=hgap(T2(L2)-T3(L2))/2
% <==> L2*hgap*B2+hgap*E2+(2k3-L2*hgap)*B3-hgap*E3=-hgap*Y2*(L2^2)
mat(5,1:6) =[0,0,L(2)*hgap,hgap,2*cd(L(3))-L(2)*hgap,-hgap];
b(5) =-hgap*Y2*L(2)*L(2);

% get coefficient for the analytical solution
a=mat\b';
x1=linspace(0,L(1));
x2=linspace(L(1),L(2));
x3=linspace(L(2),L(3));
y1=a(1)*x1+a(2);
y2=Y2*(x2.^2)+a(3)*x2+a(4);
y3=a(5)*x3+a(6);

plot(x1,y1,x2,y2,x3,y3); hold all;
title('1D heat conduction problem')
xlabel('Width')
ylabel('Temperature')

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=condu(x)
global dat
if x<=dat.width(2)
    y=18;
else
    y=16;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=esrc(x)
global dat
if (x>dat.width(1) || x<=dat.width(2)) 
    y=5000000;
else
    y=0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%