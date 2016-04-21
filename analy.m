function F=driver

% clear the console screen
clc; close all;clf
% load the data structure with info pertaining to the physical problem
dat.condu=@condu;
dat.esrc=@esrc;
dat.hgap=15764;
dat.width=[0.003175 0.034823 0.036];
bc.left.type=2; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=400; % (that data is C in: +Ddu/dn=C // u/4+D/2du/dn=C // u=C)
bc.rite.type=2;
bc.rite.C=2;
dat.bc=bc; clear bc;

verif_hc_eq(dat)

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function verif_hc_eq(dat)

cd=dat.condu; src=dat.esrc; hgap=dat.hgap; L=dat.width; 

switch dat.bc.left.type
    case 0 % Neumann
        % +Ddu/dn=C on the left becomes: -Ddu/dx=C <==> -D.B=C <==> B=-C/D
        mat(1,1:6) =[1,0,0,0,0,0];
        b(1) = -dat.bc.left.C / cd(L(1));
    case 1 % Robin
        % u/4+D/2du/dn=C on the left becomes: u/4-D/2du/dx=C 
        % <==> -(D/2)*B+E/4=C 
        mat(1,1:6) =[-cd(L(1))/2,1/4,0,0,0,0];
        b(1) = dat.bc.left.C;
    case 2 % Dirichlet
        mat(1,1:6) =[0,1,0,0,0,0];
        b(1) = dat.bc.left.C;
end
switch dat.bc.rite.type
    case 0 % Neumann
        % +Ddu/dn=C on the right becomes: Ddu/dx=C <==> D.B=C <==> B=C/D
        mat(6,1:6) =[0,0,0,0,1,0];
        b(6) = dat.bc.rite.C / cd(L(3));
    case 1 % Robin
        % u/4+D/2du/dn=C on the right becomes: u/4+D/2du/dx=C 
        % <==> (BL+E)/4 + (D/2)*B = C 
        mat(6,1:6) =[0,0,0,0,L(3)/4+cd(L(3))/2,1/4];
        b(6) = dat.bc.rite.C;
    case 2 % Dirichlet
        mat(6,1:6) =[0,0,0,0,L(3),1];
        b(6) = dat.bc.rite.C;
end

Y1=-src(L(2))/cd(L(1));
Y2=-src(L(2))/(2*cd(L(2)));

mat(2,1:6) =[cd(L(1))/cd(L(2)),0,-1,0,0,0];
mat(3,1:6) =[L(1),1,-L(1),-1,0,0];
mat(4,1:6) =[0,0,2*cd(L(2))+L(2)*hgap,hgap,-L(2)*hgap,-hgap];
mat(5,1:6) =[0,0,L(2)*hgap,hgap,2*cd(L(3))-L(2)*hgap,-hgap];
b(2) =Y1*L(1);
b(3) =Y2*L(1)*L(1);
b(4) =-hgap*Y2*L(2)*L(2)+2*src(L(2))*L(2);
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
w=[0.003175 0.034823 0.036];
if x<=w(2)
    y=18;
else
    y=16;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=esrc(x)
w=[0.003175 0.034823 0.036];
if (x>w(1) | x<=w(2)) 
    y=5000000;
else
    y=0;
end
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%