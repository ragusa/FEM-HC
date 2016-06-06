function F=fem3
% clean heat conduction fem code with analytical solution

% clear the console screen
clc; close all;clf
% load the data structure with info pertaining to the physical problem
dat.condu=@condu;
dat.esrc=@esrc;
dat.width=10;
bc.left.type=2; %0=neumann, 1=robin, 2=dirichlet
bc.left.C=7; % (that data is C in: +Ddu/dn=C // u/4+D/2du/dn=C // u=C)
bc.rite.type=2;
bc.rite.C=2;
dat.bc=bc; clear bc;

% load the numerical parameters, npar, structure pertaining to numerics
% number of elements
npar.nel = 20;
% domain
npar.x = linspace(0,dat.width,npar.nel+1);
% polynomial degree
npar.porder=1;
% nbr of dofs per variable
npar.ndofs = npar.porder*npar.nel+1;
% connectivity
gn=zeros(npar.nel,npar.porder+1);
gn(1,:)=linspace(1,npar.porder+1,npar.porder+1);
for iel=2:npar.nel
    gn(iel,:)=[gn(iel-1,end) , gn(iel-1,2:end)+npar.porder ];
end
npar.gn=gn; clear gn;

% solve system
F = solve_fem3(dat,npar);
% plot
figure(1)

% % verification is always good
a=verif_hc_eq(dat);

cd=dat.condu; src=dat.esrc; L=dat.width;

x=linspace(0,L);
y=-src(1)/(2*cd(1))*(x.^2)+a(1)*x+a(2);
plot(npar.x,F,'.-',x,y,'r-');hold all
title('1D heat conduction problem, 1 zone, Cartesian coordinates')
xlabel('Width')
ylabel('Temperature')
legend('FEM','Analytical','Location','northoutside','Orientation','horizontal')

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function u=solve_fem3(dat,npar)

% initial guess
% (not needed if a direct method is used to solve the linear system)
u=ones(npar.ndofs,1);

% assemble the matrix and the rhs
[A,b]=assemble_system(npar,dat);

% solve
u=A\b;

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [A,rhs]=assemble_system(npar,dat)
% assemble the matrix, the rhs, apply BC

% shortcuts
porder= npar.porder;
gn    = npar.gn;
nel   = npar.nel;
% ideally, we would analyze the connectivity to determine nnz
nnz=(porder+1)*nel; %this is an upperbound, not exact
% n: linear system size
n=nel*porder+1;
% allocate memory
A=spalloc(n,n,nnz);
rhs=zeros(n,1);

% compute local matrices
% load Gauss Legendre quadrature (GLQ is exact on polynomials of degree up to 2n-1,
% while using only integrand evaluations (n-point quadrature).
% estimate the max poly order (it comes from the mass matrix  when coef. are
% piecewise constant and the material mesh matches the computational mesh
poly_max=2*porder;
[xq,wq] = GLNodeWt(porder+1);
% initialize local matrices/vectors
m=zeros(porder+1,porder+1);
k=m;
f=zeros(porder+1,1);
% store shapeset
[b,dbdx] =feshpln(xq,porder);

% definition of the weak form:
% int_domain (grad u D grad b) - int_bd_domain (b D grad u n) ...
%      + int_domain( XSa u b) = int_domain (b rhs)

% loop over elements
for iel=1:npar.nel
    % element extremities
    x0=npar.x(iel);
    x1=npar.x(iel+1);
    % jacobian of the transformation to the ref. element
    Jac=(x1-x0)/2;
    % get x values in the interval
    x=(x1+x0)/2+xq*(x1-x0)/2;
    d=dat.condu(x);
    q=dat.esrc(x);
    % compute local matrices + load vector
    for i=1:porder+1
        for j=1:porder+1
            m(i,j)= dot(0.*wq.*b(:,i)    , b(:,j));
            k(i,j)= dot(d.*wq.*dbdx(:,i) , dbdx(:,j));
        end
        f(i)= dot(q.*wq, b(:,i));
    end
    % assemble
    A(gn(iel,:),gn(iel,:)) = A(gn(iel,:),gn(iel,:)) + ...
       m*Jac + k/Jac;
    rhs(gn(iel,:)) = rhs(gn(iel,:)) + f*Jac;
end

% apply natural BC
Dirichlet_nodes=[];
Dirichlet_val=[];
switch dat.bc.left.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(1)=rhs(1)+dat.bc.left.C;
    case 1 % Robin
        A(1,1)=A(1,1)+1/2;
        rhs(1)=rhs(1)+2*dat.bc.left.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes 1];
        Dirichlet_val=[Dirichlet_val dat.bc.left.C];
end
switch dat.bc.rite.type
    case 0 % Neumann, int_bd_domain (b D grad u n) is on the RHS
        rhs(n)=rhs(n)+dat.bc.rite.C;
    case 1 % Robin
        A(n,n)=A(n,n)+1/2;
        rhs(n)=rhs(n)+2*dat.bc.rite.C;
    case 2 % Dirichlet
        Dirichlet_nodes=[Dirichlet_nodes n];
        Dirichlet_val=[Dirichlet_val dat.bc.rite.C];
end
% apply Dirichlet BC
for i=1:length(Dirichlet_nodes);% loop on the number of constraints
    id=Dirichlet_nodes(i);      % extract the dof of a constraint
    bcval=Dirichlet_val(i);
    rhs=rhs-bcval*A(:,id);  % modify the rhs using constrained value
    A(id,:)=0; % set all the id-th row to zero
    A(:,id)=0; % set all the id-th column to zero (symmetrize A)
    A(id,id)=1;            % set the id-th diagonal to unity
    rhs(id)=bcval;         % put the constrained value in the rhs
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [shapefun,dhdx]=feshpln (xv,p)
% computes the finite element basis functions and their first derivatives
% for any order p
% here, we use Lagrangian FE functions

xd=linspace(-1,1,p+1);

shapefun=zeros(length(xv),p+1);
dhdx    =zeros(length(xv),p+1);

% shape function
for i=1:p+1
    num=1.;
    den=1.;
    for j=1:p+1
        if(j~=i)
            num=num.*(xv-xd(j));
            den=den.*(xd(i)-xd(j));
        end
    end
    shapefun(:,i)=num./den;
end

% derivative of the shape function
for i=1:p+1
    sum=0.;
    den=1.;
    for j=1:p+1
        if(j~=i)
            num=1;
            for k=1:p+1
                if((k~=i)&&(k~=j))
                    num=num.*(xv-xd(k));
                end
            end
            sum=sum+num;
            den=den.*(xd(i)-xd(j));
        end
    end
    dhdx(:,i)=sum./den;
end

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [x,w] = GLNodeWt(n)
% GLNodeWt  Nodes and weights for Gauss-Legendre quadrature of arbitrary order
%           obtained by solving an eigenvalue problem
%
% Synopsis:  [x,w] = GLNodeWt(n)
%
% Input:     n = order of quadrature rule
%
% Output:    x = vector of nodes
%            w = vector of weights

%  Algorithm based on ideas from Golub and Welsch, and Gautschi.  For a
%  condensed presentation see H.R. Schwarz, "Numerical Analysis: A
%  Comprehensive Introduction," 1989, Wiley.  Original MATLAB
%  implementation by H.W. Wilson and L.H. Turcotte, "Advanced Mathematics
%  and Mechanics Applications Using MATLAB," 2nd ed., 1998, CRC Press

beta   = (1:n-1)./sqrt(4*(1:n-1).^2 - 1);
J      = diag(beta,-1) + diag(beta,1);    % eig(J) needs J in full storage
[V,D]  = eig(J);
[x,ix] = sort(diag(D));  %  nodes are eigenvalues, which are on diagonal of D
w      = 2*V(1,ix)'.^2;  %  V(1,ix)' is column vector of first row of sorted V

return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function a=verif_hc_eq(dat)

cd=dat.condu; src=dat.esrc; L=dat.width;

% general form of the solution:
% T = Y x x + B x + E
% dT/dx= 2Yx + B
Y=-src(1)/(2*cd(1));
switch dat.bc.left.type
    case 0 % Neumann
        % +Ddu/dn=C on the left becomes: -Ddu/dx=C <==> -D.B=C <==> B=-C/D
        mat(1,1:2) =[1,0];
        b(1) = -dat.bc.left.C / cd(0);
    case 1 % Robin
        % u/4+D/2du/dn=C on the left becomes: u/4-D/2du/dx=C 
        % <==> E/4-D/2(B)=C 
        mat(1,1:2) =[-cd(0)/2,1/4];
        b(1) = dat.bc.left.C;
    case 2 % Dirichlet
        mat(1,1:2) =[0,1];
        b(1) = dat.bc.left.C;
end
switch dat.bc.rite.type
    case 0 % Neumann
        % +Ddu/dn=C on the right becomes: Ddu/dx=C <==> D.(2YL+B)=C <==> B=C/D-2YL
        mat(2,1:2) =[1,0];
        b(2) = dat.bc.rite.C / cd(L) - 2*Y*L;
    case 1 % Robin
        % u/4+D/2du/dn=C on the right becomes: u/4+D/2du/dx=C 
        % <==> (YL^2+BL+E)/4 + D/2(2YL+B) = C 
        mat(2,1:2) =[L/4+cd(L)/2,1/4];
        b(2) = dat.bc.rite.C - Y*L*L/4 -cd(L)*Y*L;
    case 2 % Dirichlet
        mat(2,1:2) =[L,1];
        b(2) = dat.bc.rite.C - Y*L*L;
end
% get coefficient for the analytical solution
a=mat\b';

return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=condu(x)
y=.5;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function y=esrc(x)
y=1;
return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%