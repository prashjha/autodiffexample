function [x,xn,xm,w,wn] = GaussLegendreNDGauss(n,lb,ub)

% This function determines the abscisas (x) and weights (w) for the
% Gauss-Hermite quadrature of order n>1, on the interval [-INF, +INF].
    % This function is valid for any degree n>=2, as the companion matrix
    % (of the n'th degree Hermite polynomial) is constructed as a
    % symmetrical matrix, guaranteeing that all the eigenvalues (roots)
    % will be real.

% Input:
    % n = Number of quadrature points, i.e. degree of quadrature.
    % lb = LB of Uniform  distribution. 
    % ub = UB of Uniform  distribution. 
    
% Output:
    % x = Quadrature abscisas for standard normal distribution.
    % xn = Quadrature abscisas in ndgrid format. Each dimension is in a
    % separate cell.
    % xm = Quadrature abscisas for each dimension. The matrix will be n x length(lb).
    % w = 1-D quadrature weights.
    % wn = Normalized N-D quadrature weights.
    
% Error Handling
if  ~isvector(lb) || ~isvector(ub) ||  length(lb)~=length(ub)
    ~sum(lb<ub)
    disp('Input error.');
    return;
end
if size(lb,1)~=1
    lb=lb';
end
if size(ub,1)~=1
    ub=ub';
end

% 1-D Gauss-Legendre Quadrature Points and Weights Solve
% https://keisan.casio.com/exec/system/1280624821
switch (n)
   case(1)
    x = [ 0];
    w = [ 2];
   case(2)
    x = [ -0.5773502691896257645092; 0.5773502691896257645092];
    w = [ 1;1];
   case(3)
    x = [ -0.7745966692414833770359; 0; 0.7745966692414833770359];
    w = [  0.5555555555555555555556; 0.8888888888888888888889; 0.555555555555555555556];
   case(4)
    x = [ -0.861136311594052575224;-0.3399810435848562648027; 0.3399810435848562648027; 0.861136311594052575224];
    w = [ 0.3478548451374538573731; 0.6521451548625461426269;	0.6521451548625461426269;0.3478548451374538573731];
   case(5)
    x = [ -0.9061798459386639927976; -0.5384693101056830910363; 0; 0.5384693101056830910363; 0.9061798459386639927976];
    w = [ 0.2369268850561890875143; 0.4786286704993664680413; 0.5688888888888888888889; 0.4786286704993664680413; 0.2369268850561890875143];
end

% Generate N-D Quadrature Point Grid 
ndim=length(lb);
xn=cell(1,ndim);
xn{1}=.5*repmat(ub-lb,n,1).*x+.5*repmat(lb,n,1)+.5*repmat(ub,n,1);
xm=xn{1};
xn=mat2cell(xn{1},n,ones(1,ndim));
[xn{:}]=ndgrid(xn{:});

% Generate N-D Weighting Grid for Gaussian Function
wn=w;
for iii=1:ndim-1
    wn=kron(wn,w);
end
if ndim>1
    wn=reshape(wn,repmat(n,1,ndim));
end
wn=wn/2^ndim;
