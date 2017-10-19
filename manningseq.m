function [Q,h,err] = manningseq(V,S,n,d,z,plt)
 
% MANNINGSEQ Solve Mannings for discharge given an irregular cross section
%
% Syntax
%
%     [Q,h,err] = manningseq(V,S,n,d,z)
%     [Q,h,err] = manningseq(V,S,n,d,z,plt)
%
% Description
%
%     manningseq solves for discharge and maximum flow depth given flow
%     velocity V, channel bed gradient S, Manning's n, and an irregular
%     river cross-section defined by the vertices in the vector d=distance
%     and z=elevation. The function finds the zero of the function
%
%        Q = 1/n * (A/P)^(2/3) * S^(1/2) * A
%
%     where A is the channel cross-sectional area and P is the wetted
%     perimeter. A/P is the hydraulic radius R.
%
%     The function requires fminsearchbnd by John D'Errico available on the
%     File Exchange (http://www.mathworks.com/matlabcentral/fileexchange/8277)
%
% Input arguments
%
%     V        flow velocity [m/s]
%     S        slope [m/m]
%     n        Manning's n [s/m^(1/3)]
%     d        horizontal distance vector of channel profile [m]
%     z        elevation vector of channel profile [m]
%     plt      plot results (true or false, default = true)
% 
% Output arguments
%
%     Q        discharge [m^3/s]
%     h        maximum depth
%     err      deviation between Q/A and V. Large errors indicate that
%              maximum flow depth may exceed the profile's maximum height
%
% Example
%     
%     d = [0 10 20 30 40 50 60 70];
%     z =[20 10 5  3  4  8  10 20];
%     [Q,h,err] = manningseq(5,0.01,0.03,d(:),z(:));
%
%
% See also: fminsearchbnd 
% 
% Author: Wolfgang Schwanghart (w.schwanghart[at]geo.uni-potsdam.de)
% Date: 18. December, 2015
 
if nargin == 5;
    plt = true;
end
 
% force column vectors
d = d(:);
z = z(:);
z = z-min(z);
 
h0 = 10;
h = fminsearchbnd(@(h) mann(h),h0,1,min(z([1 end])));
% h = fminsearch(@(h) mann(h),h0);
[A,U] = getprof(d,z,h);
 
% check flow velocity
Vhat  = 1/n*(A/U)^(2/3)*S^.5;
err   = Vhat-V;
 
Q = A*V;
 
if plt
    plot(d,z,'k-*');
    [~,~,dnew,znew] = getprof(d,z,h);
    hold on
    p = patch(dnew,znew,zeros(size(dnew)));
    set(p,'FaceColor',[0.3 0.3 1])
    hold off
end
 
function dV = mann(h)
 
[A,U] = getprof(d,z,h);
R     = hydradius(A,U);
Vhat  = 1/n *R.^(2/3) * S.^(.5);
dV    = abs(Vhat-V);
end
end
 
function R = hydradius(A,U)
 
% R = A/U
 
R = A./U;
end
 
function [A,U,dnew,znew] = getprof(d,z,h)
 
[xint,yint,seg1] = wpolyxpoly(d,z,[min(d)-1 max(d)+1],[h h]);
 
dnew = [xint(1); [d((seg1(1)+1) : seg1(2))]; xint(2)];
znew = [yint(1); [z(seg1(1)+1 : seg1(2))]; yint(2)];
 
A = polyarea(dnew,znew);
U = max(getdistance(dnew,znew));
end
 
 
 
function varargout = wpolyxpoly(varargin)
 
% return intersections of two polylines
%
% [xint,yint] = wpolyxpoly(x1,y1,x2,y2)
% [xint,yint,seg1,seg2] = wpolyxpoly(x1,y1,x2,y2)
% [...] = wpolyxpoly(x1,y1,x2,y2,maxdist)
% [...] = wpolyxpoly(x1,y1,x2,y2,maxdist,tol)
%
% Wolfgang Schwanghart
% w.schwanghart@unibas.ch (28. January 2008)
 
 
% _______________________________
% check number of input arguments
if nargin==4;
    [x1,y1,x2,y2] = varargin{:};
    flagmaxdist=0;
    tol = 1e-7;
elseif nargin==5;
    [x1,y1,x2,y2,maxdist] = varargin{:};
    flagmaxdist=1;
    tol = 1e-7;
    if ~isscalar(maxdist)
        error('maxdist must be a scalar')
    end
elseif nargin==6;
    [x1,y1,x2,y2,maxdist,tol] = varargin{:};
    flagmaxdist=1;
    if ~isscalar(maxdist)
        error('maxdist must be a scalar')
    end    
    if ~isscalar(tol)
        error('tol must be a scalar')
    end
else
    error('wrong number of input arguments')
end
 
% check input arguments
if size(x1) ~= size(y1);
    error('x1 and y1 must have same size')
elseif size(x2) ~= size(y2);
    error('x2 and y2 must have same size')
end
 
% force column vectors
x1 = x1(:);
y1 = y1(:);
x2 = x2(:);
y2 = y2(:);
 
% number of segments in polyline 1 and 2
nrseg1 = length(x1)-1;
nrseg2 = length(x2)-1;
 
% start of each segment vector
stvecx1 = x1(1:end-1);
stvecy1 = y1(1:end-1);
 
stvecx2 = x2(1:end-1);
stvecy2 = y2(1:end-1);
 
% create b
srcx = bsxfun(@minus,stvecx2,stvecx1');
srcy = bsxfun(@minus,stvecy2,stvecy1');
 
srcx = srcx(:);
srcy = srcy(:);
 
% create Index for segments
IXseg1 = repmat(1:nrseg1,nrseg2,1);
IXseg1 = IXseg1(:);
IXseg2 = repmat((1:nrseg2)',1,nrseg1);
IXseg2 = IXseg2(:);
 
% exclude segments with starting points with a distance
% farther than maxdist
if flagmaxdist    
    dis    = hypot(srcx,srcy);    
    i      = logical(dis<maxdist);            
    srcx   = srcx(i);
    srcy   = srcy(i);
    IXseg1 = IXseg1(i);
    IXseg2 = IXseg2(i);
    
    if isempty(srcx)
        varargout{1}=[];
        varargout{2}=[];
        varargout{3}=[];
        varargout{4}=[];
        return
    end     
end
 
% list end of each segment vector
endvecx1 = x1(2:end)-stvecx1;
endvecy1 = y1(2:end)-stvecy1;
 
endvecx2 = x2(2:end)-stvecx2;
endvecy2 = y2(2:end)-stvecy2;
 
% create equation matrix
x1enh = repmat(endvecx1,1,nrseg2)';
x1enh = x1enh(:);
y1enh = repmat(endvecy1,1,nrseg2)';
y1enh = y1enh(:);  
 
x2enh = repmat(endvecx2,nrseg1,1);
y2enh = repmat(endvecy2,nrseg1,1);
 
% remove values
if flagmaxdist
    x1enh = x1enh(i);
    y1enh = y1enh(i);
    x2enh = x2enh(i);
    y2enh = y2enh(i);
    clear i
end
 
% find parallel segments using cross products
% --> determination of the area of a parallelogram
c  = x1enh.*y2enh - y1enh.*x2enh;
i3 = abs(c)<tol;
clear c
if sum(i3) ~= 0;
    x1enh  = x1enh(~i3);
    y1enh  = y1enh(~i3);
    x2enh  = x2enh(~i3);
    y2enh  = y2enh(~i3);
    IXseg1 = IXseg1(~i3);
    IXseg2 = IXseg2(~i3);
    srcx   = srcx(~i3);
    srcy   = srcy(~i3);
    clear i3 
else
    clear i3;
end
 
% create sparse block-diagonal matrix A
 
nrintmax=length(srcx);
 
b     = reshape([srcx srcy]',nrintmax*2,1);
col0  = reshape([x1enh y2enh]',nrintmax*2,1);
colm1 = reshape([y1enh zeros(nrintmax,1)]',nrintmax*2,1);
colp1 = reshape([zeros(nrintmax,1) x2enh ]',nrintmax*2,1);
 
clear x1enh y1enh x2enh y2enh
     
A = spdiags([colm1 col0 colp1],[-1 0 1],nrintmax*2,nrintmax*2);
 
% solve the set of equations Ax = b
x = A\b;
clear A b;
 
% reshape x
x = reshape(x,2,nrintmax)';
 
% find rows where alpha and beta are both between 0 and 1;
i2 = (x(:,1)>=0 & x(:,1)<1) & (x(:,2)<=0 & x(:,2)>-1);
 
% check last segments
% if any(x(end,:) == 1)
%     i2(end) = true;
% end
 
seg1=IXseg1(i2);
seg2=IXseg2(i2);
 
if isempty(seg1)
    varargout{1}=[];
    varargout{2}=[];
    varargout{3}=[];
    varargout{4}=[];
    varargout{5}=[];
else
    alpha = x(i2);
    int = [x1(seg1) y1(seg1)] + bsxfun(@times,[endvecx1(seg1) endvecy1(seg1)],alpha);
    
    varargout{1} = int(:,1);
    varargout{2} = int(:,2);
    varargout{3} = seg1;
    varargout{4} = seg2; 
    varargout{5} = x(i2);
end
end

function [cumdxy] = getdistance(x,y)
    x = x(:);
    y = y(:);
    dx = diff(x);
    dy = diff(y);
    dxy = hypot(dx, dy); % square root of sum of squares
    cumdxy = [0; cumsum(dxy,1)];
end % 

