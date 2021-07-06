function basisobj = create_FEM_basis(p, e, t, order, dl)
%  CREATE_FEM_BASIS sets up a finite element basis for the analysis
%  of spatial data.  It requires a triangular mesh as input, such as is
%  set up by functions initmesh and refinemesh in the pde toolbox. 
%  See page 5-45 in the manual for the pde toolbox for an example.  
%  The triangular mesh is defined by a Delaunay triangulation of the domain 
%  of the basis functions, along with, optionally, a decomposition of the 
%  domain into subdomains.  
%  See pde toolbox function decsg for a description of the decomposed 
%  geometry concept.
%
%  The finite elements used for functional data analysis are second
%  order Lagrangian elements.  These are triangles covering a region,
%  and the basis system is piecewise quadratic.  There is a basis
%  function associated with each node in the system.
%  Nodes are points that are either:
%    vertices of triangles or
%    midpoints of edges of triangles.
%  There are 6 nodes for each triangle, and there are six degrees of 
%  freedom in a quadratic function in X and Y.  Consequently, the 6 values
%  of the nodal basis functions uniquely define a quadratic function over
%  a triangle.
%
%  Arguments:
%  P  ...  The NBASIS by 2 matrix of vertices of triangles containing
%          the X- and Y-coordinates of the vertices.
%          P may also be provided as a 2 by NBASIS matrix.
%  E  ...  The number of edges by 7 matrix defining the segments of the 
%          boundary of the region which are also edges of the triangles
%          adjacent to the boundary.  
%          The values in matrix E as follows:
%          Columns 1 and 2:  the indices of the vertices in matrix P of
%          the starting and ending points of the edges
%          Columns 3 and 4:  the starting and ending parameter values.  
%          These are initialized to 0 and 1.
%          Column  5:        the edge segment number 
%          Columns 6 and 7:  the left and right-hand side subdomain numbers
%          E may also be provided as a 7 by number of edges matrix.
%  T  ...  The 4 by no. of triangles matrix specifying triangles and 
%          their properties:
%          Columns 1 to 3:   the indices in P of the vertices of each 
%          triangle in counter-clockwise order.
%          Column  4:        the subdomain number of the triangle
%          T may also be provided with 3 columns, in which case the
%          final column is set to ones.  It may also be provided as
%          either a 3 or 4 by number of triangles matrix.
%  ORDER . Order of elements, which may be either 1 or 2 (default)
%  DL ...  Decomposed geometry object as produced by pdetools function
%          decsg.  This defines the characteristics of the subdomains
%          and the boundary segments.  It is optional.
%  
%  Returns:
%  An object of the basis class with parameters stored in member params,
%  which in this case is a struct object with members p, e and t.

%  Last modified 20 August 2010 by Jim Ramsay

%  check for number of arguments

if nargin < 3
    error('Less than three input arguments.');
end

%  default parameter values

if nargin < 5,  dl = [];    end
if nargin < 4,  order = 2;  end

%  check t for having 3 rows or columns, and add 1's to
%  make the dimension 4

if     size(t,2) == 3
    %  add a column of 1's
    t = [t, ones(size(t,1),1)];
end

if size(t,1) == 3
    %  add a row of 1's
    t = [t; ones(size(t,1),1)'];
end

%  check dimensions of P, E and T and transpose if necessary

if isempty(e)
    if size(p,1) ~= 2 && size(p,2) ~= 2 || ...
       size(t,1) ~= 4 && size(t,2) ~= 4
        error('Dimensions of at least one of P, E and T are not correct.');
    end    
else
    if size(p,1) ~= 2 && size(p,2) ~= 2 || ...
       size(e,1) ~= 7 && size(e,2) ~= 7 || ...
       size(t,1) ~= 4 && size(t,2) ~= 4
        error('Dimensions of at least one of P, E and T are not correct.');
    end
end

if size(p,2) ~= 2 && size(p,1) == 2
    p = p';
end
if size(e,2) ~= 7 && size(e,1) == 7
    e = e';
end
if size(t,2) ~= 4 && size(t,1) == 4
    t = t';
end

type     = 'FEM';
%  Argument rangeval is not needed for an FEM basis since domain
%  boundary is defined by the triangular mesh.
rangeval = [];  
%  The number of basis functions corresponds to the number of vertices
%  for order = 1, and to vertices plus edge midpoints for order = 2

%  set up the nodal points and corresponding mesh:  
%    this is p' and t(1:3,:)' in the order 1 case, but
%    includes mid-points in the order 2 case

nodeStruct = makenodes_new(p, t(:,1:3), order);

%  The params argument is a struct object 

petstr.p         = p;
petstr.e         = e;
petstr.t         = t;
petstr.order     = order;
petstr.nodes     = nodeStruct.nodes;
petstr.nodeindex = nodeStruct.nodeindex;
petstr.J         = nodeStruct.J;
petstr.metric    = nodeStruct.metric;
petstr.dl        = dl;

params = petstr;

nbasis = size(nodeStruct.nodes,1);

basisobj = basis(type, rangeval, nbasis, params);