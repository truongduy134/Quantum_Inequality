function varargout = nnz(varargin)
%NNZ (overloaded)
%
%    n = nnz(X)
%
% The NNZ operator is implemented using the concept of nonlinear operators
% in YALMIP. NNZ(X) creates a new so called derived variable that can be
% treated as any other variable in YALMIP. When SOLVESDP is issued,
% logic constraints are added to the problem to model the NNZ operator.

% Author Johan L�fberg
% $Id: nnz.m,v 1.15 2007-08-02 19:17:36 joloef Exp $

switch class(varargin{1})

    case 'sdpvar'
        % Simple binary vector?
        if all(ismember(getvariables(varargin{1}),yalmip('binvariables')))
            if is(varargin{1},'lpcone')
                varargout{1} = sum(reshape(varargin{1},prod(size(varargin{1})),1));
                return
            end
        end
        % Nope, more advanced variable.
        varargout{1} = yalmip('define','nnz',varargin{:});

    case 'char'
        z = varargin{2};
        x = varargin{3};
        switch varargin{1}
            case 'graph'
                F = nnz_internal(z,x,0);
                properties = struct('convexity','convex','monotonicity','none','definiteness','none','model','graph');
            case {'exact','integer'}
                F = nnz_internal(z,x,1);
                properties = struct('convexity','none','monotonicity','none','definiteness','none','model','integer');
            otherwise
                error('Unexpected call to NNZ OPERATOR')
        end
        varargout{1} = F;
        varargout{2} = properties ;
        varargout{3} = x;
    otherwise
end
