function result = testFunc(x, varargin)
	TotalNumOfVar = nargin

	if isempty(varargin)
		result = x;
	else
		result = x + square(varargin{1});
	end
end

function result = square(x)
	result = x * x;
end