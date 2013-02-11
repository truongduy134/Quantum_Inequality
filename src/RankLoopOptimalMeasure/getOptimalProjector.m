% The function: getOptimalProjector(M, order, cellMono, varProperties)
%
% Important note: This function is only applicable to two-party cases only! Variables must be projectors
%
% Inputs:
%	+ M: the moment matrix.
%	+ order: the order (or level) of the semidefinite programming.
% 	+ cellMono: a cell containing monomials of degree <= the order of the semidefinite program. Note that the moment matrix is indexed by cellMono. 
%	+ varProperties: a cell of size 1 x 2. 
%			Each element in m_varProperties is a cell of size 1 x m where m is the number of inputs of
%			the measurements.
%			Each cell has elements as arrays. Each array contains integers representing variables.
%
% Outputs:
%	+ cellProjector: a cell of size N(where N is the total number of variables) where
%			cellStateAndObservable{i} (i >= 1): the matrix representation of the projector measurement whose number is i. 
%	+ opState: the optimal quantum state
function [cellProjector opState] = getOptimalProjector(momentMatrix, order, cellMono, varProperties)
	cholDecompose = choleskyDecompose(momentMatrix, 1e-10);

	% Pre-processing:
	% Enumerate two lists containing variables in each partition.
	varParOne = cell2mat(varProperties{1});
	numParOne = length(varParOne);
	varParTwo = cell2mat(varProperties{2});
	numParTwo = length(varParTwo);
	numVar = numParOne + numParTwo;

	cellProjector = cell(1, numVar);

	epsilon = 1e-10;
	% Extract optimal projector measurements
	for index = 1 : numParOne
		indexSubMat = filterMonoByLeadTerm(cellMono, varParOne(index), 1);
		cellProjector{varParOne(index)} = findProjection(cholDecompose(:, indexSubMat), epsilon);
	end

	for index = 1 : numParTwo
		indexSubMat = filterMonoByLeadTerm(cellMono, varParTwo(index), 2);
		cellProjector{varParTwo(index)} = findProjection(cholDecompose(:, indexSubMat), epsilon);
	end

	% Extract optimal state
	opState = cholDecompose(:, 1);
end