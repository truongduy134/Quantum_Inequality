% The function: hasRankLoop(M, L, order, varProperties)
%
% This function is only applicable to projectors. Also, there should be only 2 parties.
%
% Input:
%	+ M: the moment matrix
% 	+ L: the list of monomials that the moment matrix M is indexed by (the order in L matters)
%	+ order: the order (or level) of the semidefinite programming.
%	+ varProperties: a cell of size 1 x 2. 
%			Each element in m_varProperties is a cell of size 1 x m where m is the number of inputs of
%			the measurements.
%			Each cell has elements as arrays. Each array contains integers representing variables.
%
% Output:
%	+ rankLoop = 1 if there is a rank loop. Otherwise, rankLoop = 0.
%	+ rankMoment: the rank of the moment matrix in case there is a rank loop. Otherwise, the value of rankMoment is ignored!
%	+ epsilon: a cut-off threshold so that rankLoop can occur.
function [rankLoop rankMoment epsilon] = hasRankLoop(momentMatrix, order, listMono, varProperties)
	if order < 2
		% Rank loop is not applicable to the semidefinite program at level 1
		rankLoop = 'Rank loop is not applicable to the semidefinite program at level 1';
		rankMoment = 0;
		return;
	end

	% Pre-processing:
	% Enumerate two lists containing variables in each partition.
	varParOne = cell2mat(varProperties{1});
	numParOne = length(varParOne);
	varParTwo = cell2mat(varProperties{2});
	numParTwo = length(varParTwo);

	epsilon = 1e-12;

	while epsilon < 1e-5
		rankMoment = computeRank(momentMatrix, epsilon);

		rankLoop = 1;		% Initialization
		for indexOne = 1 : numParOne
			for indexTwo = 1 : numParTwo
				indexSubMat = chooseSubMatForRankLoop(listMono, varParOne(indexOne), varParTwo(indexTwo), order);

				% Since the moment matrix M is Hermitian, submatrices of the form M(indexArr, indexArr) are Hermitian too!
				rankSubMat = computeRank(momentMatrix(indexSubMat, indexSubMat), epsilon);
		
				if rankSubMat ~= rankMoment
					% Rank loop does not occur
					rankLoop = 0;
					break;
				end
			end

			if ~rankLoop
				break;
			end
		end

		if rankLoop
			return;
		end

		epsilon = epsilon * 10;
	end

	% If the function reaches this point, rank loop does not occur
	rankLoop = 0;
end
