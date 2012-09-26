% Function: compareMono  % This method compares two monomials according to the degree - lexicography criteria.% It returns 1 if monoOne '>' monoTwo, -1 if monoOne '<' monoTwo, and 0 if monoOne '==' monoTwo% Assumption: The two monomials have the same set of variables and partitions.function compareResult = compareTo(monoOne, monoTwo)	if monoONe.m_degree < monoTwo.m_degree		compareResult = -1;		return;	end				if monoOne.m_degree > monoTwo.m_degree		compareResult = 1;		return;	end				% Now we consider the lexicographical criterion.	% Two nonzero monomials must have the same number of partitions in m_Ordering.	for i = 1 : length(monoOne.m_Ordering)		emptyMonoOne = isempty(monoOne.m_varOrdering{i};		emptyMonoTwo = isempty(monoTwo.m_varOrdering{i};				if emptyMonoOne && ~emptyMonoTwo			compareResult = 1;			return;		end						if ~emptyMonoOne && emptyMonoTwo			compareResult = -1;			return;		end						if ~emptyMonoOne && ~emptyMonoTwo			% Compare each variable number in each partition.			len = min(length(monoOne.m_Ordering{i}), length(monoTwo.m_Ordering{i}));			for j = 1 : len				if monoOne.m_Ordering{i}(j) < monoTwo.m_Ordering{i}(j)					compareResult = -1;					return;				end										if monoOne.m_Ordering{i}(j) > monoTwo.m_Ordering{i}(j)					compareResult = 1;					return;				end										if length(monoOne.m_Ordering{i}) < length(monoTwo.m_Ordering{i})					compareResult = -1;					return;				end										if length(monoOne.m_Ordering{i}) > length(monoTwo.m_Ordering{i})					compareResult = 1;					return;				end			end		end	end		% If the program reach this point, two monomials are the same in terms of variable ordering	compareResult = 0;end