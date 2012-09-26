% The function: testModTwoGame(sdpLevel)
% Input parameter:
%	+ sdpLevel: the level of the semidefinite program which we want the solver to run at.
% If there is no input, we assume sdpLevel = 1
%
% The function runs the quantum bound problem solver to approximate the bound of the inequality of the mod-2 game and prints the results which comprise:
%		1) The bound value.
%		2) Extra information provided by underlying solvers (such as SeDuMi, YALMIP, etc.)
%
% Problem description:
%	+ Find an upper bound of:
%		P = 1/4 * (<A^1_0 * B^1_0> + <A^0_0 * B^0_0> + <A^0_0 * B^0_1> + <A^1_0 * B^1_1> +
%				   <A^0_1 * B^0_0> + <A^0_1 * B^1_1> + <A^1_1 * B^0_1> + <A^1_1 * B^1_0>)    
%	+ All variables are observables.
%	+ There are 2 partitions (parties) A = {A_1, A_2} and B = {B_1, B_2}
%	
% We map A^0_0 to 1, 
%		 A^1_0 to 2,
%		 A^0_1 to 3,
%		 A^1_1 to 4,
%		 B^0_0 to 5, 
%		 B^1_0 to 6,
%		 B^0_1 to 7,
%		 B^1_1 to 8,
%
% We rewrite P as P = 1/4 * <2, 6> + 1/4 * <1, 5> + 1/4 * <1, 7> + 1/4 * <2, 8> + 1/4 * <3, 5> + 1/4 * <3, 8> + 1/4 * <4, 7> + 1/4 * <4, 6>
function opResult = testModTwoGame(varargin)
	if nargin == 0
		sdpLevel = 1;
		hashGenMonoInfo = {};
	else
		sdpLevel = varargin{1};
		hashGenMonoInfo = varargin{2};
	end
	
	% Declare monomial
	monoOne = NonCommuteMonomial({[2], [6]}, 1/4);			% Monomial 1/4 * <2, 6> = 1/4 * <A^1_0 * B^1_0>
	monoTwo = NonCommuteMonomial({[1], [5]}, 1/4);			% Monomial 1/4 * <1, 5> = 1/4 * <A^0_0 * B^0_0>
	monoThree = NonCommuteMonomial({[1], [7]}, 1/4);		% Monomial 1/4 * <1, 7> = 1/4 * <A^0_0 * B^0_1>
	monoFour = NonCommuteMonomial({[2], [8]}, 1/4);			% Monomial 1/4 * <2, 8> = 1/4 * <A^1_0 * B^1_1>
	monoFive = NonCommuteMonomial({[3], [5]}, 1/4);			% Monomial 1/4 * <3, 5> = 1/4 * <A^0_1 * B^0_0>
	monoSix = NonCommuteMonomial({[3], [8]}, 1/4);			% Monomial 1/4 * <3, 8> = 1/4 * <A^0_1 * B^1_1>
	monoSeven = NonCommuteMonomial({[4], [7]}, 1/4);		% Monomial 1/4 * <4, 7> = 1/4 * <A^1_1 * B^0_1>	
	monoEight = NonCommuteMonomial({[4], [6]}, 1/4);		% Monomial 1/4 * <4, 6> = 1/4 * <A^1_1 * B^1_0>
	
	% Declare the polynomial whose monomials are monoOne, ..., monoEight
	% The information about the variables is represented as {{[1 2] [3 4]}, {[5 6] [7 8]}} which means:
	%		1) There are 3 partitions (Partition 1 contains variables who labels are 1, 2, 3, 4. 
	%								   Partition 2 contains variables who labels are 5, 6, 7, 8).
	%		   Note that a * b = b * a if a, b are in different partitions.
	%		2) There are 4 input groups (Input group 1 contains variables who labels are 1, 2.
	%									 Input group 2 contains variables who labels are 3, 4.
	%									 Input group 3 contains variables who labels are 5, 6.
	%									 Input group 4 contains variables who labels are 7, 8).
	%		   The sum of all variables in an input group equals the identity.
	%
	% The variable type is PROJECTOR which is mapped to 0.
	polyOp = NonCommutePolynomial({monoOne, monoTwo, monoThree, monoFour, monoFive, monoSix, monoSeven, monoEight}, {{[1 2] [3 4]}, {[5 6] [7 8]}}, 0);
	
	% Call the solver
	% We check the identity constraint
	checkIdentityConstraint = 1;		% 1 means CHECK
	result = findQuantumBound(polyOp, sdpLevel, checkIdentityConstraint, hashGenMonoInfo);
	
	% Print result
	disp('Upper bound value = ');
	disp(result{1});
	disp('Solver message = ');
	disp(result{2});
	
	disp('Rank Loop Result = ')
	epsilon = 1e-10;
	rankLoop = 0;	% Dummy
	result{3} = standardizeMatrix(result{3}, 1e-10);
	while(epsilon < 1e-5)
		rankLoop = hasRankLoop(result{3}, sdpLevel, result{5}, polyOp.m_varProperties, epsilon);
		if rankLoop
			break;
		else
			epsilon = epsilon * 10;
		end
	end
	
	if rankLoop
		disp('Rank loop occurs');
		disp('The threshold epsilon = ');
		disp(epsilon);

		opResult = getOptimalProjector(result{3}, sdpLevel, result{5}, polyOp.m_varProperties);
	else
		disp('Rank loop does not occur');
		disp(epsilon);
	end

	cholDecompose = choleskyDecompose(result{3}, 1e-10);
	cholDecompose = standardizeMatrix(cholDecompose, 1e-6);
	%cholDecompose = standardizeMatrix(cholDecompose, 1e-3);
	%diagonal = diag(cholDecompose)

	%rankMoment = rankHermitian(result{3}, 1e-8)
	%rankChol = rankHermitian(cholDecompose, 1e-8)
	%deter = det(cholDecompose)
	numMono = length(result{5});

	firstMono = NonCommuteMonomial({[3], []}, 1);
	for index = 1 : numMono
		product = mulMonomial(firstMono, result{5}{index}, polyOp.m_varType, polyOp.m_inputGroup);

		monoVarOrder = result{5}{index}.m_varOrdering
		productVarOrder = product.m_varOrdering

		if product.m_degree >= 0 && product.m_degree <= sdpLevel
			orderStr = getMonoVarOrderStr(result{5}{index});
			coord = result{4}(orderStr);

			if coord(1) ~= 1
				error('Something wrong with findQuantumBound')
			end

			column = coord(2);

			orderStrP = getMonoVarOrderStr(product);
			coordP = result{4}(orderStr);
			if coordP(1) ~= 1
				error('Something wrong with findQuantumBound')
			end

			columnP = coordP(2);

			LHS = opResult{3} * cholDecompose(:, column);
			RHS = cholDecompose(:, columnP);

			differ = LHS - RHS;
			disp('Difference = ')
			differ = standardizeMatrix(differ, 1e-10)
		elseif product.m_degree < 0
			orderStr = getMonoVarOrderStr(result{5}{index});
			coord = result{4}(orderStr);

			if coord(1) ~= 1
				error('Something wrong with findQuantumBound')
			end

			column = coord(2);

			LHS = opResult{3} * cholDecompose(:, column);

			disp('Expect zero = ')
			LHS = standardizeMatrix(LHS, 1e-10)
		end
	end

	numVarPlusOne = length(opResult);
	disp(evaluatePolynomial(polyOp, opResult(1 : (numVarPlusOne - 1)), opResult{numVarPlusOne}));
	disp('Check Moment Var Value: ');
	disp(checkMomentVarValue(result{3}, polyOp, result{5}, opResult(1 : (numVarPlusOne - 1)), opResult{numVarPlusOne}, 1e-10));
end
	