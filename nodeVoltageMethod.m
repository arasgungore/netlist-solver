function nodeVoltages = nodeVoltageMethod(filename)
%nodeVoltageMethod(filename):
%   filename:   directory/name of the input file which contains the circuit
%               data in netlist format
%
%   Reads data from the input file and then calculates the node voltage values using an algorithm for Modified Node Analysis
%   For more information about the algorithm: http://www.swarthmore.edu/NatSci/echeeve1/Ref/mna/MNA3.html
%   Returns a column vector containing the voltages of nodes 1-N in order


fileID = fopen(filename);       % opens the input file
assert(fileID ~= -1, 'Could not open file ''%s''', filename)    % throws exception if unable to open the input file
input_data = textscan(fileID, '%s %d %d %f');                   % reads data from the input file and stores it in a cell array called "input_data"
closed = fclose(fileID);        % closes the input file
assert(closed ~= -1, 'Could not close file %d: ''%s''', fileID, filename)   % throws exception if unable to close the input file

% gets the columns of the "input_data", categorizes and stores it in column vectors
element_names = input_data{1};      % first column of the input file
first_nodes = input_data{2}(:);     % second column of the input file
second_nodes = input_data{3}(:);    % third column of the input file
values = input_data{4}(:);          % fourth column of the input file

% gets the logical arrays of the elements where 1's represent that an element starting with V, I or R exists on those indices
is_voltage = startsWith(element_names, 'V');
is_current = startsWith(element_names, 'I');
is_resistance = startsWith(element_names, 'R');

% exception handling for incorrect file formatting
assert(all(is_voltage | is_current | is_resistance), 'Unidentified element names exist at ''%s''', filename)
%assert(all(first_nodes < second_nodes), 'Some first nodes are bigger than second nodes at ''%s''', filename)
assert(all(values(is_resistance) >= 0), 'Negative resistance values exist at ''%s''', filename)

m = nnz(is_voltage);            % gets the number of independent voltage sources and stores it in "m"
n = max(second_nodes);          % gets the number of nodes and stores it in "n"


% the G matrix is nxn and is determined by the interconnections between the resistors
G = zeros(n);       % creates a nxn zero matrix
for node = 1:n      % each element in the diagonal of G is equal to the sum of the conductance (1/R) of each element connected to the corresponding node
    G(node, node) = sum(1 ./ values((first_nodes == node | second_nodes == node) & is_resistance));
end
for r_ind = find(is_resistance).'
    R_first_nodes = first_nodes(r_ind);
    R_second_nodes = second_nodes(r_ind);
    if R_first_nodes > 0 && R_second_nodes > 0      % ground (0) is excluded
        G(R_first_nodes, R_second_nodes) = G(R_first_nodes, R_second_nodes) - 1/values(r_ind);
        G(R_second_nodes, R_first_nodes) = G(R_second_nodes, R_first_nodes) - 1/values(r_ind);
    end     % each off diagonal element of G is equal to the sum of the negative conductance of each element connected to the pair of corresponding node
end

% the B matrix is nxm and is determined by the connection of the voltage sources
B = zeros(n,m);
voltage_no = 1;
for v_ind = find(is_voltage).'  % if the positive terminal of the i'th voltage source is connected to node k, then B(i,k) is 1
    if first_nodes(v_ind) > 0   % if the negative terminal of the i'th voltage source is connected to node k, then B(i,k) is -1
        B(first_nodes(v_ind), voltage_no) = -1; %if no voltage source is connected to node k, then B(:,k) is 0
    end         % in our input file, negative terminal corresponds to the first node and positive terminal corresponds to the second node
    if second_nodes(v_ind) > 0
        B(second_nodes(v_ind), voltage_no) = 1;     % ground (0) is excluded again, since 0 is not a valid index
    end     % normally, we would expect second node to be always greater than first node, and hence second node > 0 (ground)
    voltage_no = voltage_no + 1;    % however, I included this extra if condition (line 57) to avoid errors for circuit elements
end                                 % where second nodes may come in smaller than first nodes, just to be on the safe side

% the C matrix is mxn and is determined by the connection of the voltage sources
C = B.';        % the C matrix equals to the transpose of the B matrix

% the D matrix is mxm and is zero if only independent sources are considered
D = zeros(m);   % D is a mxm zero matrix, since there is no dependent voltage sources

A = [G, B; C, D];   % thus we create the A coefficient matrix which is (n+m)x(n+m)


% the i matrix is nx1 and contains the sum of the currents through the passive elements into the corresponding node
i = zeros(n,1);     % the value of each element of i is determined by the sum of current sources into the corresponding node
for node = 1:n      % if k'th node is a second node of a current source, we add that value to i(k)
    i(node) = sum(values(second_nodes == node & is_current)) - sum(values(first_nodes == node & is_current));
end                 % if k'th node is a first node of a current source, we subtract that value from i(k)

% the e matrix is mx1 and holds the values of the independent voltage sources
e = values(is_voltage);

z = [i; e];         % thus we create the z right-hand side matrix which is (n+m)x1

% solves the matrix equation A * x = z and stores the solution in a (m+n)x1 matrix "x"
x = A \ z;              % A * x = z  -->  x = A-1 * z   where A-1 is the inverse of A matrix   (apparently A \ z is faster than inv(A) * z)
nodeVoltages = x(1:n);  % since the x(1:n) contains the node voltage values, and the remaining x(n+1:end) contains the currents
                        % through the voltage sources, we take x(1:n) as output
end