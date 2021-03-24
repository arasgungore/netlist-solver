function plotPowerOverLoadResistance(filename, R_loadname, R_values)
%plotPowerOverLoadResistance(filename, R_loadname, R_values):
%   filename:   directory/name of the input file which contains the circuit
%               data in netlist format
%   R_loadname: name of the resistor in the input file we want to pick as
%               load resistor
%   R_values:   a vector containing the load resistance values we want to plot
%
%   Reads data from the input file, calculates the node voltage values using an algorithm for Modified Node Analysis,
%   and then plots the power dissipated by the load resistor over load resistance, with a marker pointing the peak of the graph


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

% gets the load resistor's index and which nodes the load resistor is connected to
is_load = find(strcmp(element_names, R_loadname));
node1 = first_nodes(is_load);
node2 = second_nodes(is_load);

% exception handling for incorrect file formatting
assert(all(is_voltage | is_current | is_resistance), 'Unidentified element names exist at ''%s''', filename)
%assert(all(first_nodes < second_nodes), 'Some first nodes are bigger than second nodes at ''%s''', filename)
assert(all(values(is_resistance) >= 0), 'Negative resistance values exist at ''%s''', filename)
assert(nnz(is_load) == 1, 'There is no/more than one variable named %s at ''%s''', R_loadname, filename)

m = nnz(is_voltage);            % gets the number of independent voltage sources and stores it in "m"
n = max(second_nodes);          % gets the number of nodes and stores it in "n"
R_len = length(R_values);
P_values = zeros(1, R_len);     % preallocating and padding "P_values" vector where we are going to store dissipated powers for corresponding R_values

for j = 1:R_len                     % for each individual resistance value in vector "R_values"
    values(is_load) = R_values(j);  % since we do our calculations with data stored in "values"
    % we change the load resistance value stored in "values(is_load)" to a new value     (we disregard its resistance value written in the input file)

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
            B(first_nodes(v_ind), voltage_no) = -1; % if no voltage source is connected to node k, then B(:,k) is 0
        end         % in our input file, negative terminal corresponds to the first node and positive terminal corresponds to the second node
        if second_nodes(v_ind) > 0
            B(second_nodes(v_ind), voltage_no) = 1;     % ground (0) is excluded again, since 0 is not a valid index
        end
        voltage_no = voltage_no + 1;
    end

    % the C matrix is mxn and is determined by the connection of the voltage sources
    C = B.';            % the C matrix equals to the transpose of the B matrix

    % the D matrix is mxm and is zero if only independent sources are considered
    D = zeros(m);       % D is a mxm zero matrix, since there is no dependent voltage sources

    A = [G, B; C, D];   % thus we create the A coefficient matrix which is (n+m)x(n+m)


    % the i matrix is nx1 and contains the sum of the currents through the passive elements into the corresponding node
    i = zeros(n,1);     % the value of each element of i is determined by the sum of current sources into the corresponding node
    for node = 1:n      % if k'th node is a second node of a current source, we add that value to i(k)
        i(node) = sum(values(second_nodes == node & is_current)) - sum(values(first_nodes == node & is_current));
    end                 % if k'th node is a first node of a current source, we subtract that value from i(k)

    % the e matrix is mx1 and holds the values of the independent voltage sources
    e = values(is_voltage);
    z = [i; e];

    % solves the matrix equation A * x = z and stores the solution in a (m+n)x1 matrix "x"
    x = A \ z;          % A * x = z  -->  x = A-1 * z   where A-1 is the inverse of A matrix   (apparently A \ z is faster than inv(A) * z)
    if node1 == 0
        V1 = 0;         % if node1 is connected to the ground, then V1 = 0
    else
        V1 = x(node1);  % if node1 is not ground, then we get the voltage of node1 from x(node1)
    end
    if node2 == 0
        V2 = 0;         % if node2 is connected to the ground, then V2 = 0
    else
        V2 = x(node2);  %if node2 is not ground, then we get the voltage of node2 from x(node2)
    end

    P_values(j) = (V2 - V1)^2 / values(is_load);    % calculates and assigns the power dissipated by the load resistor to P_values(j)
end             % Power = Voltage^2 / Resistance


% creates a new figure and plots power (y-axis) over load resistance (x-axis)
figure
plot(R_values, P_values)
title('Power Dissipation vs Load Resistance', 'Color', 'b', 'FontSize', 14)
subtitle(['$P_{Output} = \frac{(V_' int2str(node2) ' - V_' int2str(node1) ')^2}{R_{Load}}$'], 'FontSize', 14, 'Interpreter', 'latex')
xlabel('{\itR_{Load}} (\Omega): Load Resistance')   % I used LaTeX interpreter to write the power formula as subtitle
ylabel('{\itP_{Output}} (W): Power Dissipated')

% finds the maximum power dissipated by the load resistor, which corresponds to the maximum of the power-resistance curve
imax = max(P_values) == P_values;
P_peak = P_values(imax);
R_peak = R_values(imax);            % writes the corresponding P and R values of the point that is maximum of the curve to the plot
text(R_peak, P_peak, ['{\itP_{Max}} = ' num2str(P_peak, '%f') ' W'], 'VerticalAlignment', 'bottom', 'FontSize', 8);
text(R_peak, P_peak, ['{\itR_{Load}} = ' num2str(R_peak, '%f') ' \Omega'], 'VerticalAlignment', 'top', 'FontSize', 8);
hold on
plot(R_peak, P_peak, 'r*')          % puts a marker on the maximum of the power-resistance curve
%legend('Power-Resistance Curve', 'Maximum Power Dissipated')

fprintf('Max Power Dissipated by the Load Resistor %s: %f W\n', R_loadname, P_peak)
fprintf('Load Resistance of %s at Max Power Output: %f Ohm\n\n', R_loadname, R_peak)

end