%% Main Program
clear, clc, close all       % clear the workspace, clear the command window, and then close any open figures

help nodeVoltageMethod      % get information about the "nodeVoltageMethod" function


filename = 'input_1.txt';   % name of the input file that contains the circuit data in netlist format
node_voltages = nodeVoltageMethod(filename)


filename = input('Enter directory/name of the input file including extension (''.txt''): ', 's');
node_voltages = nodeVoltageMethod(filename)

%% Bonus Program
clear, clc, close all       % clear the workspace, clear the command window, and then close any open figures

help plotPowerOverLoadResistance    % get information about the "plotPowerOverLoadResistance" function


filename = 'input_1.txt';           % name of the input file that contains the circuit data in netlist format
R_loadname = 'R2';                  % chooses R2 as load resistor
R_values = linspace(0, 20, 200);    % assigns 200 resistance values ranging from 0 to 20 Ohm to R2 which we are going to plot
plotPowerOverLoadResistance(filename, R_loadname, R_values)


filename = input('Enter directory/name of the input file including extension (''.txt''): ', 's');
R_loadname = input('Enter the name of the resistor you want to pick as load resistor: ', 's');
R_begin = input('Enter the initial resistance value (in Ohm): ');
R_end = input('Enter the final resistance value (in Ohm): ');
R_points = input('Enter the number of points to plot: ');
plotPowerOverLoadResistance(filename, R_loadname, linspace(R_begin, R_end, R_points))