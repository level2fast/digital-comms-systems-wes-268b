% Define LDPC code parameters
numVariableNodes = 8;
numCheckNodes = 4;

% Define the edges in the Tanner graph
edges = [
    1 5;
    1 6;
    2 5;
    2 6;
    3 6;
    3 7;
    4 7;
    4 8;
    5 1;
    5 2;
    5 3;
    5 4;
    6 1;
    6 2;
    6 3;
    6 4;
    7 3;
    7 4;
    8 4;
    8 1;
];

% Create the Tanner graph using the graph function
G = graph(edges(:, 1), edges(:, 2));

% Plot the Tanner graph
figure;
plot(G, 'Layout', 'force', 'NodeLabel', {});

% Customize the plot for better visualization (optional)
title('Tanner Graph for LDPC Code');
xlabel('Variable Nodes (Bits)');
ylabel('Check Nodes (Parity-Check Equations)');
