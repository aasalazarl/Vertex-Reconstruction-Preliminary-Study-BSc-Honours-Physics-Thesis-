function referenceNcenteredPoints
% Input data from the Python code: 
% threePlaneSimulation_referencePointsGenerator.py.
% Focused on Gaussian smearing along the x-axis.
%
% Known vertex and points on different planes:
%   vertex: (1.5, 1.5, -2.5)
%   planeA:
%       trace1: (2.3899, 1.5000, -5.0000)
%       trace2: (0.5701, 1.5000, -5.0000)
%   planeB:
%       trace1: (2.4991, 1.5000, -5.3000)
%       trace2: (0.4609, 1.5000, -5.3000)
%   planeC:
%       trace1: (2.6083, 1.5000, -5.6000)
%       trace2: (0.3517, 1.5000, -5.6000)

%-----------------------------------------------------------
%-----------------------------------------------------------
% Distances between reference and generated vertices.
vertices_distances = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Master loop.

for j = linspace(1, 1, 1) % linspace(1, N, N) looks at N decays.

%-----------------------------------------------------------
%-----------------------------------------------------------
% Tracing.

%-----------------------------------------------------------
% reference trace.

% Plane A.
trace1A = [2.3899 1.5000 -5.0000];
trace2A = [0.5701 1.5000 -5.0000];
planeA = [trace1A; trace2A];

% Plane B.
trace1B = [2.4991 1.5000 -5.3000];
trace2B = [0.4609 1.5000 -5.3000];
planeB = [trace1B; trace2B];

% Plane C.
trace1C = [2.6083 1.5000 -5.6000];
trace2C = [0.3517 1.5000 -5.6000];
planeC = [trace1C; trace2C];

planes = [planeA; planeB; planeC];

%-----------------------------------------------------------
% generated trace.

% Generate mu values.
mu = [];
for i = linspace(1, 6, 6)
    mu_dummy = planes(i);
    mu = [mu mu_dummy];
end

% Generate standard deviation, sigma, for Gaussian smearing.
sigma = 0.09; % Make sigma largen than 0.005 for illustrative purposes.

% Generate generated points.
planes_centered = planes;

for i = linspace(1, 6, 6)   % Gaussian smearing.
    planes_centered(i) = GaussianSmearing(mu(i), sigma);
end

%-----------------------------------------------------------
%-----------------------------------------------------------
% Least-squares fit.

%-----------------------------------------------------------
% Organize reference data.

z_reference1 = [];
z_reference2 = [];

x_reference1 = [];
x_reference2 = [];

for i = linspace(1, 6, 6)
    
    % Separate traces. Trace 1 in row 1 and trace 2 in row 2.
    
    if rem(i, 2) ~= 0    % Trace 1.
        % Actualize trace 1 in z_reference.
        z_dummy = planes(i, 3);
        z_reference1 = [z_reference1 z_dummy];
    
        % Actualize trace 1 in x_reference.
        x_dummy = planes(i);
        x_reference1 = [x_reference1 x_dummy];
        
    else     % Trace 2.
        % Actualize trace 2 in z_reference.
        z_dummy = planes(i, 3);
        z_reference2 = [z_reference2 z_dummy];
        
        % Actualize trace 2 in x_reference.
        x_dummy = planes(i);
        x_reference2 = [x_reference2 x_dummy];
    end
end

% Summarize generated coordinates.

z_reference = [z_reference1; z_reference2];
x_reference = [x_reference1; x_reference2];

% Least-squares fit for generated data.
c1 = linefit(z_reference(1, :), x_reference(1, :));
c2 = linefit(z_reference(2, :), x_reference(2, :));

zFitR = linspace(0, -7.0000, 100);

xFitR1 = c1(1) .* zFitR + c1(2);
xFitR2 = c2(1) .* zFitR + c2(2);

% Summarize linear fit.
xFitR = [xFitR1; xFitR2];


%-----------------------------------------------------------
% Organize generated data.

z_centered1 = [];
z_centered2 = [];

x_centered1 = [];
x_centered2 = [];

for i = linspace(1, 6, 6)
    
    % Separate traces. Trace 1 in row 1 and trace 2 in row 2.
    
    if rem(i, 2) ~= 0    % Trace 1.
        % Actualize trace 1 in z_centered.
        z_dummy = planes_centered(i, 3);
        z_centered1 = [z_centered1 z_dummy];
    
        % Actualize trace 1 in x_centered.
        x_dummy = planes_centered(i);
        x_centered1 = [x_centered1 x_dummy];
        
    else     % Trace 2.
        % Actualize trace 2 in z_centered.
        z_dummy = planes_centered(i, 3);
        z_centered2 = [z_centered2 z_dummy];
        
        % Actualize trace 2 in x_centered.
        x_dummy = planes_centered(i);
        x_centered2 = [x_centered2 x_dummy];
    end
end

% Summarize generated coordinates.

z_centered = [z_centered1; z_centered2];
x_centered = [x_centered1; x_centered2];

% Least-squares fit for generated data.
c1 = linefit(z_centered(1, :), x_centered(1, :));
c2 = linefit(z_centered(2, :), x_centered(2, :));

zFitC = linspace(0, -7.000, 100);

xFitC1 = c1(1) .* zFitC + c1(2);
xFitC2 = c2(1) .* zFitC + c2(2);

% Summarize linear fit.
xFitC = [xFitC1; xFitC2];

%-----------------------------------------------------------
%-----------------------------------------------------------
% Vertices.

%-----------------------------------------------------------
% Reference vertex.
vR = [1.5 1.5 -2.5];

%-----------------------------------------------------------
% generated vertex.
[zG, xG] = polyxpoly(zFitC, xFitC(1, :), zFitC, xFitC(2, :));
vG = [xG, vR(2), zG];

%-----------------------------------------------------------
% Distances between vertices.
vertices_dummy = sqrt(sum(((vG - vR).^2)));
vertices_distances = [vertices_distances vertices_dummy];

%----------------------------------------------------------
%----------------------------------------------------------
% Generate plot.

if j == 1
    
    figure(1)

    % Reconstruction from Guassian smearing.
    % Centered points.
    p1 = plot(z_centered(1, :), x_centered(1, :), '.r', ...
        z_centered(2, :), x_centered(2, :), '.r');
    hold on
    % Traces for centered points.
    p2 = plot(zFitC, xFitC(1, :), 'red', zFitC, xFitC(2, :), 'red');
    hold on
    
    % Reference points.    
    p3 = plot(zFitC(72), xFitC(1, 72), '.b', zFitC(76), ...
        xFitC(1, 76), '.b', zFitC(80), xFitC(1, 80), '.b');
    hold on
    p4 = plot(zFitC(72), xFitC(2, 72), '.b', zFitC(76), ...
        xFitC(2, 76), '.b', zFitC(80), xFitC(2, 80), '.b');
    hold off
    
    xlabel('z-position [m]')
    ylabel('x-position [m]')
    
    legend([p1(1), p2(1), p3(1)], 'Centered points (reality) / generated points (simulation)', ...
        'Real decay (reality) / reference decay (simulation)', 'Real points (reality) / reference points (simulation')
    
    print(gcf, 'referenceNcenteredPoints', '-dpng', '-r300')
    end
end

end