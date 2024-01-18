function vertexReconstruction_3D
% Input data from the Python code: 
% threePlaneSimulation_referencePointsGenerator.py.
% Focused on Gaussian smearing along the x-axis and y-axis.
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

for j = linspace(1, 1000, 1000) % linspace(1, N, N) looks at N decays.

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
mu_x = [];
mu_y = [];
for i = linspace(1, 6, 6)
    % Acutalize mu_x.
    mu_x_dummy = planes(i);
    mu_x = [mu_x mu_x_dummy];
    
    % Actualize mu_y.
    mu_y_dummy = planes(i, 2);
    mu_y = [mu_y mu_y_dummy];
end

% Generate standard deviation, sigma, for Gaussian smearing.
sigma = 0.005;

% Generate generated points.
planes_generated = planes;

for i = linspace(1, 6, 6)   % Gaussian smearing.
    planes_generated(i) = GaussianSmearing(mu_x(i), sigma);
    planes_generated(i, 2) = GaussianSmearing(mu_y(i), sigma);
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

y_reference1 = [];
y_reference2 = [];

for i = linspace(1, 6, 6)
    
    % Separate traces. Trace 1 in row 1 and trace 2 in row 2.
    
    if rem(i, 2) ~= 0    % Trace 1.
        % Actualize trace 1 in z_reference.
        z_dummy = planes(i, 3);
        z_reference1 = [z_reference1 z_dummy];
    
        % Actualize trace 1 in x_reference.
        x_dummy = planes(i);
        x_reference1 = [x_reference1 x_dummy];
        
        % Actualize trace 1 in y_reference.
        y_dummy = planes(i, 2);
        y_reference1 = [y_reference1 y_dummy];
        
    else     % Trace 2.
        % Actualize trace 2 in z_reference.
        z_dummy = planes(i, 3);
        z_reference2 = [z_reference2 z_dummy];
        
        % Actualize trace 2 in x_reference.
        x_dummy = planes(i);
        x_reference2 = [x_reference2 x_dummy];
        
        % Actualize trace 2 in y_reference.
        y_dummy = planes(i, 2);
        y_reference2 = [y_reference2 y_dummy];
    end
end

% Summarize generated coordinates.

z_reference = [z_reference1; z_reference2];
x_reference = [x_reference1; x_reference2];
y_reference = [y_reference1; y_reference2];

% Least-squares fit for generated data.
c1 = linefit3D(z_reference(1, :), x_reference(1, :), y_reference(1, :));
c2 = linefit3D(z_reference(2, :), x_reference(2, :), y_reference(2, :));

% Generate zFitR and xFitR.
zFitR = linspace(0, -7.0000, 100);

d1 = linefit(z_reference(1, :), x_reference(1, :));
d2 = linefit(z_reference(2, :), x_reference(2, :));

xFitR1 = d1(1) .* zFitR + d1(2);
xFitR2 = d2(1) .* zFitR + d2(2);

xFitR = [xFitR1; xFitR2];

% Generate yFitR.

yFitR1 = c1(1) .* zFitR + c1(2) .* xFitR(1, :) + c1(3);
yFitR2 = c2(1) .* zFitR + c2(2) .* xFitR(2, :) + c2(3);

% Summarize linear fit.
yFitR = [yFitR1; yFitR2];


%-----------------------------------------------------------
% Organize generated data.

z_generated1 = [];
z_generated2 = [];

x_generated1 = [];
x_generated2 = [];

y_generated1 = [];
y_generated2 = [];

for i = linspace(1, 6, 6)
    
    % Separate traces. Trace 1 in row 1 and trace 2 in row 2.
    
    if rem(i, 2) ~= 0    % Trace 1.
        % Actualize trace 1 in z_generated.
        z_dummy = planes_generated(i, 3);
        z_generated1 = [z_generated1 z_dummy];
    
        % Actualize trace 1 in x_generated.
        x_dummy = planes_generated(i);
        x_generated1 = [x_generated1 x_dummy];
        
        % Actualized trace 1 in y_generated.
        y_dummy = planes_generated(i, 2);
        y_generated1 = [y_generated1 y_dummy];
        
    else     % Trace 2.
        % Actualize trace 2 in z_generated.
        z_dummy = planes_generated(i, 3);
        z_generated2 = [z_generated2 z_dummy];
        
        % Actualize trace 2 in x_generated.
        x_dummy = planes_generated(i);
        x_generated2 = [x_generated2 x_dummy];
        
        % Actualize trace 2 in y_generated.
        y_dummy = planes_generated(i, 2);
        y_generated2 = [y_generated2 y_dummy];
    end
end

% Summarize generated coordinates.

z_generated = [z_generated1; z_generated2];
x_generated = [x_generated1; x_generated2];
y_generated = [y_generated1; y_generated2];

% Least-squares fit for generated data.
c1 = linefit3D(z_generated(1, :), x_generated(1, :), y_generated(1, :));
c2 = linefit3D(z_generated(2, :), x_generated(2, :), y_generated(2, :));

% Generate zFitG and xFitG.
zFitG = linspace(0, -7.0000, 100);

d1 = linefit(z_generated(1, :), x_generated(1, :));
d2 = linefit(z_generated(2, :), x_generated(2, :));

xFitG1 = d1(1) .* zFitG + d1(2);
xFitG2 = d2(1) .* zFitG + d2(2);

xFitG = [xFitG1; xFitG2];

% Generate yFitR.

yFitG1 = c1(1) .* zFitG + c1(2) .* xFitG(1, :) + c1(3);
yFitG2 = c2(1) .* zFitG + c2(2) .* xFitG(2, :) + c2(3);

% Summarize linear fit.
yFitG = [yFitG1; yFitG2];

%-----------------------------------------------------------
%-----------------------------------------------------------
% Vertices.

%-----------------------------------------------------------
% reference vertex.
vE = [1.5 1.5 -2.5];

%-----------------------------------------------------------
% generated vertex.
[zO, xO] = polyxpoly(zFitG, xFitG(1, :), zFitG, xFitG(2, :));
yO = ((c1(1) * zO + c1(2) * xO + c1(3)) + ...
    (c2(1) * zO + c2(2) * xO + c2(3))) / 2;
vO = [xO, yO, zO];

%-----------------------------------------------------------
% Distances between vertices.
vertices_dummy = sqrt(sum(((vO - vE).^2)));
vertices_distances = [vertices_distances vertices_dummy];

%----------------------------------------------------------
%----------------------------------------------------------
% Generate plot.

if j == 1
    
    figure(1)

    % Expectation.
    % Plot reference detection coordinates.
    p1 = plot3(z_reference(1, :), x_reference(1, :), ...
        y_reference(1, :), '.b', z_reference(2, :), ...
        x_reference(2, :), y_reference(2, :), '.b');
    hold on
    % Plot reference traces.
    p2 = plot3(zFitR, xFitR(1, :), yFitR(1, :), ...
        'blue', zFitR, xFitR(2, :), yFitR(2, :), 'blue');
    hold on

    % Reconstruction from Guassian smearing.
    % Plot generated detection coordinates.
    p3 = plot3(z_generated(1, :), x_generated(1, :), ...
        y_generated(1, :), '.r', z_generated(2, :), ...
        x_generated(2, :), y_generated(2, :), '.r');
    hold on
    % Plot generated traces.
    p4 = plot3(zFitG, xFitG(1, :), yFitG(1, :), 'red', ...
        zFitG, xFitG(2, :), yFitG(2, :), 'red');
    legend([p1(1), p2(1), p3(1), p4(1)], ...
        'Reference detection points' ,'Reference decay', ...
        'Detection points from Gaussian smearing', ...
        'generated decay from Guassian smearing')
    hold off
    
    xlabel('z-position [m]')
    ylabel('x-position [m]')
    zlabel('y-position [m]')
    zlim([1, 2])
    
    print(gcf, '3DvertexGeneration', '-dpng', '-r300')
    
    end

end

%-----------------------------------------------------------
%-----------------------------------------------------------
% Generate a histogram.
figure(2)
hist(vertices_distances, 10)
H = hist(vertices_distances, 10)

xlabel('Vertex-to-vertex distance [m]')
ylabel('Count')

%print(gcf, 'vertex-to-vertex_distances3D', '-dpng', '-r300')
print(gcf, '3Dvertex2vertexDistances', '-dpng', '-r300')

poissonParameter = poissfit(vertices_distances)

bins = 0:9;
obsCounts = H;
n = sum(obsCounts);
pd = fitdist(bins', 'Poisson', 'Frequency', obsCounts');
expCounts = n *pdf(pd, bins);

[h,p,st] = chi2gof(bins,'Ctrs',bins,...
                        'Frequency',obsCounts, ...
                        'Expected',expCounts,...
                        'NParams',1)

end