function vertexReconstruction_2D
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

for j = linspace(1, 10000, 10000) % linspace(1, N, N) looks at N decays.

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
sigma = 0.005;

% Generate generated points.
planes_generated = planes;

for i = linspace(1, 6, 6)   % Gaussian smearing.
    planes_generated(i) = GaussianSmearing(mu(i), sigma);
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

z_generated1 = [];
z_generated2 = [];

x_generated1 = [];
x_generated2 = [];

for i = linspace(1, 6, 6)
    
    % Separate traces. Trace 1 in row 1 and trace 2 in row 2.
    
    if rem(i, 2) ~= 0    % Trace 1.
        % Actualize trace 1 in z_generated.
        z_dummy = planes_generated(i, 3);
        z_generated1 = [z_generated1 z_dummy];
    
        % Actualize trace 1 in x_generated.
        x_dummy = planes_generated(i);
        x_generated1 = [x_generated1 x_dummy];
        
    else     % Trace 2.
        % Actualize trace 2 in z_generated.
        z_dummy = planes_generated(i, 3);
        z_generated2 = [z_generated2 z_dummy];
        
        % Actualize trace 2 in x_generated.
        x_dummy = planes_generated(i);
        x_generated2 = [x_generated2 x_dummy];
    end
end

% Summarize generated coordinates.

z_generated = [z_generated1; z_generated2];
x_generated = [x_generated1; x_generated2];

% Least-squares fit for generated data.
c1 = linefit(z_generated(1, :), x_generated(1, :));
c2 = linefit(z_generated(2, :), x_generated(2, :));

zFitG = linspace(0, -7.000, 100);

xFitG1 = c1(1) .* zFitG + c1(2);
xFitG2 = c2(1) .* zFitG + c2(2);

% Summarize linear fit.
xFitG = [xFitG1; xFitG2];

%-----------------------------------------------------------
%-----------------------------------------------------------
% Vertices.

%-----------------------------------------------------------
% Reference vertex.
vR = [1.5 1.5 -2.5];

%-----------------------------------------------------------
% generated vertex.
[zG, xG] = polyxpoly(zFitG, xFitG(1, :), zFitG, xFitG(2, :));
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

    % Expectation.
    % Plot reference detection coordinates.
    p1 = plot(z_reference(1, :), x_reference(1, :), '.b', ...
        z_reference(2, :), x_reference(2, :), '.b');
    hold on
    % Plot reference traces.
    p2 = plot(zFitR, xFitR(1, :), 'blue', zFitR, xFitR(2, :), 'blue');
    hold on

    % Reconstruction from Guassian smearing.
    % Plot generated detection coordinates.
    p3 = plot(z_generated(1, :), x_generated(1, :), '.r', ...
        z_generated(2, :), x_generated(2, :), '.r');
    hold on
    % Plot generated traces.
    p4 = plot(zFitG, xFitG(1, :), 'red', zFitG, xFitG(2, :), 'red');
    legend([p1(1), p2(1), p3(1), p4(1)], 'Reference points' , ...
        'Reference decay', 'Generated points from Gaussian smearing', ...
        'Generated decay from Guassian smearing')
    hold off
    
    xlabel('z-position [m]')
    ylabel('x-position [m]')
    
    print(gcf, '2DvertexGeneration', '-dpng', '-r300')
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

%print(gcf, 'vertex-to-vertex_distances', '-dpng', '-r300')
print(gcf, '2Dvertex2vertexDistances', '-dpng', '-r300')

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