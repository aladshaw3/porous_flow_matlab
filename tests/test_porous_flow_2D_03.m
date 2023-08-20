%% Test 4 - Multi-domain properties and reactions

% Create geometry for the pde
%   2D geometry 

% Rectangle is code 3, 4 sides (in m)
%          x1,x2,x3,x4       y1,y2,y3,y4
R1 = [3,4, 0, 0.05, 0.05, 0,        0, 0, 0.025, 0.025]';
R2 = [3,4, 0.05, 0.1, 0.1, 0.05,        0, 0, 0.025, 0.025]';
geom = [R1,R2];

% Names for the two geometric objects
ns = (char('R1','R2'))';

% Set formula
sf = 'R1 + R2';

% Create geometry
% E1 - bottom, E2 - right, E3 - top, E4 - left
g = decsg(geom,sf,ns);

% Plot results
f1 = figure;
pdegplot(g,"EdgeLabels","on","FaceLabels","on")
saveas(f1,'output/Gifs/pressure_test03_subdomains.png');
close(f1);

% BCs - boundary id 1 = input, positive velocity
%       boundary id 2 = output,
%       face id 1 = subdomain 1
%       face id 2 = subdomain 2

% Rxns
%       Subdomain 1
%           A --> B
%           B --> C
%       Subdomain 2
%           A --> D
%           D --> A
%
% Mobility: B and D are immobile, A and C are mobile 