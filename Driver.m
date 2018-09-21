
[dsm,R] = geotiffread('322.tif');

% Compute MRF classification map.
tic;
radius = 20; % set to roughly the size of the largest building in the scene
del_0 = 1.5;
coherence_const = 0.75;
max_iterations = 10;
is_terrain = MRF_extract_bare_earth(dsm,radius,del_0,coherence_const,max_iterations);
t = toc

geotiffwrite('is_terrain.tif', is_terrain, R, 'CoordRefSysCode', 32617)

