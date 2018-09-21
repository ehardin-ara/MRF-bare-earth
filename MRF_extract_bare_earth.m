function [isTerrain] = MRF_extract_bare_earth(DSM,smoothing_radius,del_0,coherence_const,max_iterations)
%MRF_extract_bare_earth - This is an earth classification algorithm.
% The MRF model is based on Baillard, C., & Ma?tre, H. (1999). 3-D
% reconstruction of urban scenes from aerial stereo imagery: a focusing
% strategy. Computer Vision and Image Understanding, 76(3), 244-258. 
% Energy is minimized with MATLAB's built-in maxflow algorithm.
%
% Syntax:  [isTerrain] = MRF_extract_bare_earth(DSM,data_window,del_0,iterations)
%
% Inputs:
%    DSM - Digital surface model as a double precision matlab array
%    smoothing_radius - scalar radius in pixels for averaging smoothing kernel. 
%    Set to roughly the size of the largest buildings. 
%    del_0 - scalar height at which classification becomes biased
%    toward off-terrain classification. Same interpretation as in Baillard
%    (1999). Values less than one building story work well (~1.5m).
%    coherence_const - constant to adjust strength of unary and
%    interaction energies. Values of about 0.75 are good starting points. Values
%    more than 0.5, e.g., 0.75, should yield a little more detail. 
%    iterations - max number of iterations of the algorithm. The algorithm
%    usually converges in about 5 iterations. The algorithm will break
%    out if convergence is reached before the max number of iterations is
%    reached.
%
% Outputs:
%    isTerrain - logical matrix with 1 indicating terrain and 0 indicating
%    off-terrain, e.g., trees or buildings. 
%
% Author: Eric Hardin
% email: ehardin@ara.com
% Sept 21, 2018

    [rows,cols] = size(DSM);
    %steps = [[0,1];[-1,0];[1,0];[0,-1]]; % first order
    steps = [[-1,1];[0,1];[1,1];[-1,0];[1,0];[-1,-1];[0,-1];[1,-1]]; % second order
    clique_indexes = BuildCliqueIndexes(rows,cols,steps);
    
    del_DSM = DSM(clique_indexes(:,1)) - DSM(clique_indexes(:,2)); % elevation differences between cliques

    c = -log(0.5); % constant so that d_0 is the Gaussian half-distance in meters
    h = fspecial('disk', smoothing_radius);


    DTM = imfilter(DSM,h);
    nDSM = DSM - DTM;
    isTerrain = nDSM < 0;

    convergence_prev = Inf;
    for iteration=1:max_iterations
        
        DTM = imfilter(double(isTerrain).*DSM,h) ./ (sum(h(:))*imfilter(double(isTerrain),h));
        nDSM = DSM - DTM;

        E = coherence_const*(1-exp(-c*((nDSM-del_0)/del_0).^2));
        E_data_terrain = zeros(size(nDSM));
        E_data_terrain(nDSM > del_0) = E(nDSM > del_0);
        E_data_building = zeros(size(nDSM));
        E_data_building(nDSM < del_0) = E(nDSM < del_0);

        E_clique = zeros(size(clique_indexes,1),1);
        terrain_terrain = isTerrain(clique_indexes(:,1)) & isTerrain(clique_indexes(:,2));
        terrain_building = isTerrain(clique_indexes(:,1)) & ~isTerrain(clique_indexes(:,2));
        building_terrain = ~isTerrain(clique_indexes(:,1)) & isTerrain(clique_indexes(:,2));
        %building_building = ~isTerrain(clique_indexes(:,1)) & ~isTerrain(clique_indexes(:,2));
        E = exp(-c*(del_DSM/del_0).^2);
        E_clique(terrain_terrain) = 1-E(terrain_terrain);
        E_clique(terrain_building | building_terrain) = E(terrain_building | building_terrain);
        E_clique( ...
            (terrain_building & (del_DSM > 0)) | (building_terrain & (del_DSM < 0)) ...
        ) = 1;
        E_clique = (1-coherence_const)*E_clique/length(steps); % divide by pixels in clique to make unary and interaction potentials equal
        
        isTerrain_past = isTerrain;
        isTerrain = zeros(rows,cols);
        
        N = numel(DSM);
        index_s = [clique_indexes(:,1); (N+1)*ones(N,1); (1:N)']; % The Nth+1 and Nth+2 nodes correspond to the source and sink nodes.
        index_t = [clique_indexes(:,2); (1:N)'; (N+2)*ones(N,1)];
        weights = [E_clique; E_data_terrain(:); E_data_building(:)]; 
        A = sparse(index_s, index_t, weights, N+2, N+2);
        G = digraph(A);

        [~,~,~,ct] = maxflow(G,N+1,N+2);
        isTerrain(ct(1:end-1)) = 1; % disredard actual sink node, which is the last element

        
        convergence_curr = 100*mean(isTerrain_past(:)==isTerrain(:));
        disp(['Iteration ', num2str(iteration), ' of ', num2str(max_iterations), ...
            ': ', num2str(convergence_curr), '% converged.']);
        
        if (abs(convergence_curr-convergence_prev)<0.05)
            disp('Algorithm has finished converging.')
            break;
        end
        convergence_prev = convergence_curr;

    end
    
end

% BuildCliqueIndexes builds out an Nx2 array of clique indexes. The first
% and second column correspond to connection origins and terminations
% respectively. The values in the matrix are the linear indexes of the
% node. E.g., a row in clique_indexes with [1,2] means that the the pixel
% in the 1st row 1st column is connected to the pixel in the 1st row 2nd 
% column. 
% Inputs: 
% rows and cols are the numbers of rows and columns in the DSM. 
% steps specifies how the cliques are built. E.g., steps =
% [[-1,1];[0,1];[1,1];[-1,0];[1,0];[-1,-1];[0,-1];[1,-1]]; indicates a
% second order clique in which pixels are connected to their eight nearest
% neighbors.
function clique_indexes = BuildCliqueIndexes(rows,cols,steps)
    clique_indexes = []; 
    [Rows,Cols] = meshgrid(1:rows,1:cols);
    Rows = Rows(:);
    Cols = Cols(:);
    for i=1:size(steps,1)
        clique_indexes = [clique_indexes; [Rows, Cols, Rows+steps(i,1), Cols+steps(i,2)]];
    end
    
    % remove connections that fell off of the grid
    offGrid = (clique_indexes(:,3) < 1)    | (clique_indexes(:,4) < 1)   | ...
              (clique_indexes(:,3) > rows) | (clique_indexes(:,4) > cols);
    clique_indexes(offGrid,:) = [];
    
    % convert subscripts to linear indexes
    linearIdx_s = sub2ind([rows,cols],clique_indexes(:,1),clique_indexes(:,2));
    linearIdx_t = sub2ind([rows,cols],clique_indexes(:,3),clique_indexes(:,4));
    clique_indexes = [linearIdx_s,linearIdx_t];

end
