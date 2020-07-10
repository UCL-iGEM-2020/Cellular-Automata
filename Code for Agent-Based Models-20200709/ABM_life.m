%% life agent based model
function ABM_life(grid,ngen,pause_time)

% initialise agent and number_of_neighbours matrices
M = randi(2,grid,grid)-1;
nneigh = zeros(grid,grid);

%initialise matrix to store number of active cells
activity = zeros(1,ngen);

% plot the initial agent matrix
figure()
colormap(gray(2));
image(2*(1-M));

% initialise plot for game of life
figure()
colormap(gray(2));

% loop over number of generations
for g=1:ngen    
    % Loop to count the nearest neighbours
    for i=1:grid               % go through rows
        for j=1:grid           % go through columns
           count = 0;
           %disp([i j])
           for x=i-1:i+1        % go through horizontal neighbours
               for y=j-1:j+1    % go through horizontal neighbours
                   
                   if x >= 1 && y >= 1 && x <= grid && y <= grid   % correct for margins
                        if x == i && y == j
                            % dont count the current element
                        else
                            count = count + M(x,y);   % add 1 to count if neighbour is 1
                            %disp([ x y count])
                        end
                   end
               end
           end
        
           nneigh(i,j) = count; % update number_of_neighbours matrix
        end
    end
    
    % Show the result of neighbour counting in the command window
 %   disp('NN');
 %   disp(nneigh);
    
    % update the grid
    for i=1:grid
        for j=1:grid
            if M(i,j) == 1
                % see if this cell lives or dies
                if nneigh(i,j) < 2 || nneigh(i,j) > 3
                    M(i,j) = 0;
                end
            else
                % see if this cell is born
                if nneigh(i,j) == 3
                    M(i,j) = 1;
                end
            end
        end
    end
    activity(g) = sum(sum(M));
    
    % plot updated agent matrix
    image(2*(1-M));
    pause(pause_time);
end
figure()
plot(activity)
end
