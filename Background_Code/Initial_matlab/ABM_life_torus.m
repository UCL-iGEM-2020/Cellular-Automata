% Life: agent based model
% Uses toroidal boundary conditions
% To run use, e.g. life_torus(30,100,0.01,randi(2,30)-1);

function ABM_life_torus(grid,ngen,pause_time,M)

s = size(M);

if s(1) ~= s(2) || s(1) ~= grid
    disp('input matrix M does not have dimensions grid x grid');
    return 
end

nneigh = zeros(grid,grid);

figure()
colormap(gray(2));
image(2*(1-M) );
pause(pause_time);
for g=1:ngen
    %disp('M')
    %disp(M)
    
    % count the nearest neighbours
    for i=1:grid
        for j=1:grid
           count = 0;
           %disp([i j])
           for x=i-1:i+1
               for y=j-1:j+1
                   %if x >= 1 && y >= 1 && x <= grid && y <= grid
                   if 1
                        if x == 0
                            x = grid;
                        end
                        if x == grid+1
                            x = 1;
                        end
                        if y == 0
                            y = grid;
                        end
                        if y == grid+1
                            y = 1;
                        end
                        if x == i && y == j
                            % dont count this element
                        else
                            count = count + M(x,y);
                            %disp([ x y count])
                        end
                   end
               end
           end
        
           nneigh(i,j) = count;   
        end
    end
    %disp('NN')
    %disp(nneigh)
    
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
    image(2*(1-M));
    pause(pause_time);
end

disp('Terminated')
end
