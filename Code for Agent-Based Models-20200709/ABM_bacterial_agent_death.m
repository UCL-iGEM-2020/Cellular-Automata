% life agent based model
function agent_count=ABM_bacterial_agent_death(grid,ngen,pause_time,M,edge,ndie_low,ndie_high,nborn, p)

s = size(M);

if s(1) ~= s(2) || s(1) ~= grid
    disp('input matrix M does not have dimensions grid x grid');
    return 
end

nneigh = zeros(grid,grid);

agent_count = zeros(1,ngen);

figure()
colormap(gray(2));
image(2*(1-M) );
pause(pause_time);
for g=1:ngen
    %disp('M')
    %disp(M)
    agent_count(g) = sum(sum(M));
    
    % count the nearest neighbours
    for i=1:grid
        for j=1:grid
           count = 0;
           %disp([i j])
           for x=i-1:i+1
               for y=j-1:j+1
                   
                   if edge == 1
                        if x >= 1 && y >= 1 && x <= grid && y <= grid
                   
                            if x == i && y == j
                                % dont count this element
                            else
                                count = count + M(x,y);
                                %disp([ x y count])
                            end
                        end
                   else     
                        if x == 0
                            x = grid
                        end
                        if x == grid+1
                            x = 1
                        end
                        if y == 0
                            y = grid
                        end
                        if y == grid+1
                            y = 1
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
                if nneigh(i,j) < ndie_low || nneigh(i,j) > ndie_high
                    M(i,j) = 0;
                end
                
                
                
            else
                % see if this cell is born
                if nneigh(i,j) == nborn
                    M(i,j) = 1;
                end
            end
            
            % add random death
            u = rand();
            if u < p
                M(i,j) = 0;
            end
            
            
        end
    end
    image(2*(1-M));
    pause(pause_time);
end

end