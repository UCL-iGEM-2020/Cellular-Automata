% run_life
ABM_life(50,50,0.1)

%% Blinker 	
M = zeros(50);
M(5,4:6)=1;
ABM_life_mod(50,50,0.1,M)

%% Toad
M = zeros(50);
M(5,4:6)=1; M(6,3:5)=1;
ABM_life_mod(50,50,0.1,M)

%% Glider
M = zeros(50);
M(3,1:3)=1;M(2,3)=1;M(1,2)=1;
ABM_life_mod(50,50,0.1,M)

%% Spaceship
M = zeros(50);
%M(5,5:8)=1;M(4,9)=1;M(3:4,5)=1; M(2,6)=1;M(2,9)=1;
M(25,25:28)=1;M(24,29)=1;M(23:24,25)=1; M(22,26)=1;M(22,29)=1;
ABM_life_mod(50,80,0.1,M)

%% Spaceship shifted
M = zeros(50);
M(25,27:30)=1;M(24,31)=1;M(23:24,27)=1; M(22,28)=1;M(22,31)=1;
ABM_life_mod(50,80,0.1,M)

%% Diehard
M = zeros(50);
M(10,10:11)=1;M(11,11)=1; M(11,15:17)=1;M(9,16)=1;
ABM_life_mod(50,150,0.1,M)

%% Spaceship (with toroidal boundary conditions)
M = zeros(50);
%M(5,5:8)=1;M(4,9)=1;M(3:4,5)=1; M(2,6)=1;M(2,9)=1;
M(25,25:28)=1;M(24,29)=1;M(23:24,25)=1; M(22,26)=1;M(22,29)=1;
ABM_life_torus(50,80,0.1,M)

%% random initial conditions
M = randi(2,50,50)-1;
ABM_life_torus(50,80,0.1,M)
