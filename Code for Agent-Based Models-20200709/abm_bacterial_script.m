%% bacterial agent - standard life
M = randi(2,50,50)-1;
ndie_low = 2;
ndie_high = 3;
nborn = 3;
count = ABM_bacterial_agent(50,80,0.1,M,1,ndie_low,ndie_high,nborn);
figure()
scatter(1:80,count')

%% bacterial agent - reproductive 1
M = zeros(50);
M(25,24) = 1;
ngen = 200;
ndie_low = 2;
ndie_high = 3;
nborn = 1;
count = ABM_bacterial_agent(50,ngen,0.1,M,1,ndie_low,ndie_high,nborn);
%fig = figure('Visible','off','PaperPosition',[0 0 6 4],'PaperSize',[6 4]);
figure();
scatter(1:ngen,count')
%print(fig,'-dpdf','figures/agent_reproductive1.pdf')

%% bacterial agent - reproductive 2
M = zeros(50);
M(15,24) = 1;
M(24,15) = 1;
ngen = 100;
ndie_low = 1;
ndie_high = 4;
nborn = 1;
count = ABM_bacterial_agent(50,ngen,0.1,M,1,ndie_low,ndie_high,nborn);
figure()
scatter(1:ngen,count')

%% bacterial agent - reproductive 2 with random deaths
M = zeros(50);
M(15,24) = 1;
M(24,15) = 1;
ngen = 100;
ndie_low = 1;
ndie_high = 4;
nborn = 1;
count = ABM_bacterial_agent_death(50,ngen,0.1,M,1,ndie_low,ndie_high,nborn,0.3);
figure()
scatter(1:ngen,count')