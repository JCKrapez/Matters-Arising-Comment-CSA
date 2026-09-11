% Chameleon_VMA.m

% The first version of this function is
% Chameleon Swarm Algorithm (CSA) source codes version 1.0
%
%  Developed in MATLAB R2018a
%
%  Author and programmer: Malik Braik
%
%         e-Mail: m_fjo@yahoo.com
%                 mbraik@bau.edu.au

%   Main paper:
%   Malik Sh. Braik,
%   Chameleon Swarm Algorithm: A Bio-inspired Optimizer for Solving Engineering Design Problems
%   Expert Systems with Applications
%   DOI: https://doi.org/10.1016/j.eswa.2021.114685


% A version Chameleon.m
% was obtained from https://github.com/KirinShi/Temperature-reconstruction-using-CSA
% He, Y., Chen, M. K., Huang, M. et al. Dispersive Meta-lens Thermometry for High-temperature Measurements.
% Nature Communications 16, 10090 (2025). https://doi.org/10.1038/s41467-025-65171-7

% The analysis of the paper, next the code and data made available by the
% authors, led to the submission to Nature Communications of a document "Matters Arising"

% The changes made to the function that led to Chameleon_VMA.m made by:
% J.-C. Krapez
% 2026-09-07
% DOTA, ONERA, 13300 Salon de Provence, France
% jean-claude.krapez@onera.fr

% VMA is the flag used to highlight the changes (VMA for Version Matters
% Arising):
% - modifications to enable intermediate plots
% - modification to introduce an offset to "time" used in the formula of
% the acceleration alpha and thus avoid a division by 0 at the first
% iteration.
% - introduction of global variables get information on variables necessary for the plots

function [fmin0,gPosition,cg_curve]=Chameleon_VMA(searchAgents,iteMax,lb,ub,dim,fobj)
global ii jj % VMA  the coordinates ii and jj will be added in the plots created in Chameleon
global stdTr  % added variable containing the best fit, i.e. lowest std of Tr (fmin0) which was surprisingly overwritten by the last instruction in the He version of Chameleon.
global avis_seed seed %  seed used to initiate the random numbers 
global t_threshold  %  introduction of a threshold to avoid that, at the first iteration, the denominator is 0 when calculating the velocity updates  !!!!
global avis_plot avis_bb_ceram
if(avis_bb_ceram==0)
lambdalist = (0.65:-0.01:0.45);  % VMA lambdalist was recalled for the plots added in Chameleon
elseif(avis_bb_ceram==1)
lambdalist = (0.45:0.01:0.65);  % VMA lambdalist was recalled for the plots added in Chameleon
end

%%%%* 1
if size(ub,2)==1
    ub=ones(1,dim)*ub;
    lb=ones(1,dim)*lb;
end

% if size(ub,1)==1
%     ub=ones(dim,1)*ub;
%     lb=ones(dim,1)*lb;
% end



% VMA : seed for reproductibility of number generation
if(avis_seed>1)
    rng(seed);
end
%  disp([' seed : ' num2str(seed)])
%  pause

%% Convergence curve
cg_curve=zeros(1,iteMax);
%
% f1 =  figure (1);
% set(gcf,'color','w');
% hold on
% xlabel('Iteration','interpreter','latex','FontName','Times','fontsize',10)
% ylabel('fitness value','interpreter','latex','FontName','Times','fontsize',10);
% grid;



%% Initial population

chameleonPositions=initialization(searchAgents,dim,ub,lb);% Generation of initial solutions

%% Evaluate the fitness of the initial population
    if(mod(avis_plot,2)==0) %VMA
disp(' After initialization:'); %VMA
    end
fit=zeros(searchAgents,1);

for i=1:searchAgents
    fit(i,1)=fobj(chameleonPositions(i,:));
    if(mod(avis_plot,2)==0) %VMA
        figure(60);clf;%VMA
        plot(lambdalist,(chameleonPositions(i,:)),'k');%hold on; %VMA
        title({'Log(emissivity)'},'interpreter','latex'); % VMA
        ylabel('ln(emissivity)');  % VMA
        
        figure(61);clf;%VMA
        plot(lambdalist,(exp(chameleonPositions(i,:))),'k');%hold on; %VMA
        title({'Emissivity'},'interpreter','latex'); % VMA
        ylabel('Emissivity');  % VMA
        disp(['agent: ' num2str(i) ', cost function: ' num2str( fit(i,1)) ]); %VMA
        pause % VMA
    end
end

%% Initalize the parameters of CSA
fitness=fit; % Initial fitness of the random positions

[fmin0,index]=min(fit);
    if(mod(avis_plot,2)==0) %VMA
disp(['min cost function : ' num2str(fmin0) ]); %VMA
disp(['best agent : ' num2str(index) ]); %VMA
    end
chameleonBestPosition = chameleonPositions; % Best position initialization
gPosition = chameleonPositions(index,:); % initial global position

if(mod(avis_plot,2)==0) %VMA
    figure(60);clf;
    plot(lambdalist,(gPosition),'b','LineWidth',3);%hold on; %VMA
    ylabel('ln(emissivity)');  % VMA
    figure(61);clf;
    plot(lambdalist,exp(gPosition),'b','LineWidth',3);%hold on; %VMA
    ylabel('emissivity');  % VMA
    pause;%VMA
end

v=0.1*chameleonBestPosition;% initial velocity

v0=0.0*v;

%% Start CSA
% Main parameters of CSA according to He
rho=1.0;
p1=0.25;
p2=1.5;
c1=2.0;
c2=1.80;
% gamma=2;
gamma=1;
alpha = 2.5;%原来为3.5
beta=3.0; %yuanlai 3

% % % VMA: 
% % % Main parameters of CSA according to Braik:
% % rho=1.0;
% % p1=2.0;  
% % p2=2.0;  
% % c1=2.0; 
% % c2=1.80;  
% % gamma=2.0; 
% % alpha = 4.0;  
% % beta=3.0; 


%% Start CSA
for t=1:iteMax
    %    a = 2590*(1-exp(-log(t)));    %2590*(1-exp(-log(t)));
    %    a = 2590*(1-exp(-log(t+t_offset)));   %VMA:  addition of a time offset to avoid that, at the first iteration, the denominator is 0 when calculating the velocity updates  !!!!
    a = 2590*(1-exp(-log(max(t,t_threshold))));   %VMA:  introduction of a threshold to avoid that, at the first iteration, the denominator is 0 when calculating the velocity updates  !!!!
    
    omega=(1-(t/iteMax))^(rho*sqrt(t/iteMax)) ;
    p1 = 2* exp(-2*(t/iteMax)^2);  %
    p2 = 2/(1+exp((-t+iteMax/2)/100)) ;
    
    mu= gamma*exp(-(alpha*t/iteMax)^beta) ;
    
    ch=ceil(searchAgents*rand(1,searchAgents));
    
    %VMA
    %disp(' After Update the position of CSA (Exploration) '); %VMA
    
    if(mod(avis_plot,3)==0) %VMA
        disp(' Exploration: '); %VMA
    end
    %% Update the position of CSA (Exploration)
    for i=1:searchAgents
        
        if(mod(avis_plot,3)==0) %VMA
            disp([' agent: ' num2str(i) ', before Exploration:']); %VMA
            %chameleonPositions(i,:) %VMA
            figure(60);clf; %VMA
            plot(lambdalist,(chameleonPositions(i,:)),'k');hold on; %VMA
            ylabel('ln(emissivity)');  % VMA
            figure(61);clf; %VMA
            plot(lambdalist,(exp(chameleonPositions(i,:))),'k');hold on; %VMA
            ylabel('emissivity');  % VMA
        end
        if rand>=0.1
            chameleonPositions(i,:)= chameleonPositions(i,:)+ p1*(chameleonBestPosition(ch(i),:)-chameleonPositions(i,:))*rand()+...
                + p2*(gPosition -chameleonPositions(i,:))*rand();
        else
            for j=1:dim
                chameleonPositions(i,j)=   gPosition(j)+mu*((ub(j)-lb(j))*rand+lb(j))*sign(rand-0.50) ;
            end
        end
        if(mod(avis_plot,3)==0) %VMA
            disp([' agent: ' num2str(i) ', after Exploration:']); %VMA
            %chameleonPositions(i,:) %VMA
            figure(60); %VMA
            plot(lambdalist,(chameleonPositions(i,:)),'r');hold on; %VMA
            ylabel('ln(emissivity)');  % VMA
            figure(61); %VMA
            plot(lambdalist,(exp(chameleonPositions(i,:))),'r');hold on; %VMA
            ylabel('emissivity');  % VMA
            pause; %VMA
        end
    end
    %% Rotation of the chameleons - Update the position of CSA (Exploitation)
    
    %%% Rotation 180 degrees in both direction or 180 in each direction
    %
    if(mod(avis_plot,5)==0) %VMA
        disp(' Rotation: '); %VMA
    end
    [chameleonPositions] = rotation(chameleonPositions, searchAgents, dim);
    if(mod(avis_plot,5)==0) %VMA
        for i=1:searchAgents
            disp([' agent: ' num2str(i) ', after Rotation:' ]); %VMA
            %chameleonPositions(i,:) %VMA
            figure(60);clf; %VMA
            plot(lambdalist,(chameleonPositions(i,:)),'m');%hold on; %VMA
            ylabel('ln(emissivity)');  % VMA
            figure(61);clf; %VMA
            plot(lambdalist,exp(chameleonPositions(i,:)),'m');%hold on; %VMA
            ylabel('emissivity');  % VMA
            pause; %VMA
        end
    end
    %VMA
    
    if(mod(avis_plot,7)==0) %VMA
        disp(' Velocity updates:  '); %VMA
    end
    %%  % Chameleon velocity updates and find a food source
    for i=1:searchAgents
        
        v(i,:)= omega*v(i,:)+ p1*(chameleonBestPosition(i,:)-chameleonPositions(i,:))*rand +....
            + p2*(gPosition-chameleonPositions(i,:))*rand;
        
        
        if(mod(avis_plot,7)==0) %VMA
            figure(60);clf; %VMA
            %chameleonPositions(i,:) %VMA
            plot(lambdalist,(chameleonPositions(i,:)),'m');hold on; %VMA
            ylabel('ln(emissivity)');  % VMA
            figure(61);clf; %VMA
            plot(lambdalist,exp(chameleonPositions(i,:)),'m');hold on; %VMA
            ylabel('emissivity');  % VMA
        end
        
        chameleonPositions(i,:)=chameleonPositions(i,:)+(v(i,:).^2 - v0(i,:).^2)/(2*a);
        
        if(mod(avis_plot,7)==0) %VMA
            disp([' acceleration a : ' num2str(a) ' agent: ' num2str(i) ', before (m) and after (c) velocity updates:' ]); %VMA
            %chameleonPositions(i,:) %VMA
            figure(60); %VMA
            plot(lambdalist,(chameleonPositions(i,:)),'c');%hold on; %VMA
            ylabel('ln(emissivity)');  % VMA
            figure(61); %VMA
            plot(lambdalist,exp(chameleonPositions(i,:)),'c');%hold on; %VMA
            ylabel('emissivity');  % VMA
            pause; %VMA
        end
        
        
    end
    
    
    
    v0=v;
    
    
    %% handling boundary violations
    if(mod(avis_plot,11)==0) %VMA
        disp(' HANDLING boundary violations: '); %VMA
    end
    for i=1:searchAgents
        if(mod(avis_plot,11)==0) %VMA
            disp([' agent: ' num2str(i) ', before 1st handling boundary violations: ' ]); %VMA
            %chameleonPositions(i,:) %VMA
            figure(60);clf;%VMA
            plot(lambdalist,chameleonPositions(i,:),'k');hold on; %VMA
            ylabel('ln(emissivity)');  % VMA
            figure(61);clf;%VMA
            plot(lambdalist,exp(chameleonPositions(i,:)),'k');hold on; %VMA
            ylabel('emissivity');  % VMA
        end
        % VMA Remplacement par le VECTEUR lb si TOUTES les valeurs < lb
        % VMA Remplacement par le VECTEUR ub si TOUTES les valeurs > ub
        
        if chameleonPositions(i,:)<lb
            chameleonPositions(i,:)=lb;
        elseif chameleonPositions(i,:)>ub
            chameleonPositions(i,:)=ub;
        end
        
        % VMA   chameleonPositions(i,:)=max(chameleonPositions(i,:),lb);
        % VMA   chameleonPositions(i,:)=min(chameleonPositions(i,:),ub);
        
        if(mod(avis_plot,11)==0) %VMA
            disp([' agent: ' num2str(i) ', after 1st handling boundary violations: ' ]); %VMA
            %chameleonPositions(i,:) %VMA
            figure(60);%VMA
            plot(lambdalist,chameleonPositions(i,:),'b');hold on; %VMA
            ylabel('ln(emissivity)');  % VMA
            figure(61);%VMA
            plot(lambdalist,exp(chameleonPositions(i,:)),'b');hold on; %VMA
            ylabel('emissivity');  % VMA
            pause; %VMA
        end %VMA
    end
    
    %%
    
    
    for i=1:searchAgents
        if(mod(avis_plot,11)==0) %VMA
            disp([' agent: ' num2str(i) ', before 2nd handling boundary violations: ' ]); %VMA
            %chameleonPositions(i,:) %VMA
            figure(60);clf;%VMA
            plot(lambdalist,(chameleonPositions(i,:)),'b');hold on; %VMA
            ylabel('ln(emissivity)');  % VMA
            figure(61);clf;%VMA
            plot(lambdalist,exp(chameleonPositions(i,:)),'b');hold on; %VMA
            ylabel('emissivity');  % VMA
        end
        
        ub_=sign(chameleonPositions(i,:)-ub)>0;
        lb_=sign(chameleonPositions(i,:)-lb)<0;
        
        chameleonPositions(i,:)=(chameleonPositions(i,:).*(~xor(lb_,ub_)))+ub.*ub_+lb.*lb_;  %%%%%*2
        
        if(mod(avis_plot,11)==0) %VMA
            disp([' agent: ' num2str(i) ', after 2nd handling boundary violations: ' ]); %VMA
            %chameleonPositions(i,:) %VMA
            figure(60);%VMA
            plot(lambdalist,(chameleonPositions(i,:)),'r');hold on; %VMA
            ylabel('ln(emissivity)');  % VMA
            figure(61);%VMA
            plot(lambdalist,exp(chameleonPositions(i,:)),'r');hold on; %VMA
            ylabel('emissivity');  % VMA
            pause; %VMA
        end
        
        fit(i,1)=fobj (chameleonPositions(i,:)) ;
        %  disp(['agent ' num2str(i) ', cost function : ' num2str( fit(i,1))]); %VMA
        if fit(i)<fitness(i)
            if(mod(avis_plot,11)==0) %VMA
                disp(['agent: ' num2str(i) ', cost function: ' num2str( fit(i,1)) ' is lower than ' num2str(fitness(i)) ' by ' num2str(fitness(i)-fit(i,1)) ]); %VMA
                figure(60);
                plot(lambdalist,(chameleonPositions(i,:)),'g','linewidth',2);%hold on; %VMA
                ylabel('ln(emissivity)');  % VMA
                figure(61);
                plot(lambdalist,exp(chameleonPositions(i,:)),'g','linewidth',2);%hold on; %VMA
                ylabel('emissivity');  % VMA
                pause; %VMA
            end
            chameleonBestPosition(i,:) = chameleonPositions(i,:); % Update the best positions
            fitness(i)=fit(i); % Update the fitness
        end
        
    end
    
    
    %% Evaluate the new positions
    
    [fmin,index]=min(fitness); % finding out the best positions
    %disp(['min cost function : ' num2str(fmin) ]); %VMA
    %disp(['best agent : ' num2str(index) ]); %VMA
    
    
    % Updating gPosition and best fitness
    if fmin < fmin0
        gPosition = chameleonBestPosition(index,:); % Update the global best positions
        fmin0 = fmin;
    end
    
    if((avis_plot==1 && t==iteMax)|| avis_plot>1)
        figure(60);
        plot(lambdalist,(gPosition),'g','linewidth',3);%hold on; %VMA
        ylabel('ln(emissivity)');  % VMA
        figure(61);
        plot(lambdalist,exp(gPosition),'g','linewidth',3);%hold on; %VMA
        ylabel('emissivity');  % VMA
        %pause(0.4);%VMA
        
    end
    %% Print the results
    %   outmsg = ['Iteration# ', num2str(t) , '  Fitness= ' , num2str(fmin0)];
    %   disp(outmsg);
    
    % VMA : j'ai décommenté:
    
    %% Visualize the results
    
    cg_curve(t)=fmin0; % Best found value until iteration t
    
    %     if t>2
    %      set(0, 'CurrentFigure', f1)
    %
    %         line([t-1 t], [cg_curve(t-1) cg_curve(t)],'Color','b');
    %         title({'Convergence characteristic curve'},'interpreter','latex','FontName','Times','fontsize',12);
    %         xlabel('Iteration');
    %         ylabel('Best score obtained so far');
    %         drawnow
    %     end
    
    % VMA : uncommented previous 9 lines, starting with t=2
    %if t>2
    if t>1
        %VMA set(0, 'CurrentFigure', f1)
        figure(82);
        if(t==2)
            clf; % VMA
        end
        line([t-1 t], [cg_curve(t-1) cg_curve(t)],'Color','b');
        title({'Convergence characteristic curve'},'interpreter','latex','FontName','Times','fontsize',12);
        xlabel('Iteration');
        %VMA ylabel('Best score obtained so far');
        ylabel('Best score obtained so far (std in K)');  % VMA
        drawnow
    end
end



% VMA : The next two lines in Braik and He versions are useless. The best
% position is already known, it is  gPosition. The best fit is known, it is
% fmin0.
ngPosition=find(fitness==min(fitness));
g_best=chameleonBestPosition(ngPosition(1),:);  % Solutin of the problem

stdTr=fmin0;  % VMA Instruction added to store the general best fit in stdTr which is then used outside Chameleon without alteration.

% VMA : The last line in Braik version :
% fmin0 =fobj(g_best);
% is also useless, since the best general fit fmin0 is already known, it has been updated after each iteration.
% disp(['fmin0 after fobj(g_best) :' num2str(fmin0)]);  % VMA

% VMA : The last line in He version was that of Braik, yet changed into :
fmin0 =g_best;
% It is now even counterproductive.
% fmin0, which contained the best fit (std of radiance temperature), has been overwritten with
% g_best, which is the best position
% Why fmin0 changes its status?

% VMA: Instruction that were used for verification :
% disp(['fmin0 after the iterations :' num2str(fmin0)]);  % VMA
% disp(['fitness after the iterations :']);  % VMA
% fitness   % VMA
% disp(['gPosition after the iterations :']);  % VMA
% gPosition
% % disp(['best agent was AGAIN CALCULATED :' num2str(ngPosition)]);  % VMA
% disp(['best position was AGAIN CALCULATED :' num2str(g_best) ]);  % VMA
% pause  % VMA

end