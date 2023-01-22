clear all;
clc
format long
SearchAgents_no= 30 ; % Number of search agents
% Dim = 50;
Function_name='F1';  
Max_iteration=500; % Maximum number of iterations
[lb,ub,dim,fobj]=Get_Functions_details_B5G(Function_name);
diary on
 [Best_score, Best_pos, Convergence_curve ]=Proposed_MPA(SearchAgents_no,Max_iteration,lb,ub,dim,fobj);
% display(['The best solution obtained by Proposed_MPA is \m ', num2str(Best_pos)]);
display(['The best optimal value of the objective function found by Proposed_MPA is : ', num2str(-Best_score,10)]);

diary off

 figure('Position',[80 80 280 230])
 semilogy(Convergence_curve(1:500),'r-','Markersize',1)%
set(gca,'yminortick','off')
% set(gca,'XTicklabel',{'0','100','200','300','400','500'});
 xlabel('\fontname{Times New Roman}\fontsize{7.5}Iterations')
 ylabel('\fontname{Times New Roman}\fontsize{7.5}Fitness value')
legend('\fontname{Times New Roman}\fontsize{7.5}Proposed')


