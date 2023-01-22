% Improved Marine Predators Algorithm source code (Developed in MATLAB R2019b)
%
%  programming: Xunbu Zhai 
 
% E-mails: zxbsyuct@163.com
_______________________________________________________________________
%  Marine Predators Algorithm source code (Developed in MATLAB R2015a)
%
%  programming: Afshin Faramarzi & Seyedali Mirjalili
%
% paper:
%  A. Faramarzi, M. Heidarinejad, S. Mirjalili, A.H. Gandomi, 
%  Marine Predators Algorithm: A Nature-inspired Metaheuristic
%  EPyrndert Systems with Applications
%  DOI: doi.org/10.1016/j.eswa.2020.113377
%  
%  E-mails: afaramar@hawk.iit.edu            (Afshin Faramarzi)
%           muh182@iit.edu                   (Mohammad Heidarinejad)
%           ali.mirjalili@laureate.edu.au    (Seyedali Mirjalili) 
%           gandomi@uts.edu.au               (Amir H Gandomi)
%_________________________________________________________________________

function [Top_predator_fit,Top_predator_pos,Convergence_curve]=Proposed_MPA(SearchAgents_no,Max_iter,lb,ub,dim,fobj )

%%
  

Top_predator_pos=zeros(1,dim);
Top_predator_fit=inf; 

Convergence_curve=zeros(1,Max_iter);
stepsize=zeros(SearchAgents_no,dim);
fitness=inf(SearchAgents_no,1);

Py=initialization(SearchAgents_no,dim,ub,lb);
  
Xmin=repmat(ones(1,dim).*lb,SearchAgents_no,1);
Xmax=repmat(ones(1,dim).*ub,SearchAgents_no,1);
         
newNopar=0.5;
Iter=0;
FADs=0.2;
P=0.5;

while Iter<Max_iter    
    

    
     %------------------- Detecting top predator -----------------    
 for i=1:size(Py,1)  
 
    Flag4ub=Py(i,:)>ub;
    Flag4lb=Py(i,:)<lb;    
    Py(i,:)=(Py(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb ;                   

     fitness(i,1)= - fobj(Py(i,:),Iter,Max_iter);%%%
                     
     if fitness(i,1)< Top_predator_fit 
       Top_predator_fit=fitness(i,1); 
       Top_predator_pos=Py(i,:);
     end          
 end
     
     %------------------- Marine Memory saving ------------------- 
    
 if Iter==0
   fit_old=fitness;    Py_old=Py;
 end
     
  Inx=(fit_old<fitness);
  Indx=repmat(Inx,1,dim);
  Py=Indx.*Py_old+~Indx.*Py;
  fitness=Inx.*fit_old+~Inx.*fitness;
        
  fit_old=fitness;    Py_old=Py;

     %------------------------------------------------------------   
     
 Elite=repmat(Top_predator_pos,SearchAgents_no,1);  

h=(1-(Iter./Max_iter).^3).^2;
CF=abs(h.*sin(pi*2/3+sin(h.*pi)));%%NCF Equ.10
                             
 RL=0.05*levyGen(SearchAgents_no,dim,1.5);   %Levy random number vector
 RB=randn(SearchAgents_no,dim);          %Brownian random number vector
           
  for i=1:size(Py,1)
     for j=1:size(Py,2)        
       R=rand();
          %------------------ Phase 1 (Equ.5) ------------------- 
       if Iter<Max_iter/3 
          stepsize(i,j)=RB(i,j)*(Elite(i,j)-RB(i,j)*Py(i,j));                    
          Py(i,j)=Py(i,j)+P*R*stepsize(i,j); 
             
          %--------------- Phase 2 (Eqs. 6 & 7)----------------
       elseif Iter>Max_iter/3 && Iter<2*Max_iter/3 
          
         if i>size(Py,1)/2
            stepsize(i,j)=RB(i,j)*(RB(i,j)*Elite(i,j)-Py(i,j));
            Py(i,j)=Elite(i,j)+P*CF*stepsize(i,j); 
         else
            stepsize(i,j)=RL(i,j)*(Elite(i,j)-RL(i,j)*Py(i,j));                     
            Py(i,j)=Py(i,j)+P*R*stepsize(i,j);  
         end  
         
         %----------------- Phase 3 (Equ. 8)-------------------
       else 
           
           stepsize(i,j)=RL(i,j)*(RL(i,j)*Elite(i,j)-Py(i,j)); 
           Py(i,j)=Elite(i,j)+P*CF*stepsize(i,j);  
    
       end  
      end                                         
  end    
  
  %------------------ Detecting top predator ------------------        
		  for i=1:size(Py,1)  
			
			Flag4ub=Py(i,:)>ub;  
			Flag4lb=Py(i,:)<lb;  
			Py(i,:)=(Py(i,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
		  
			  fitness(i,1)= - fobj(Py(i,:),Iter,Max_iter);
				
			  if fitness(i,1)<Top_predator_fit 
				 Top_predator_fit=fitness(i,1);
				 Top_predator_pos=Py(i,:);
			  end     
		  end
				
			 %---------------------- Marine Memory saving ----------------
		 
			 if Iter==0
				fit_old=fitness;    Py_old=Py;
			 end
				 
				Inx=(fit_old<fitness);
				Indx=repmat(Inx,1,dim);
				Py=Indx.*Py_old+~Indx.*Py;
				fitness=Inx.*fit_old+~Inx.*fitness;
					
				fit_old=fitness;    Py_old=Py;
	 	
				
 
     %---------- Eddy formation and FADs effect (Equ. 9) ----------- 
                             
  if rand()<FADs
     U=rand(SearchAgents_no,dim)<FADs;                                                                                              
     Py=Py+CF*((Xmin+rand(SearchAgents_no,dim).*(Xmax-Xmin)).*U);

  else
     r=rand();  Rs=size(Py,1);
     stepsize=(FADs*(1-r)+r)*(Py(randperm(Rs),:)-Py(randperm(Rs),:));
     Py=Py+stepsize;
  end
                                                        
  Iter=Iter+1;  
  Convergence_curve(Iter)= -Top_predator_fit; 
       
% 
	 	
	solutionsX= Py ;  
	solutionsY= fitness;		
			 
	newSolutions=[];
	orderY=[];
	togetherY=solutionsY;
	together=solutionsX;
											 
	if size(togetherY,1 )==1						 
        togetherY=togetherY'; 
	end
	  
	[~,I]=sort(togetherY ,'ascend');   
	togetherS=together(I,:); 

    Pymin=min(togetherS);
    Pymax=max(togetherS);  
 			
	Best_X=togetherS(1,:);		
	numV = length(Best_X);
	  	
    dm=sqrt(mahal(Py(:),Py(:))); %% Equ.11
    Pyw = max(dm);
    Pyb = min(dm);
	 
	for lec=1: ceil(length(solutionsY)*newNopar)    	
		
		Pym = mean(togetherS);    % Equ.12              
        	
        Nrd=(1-(Iter/Max_iter)) ;
        std=normrnd(0 ,Nrd);
        Best_X = min(Best_X, Pyb);
		Pyc =    Best_X +std*(Best_X-Pym); 		% Equ.13
		PyLEO = Pyc;
        
        if rand < .66      
			 			  
				Pyrnd =  Pymin+(Pymax-Pymin).*rand(1,numV);  % Equ.14
				Py1c=Pyc;
                Py2c= rand*Py1c+(1-rand)*Pyrnd+ std *(Best_X-Pyw); % Equ.15
			  
             if rand < 0.5						 
										 
					PyLEO =   Py1c+(Best_X-Pyrnd)+rand*(Py2c-Py1c) ;  %		Equ.16		 
				else
					 
					PyLEO =   Best_X+(Best_X-Pyrnd)+rand*(Py2c-Py1c) ;	%	Equ.17	
             end			
        end
        		 
		 newSolutions =[newSolutions; PyLEO ];
				
	    if size(newSolutions,1)>= ceil(length(solutionsY)*newNopar) 	
			 break; 
	    end
	end  	
 
orderY= flip(I) ;
               
		for noN=1: size(newSolutions,1) 
			whichOne=orderY(noN) ;	
			whichOneX=newSolutions(noN,:);			
							
			Py(whichOne,:)=whichOneX;  	
			
			Flag4ub=Py(whichOne,:)>ub;     
			Flag4lb=Py(whichOne,:)<lb;  
			Py(whichOne,:)=(Py(whichOne,:).*(~(Flag4ub+Flag4lb)))+ub.*Flag4ub+lb.*Flag4lb;
			 
		end	     	 
	   
end 

