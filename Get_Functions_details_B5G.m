%% ���Ժ���
function [lb, ub, dim, fobj] = Get_Functions_details_B5G(F)

switch F
   %% power allocation     
    case 'F1'
        fobj = @F1;
        %%
        lb = [ 0.0000001; 0.0000001; 0.0000001; 0.00000001; 0.0000001;0.00000001; 0.00000001; 0.00000001; 0.00000001; 0.000000001; ...
             0.00000001; 0.0000001; 0.0000001; 0.00000001; 0.0000001;0.0000001; 0.00000001; 0.00000001; 0.00000001; 0.000000001 ; ...
             0.00000001; 0.00000001; 0.00000001; 0.0000001; 0.00000001;0.00000001; 0.00000001; 0.0000001; 0.0000001; 0.00000001 ; ... 
             0.00000001; 0.00000001; 0.00000001; 0.0000001; 0.0000001;0.0000001; 0.0000001; 0.0000001; 0.0000001; 0.00000001 ; ...
             0.00000001; 0.00000001; 0.00000001; 0.0000001; 0.0000001;0.0000001; 0.0000001; 0.0000001; 0.0000001; 0.00000001 ]';
        ub = [0.0064; 0.0064; 0.0064; 0.0064; 0.0064; 0.0064;  0.0064; 0.0064; 0.0064;  0.0064; ...
              0.0064; 0.0064; 0.0064; 0.0064; 0.0064; 0.0064;  0.0064; 0.0064; 0.0064;  0.0064; ...
              0.0064; 0.0064; 0.0064; 0.0064; 0.0064; 0.0064;  0.0064; 0.0064; 0.0064;  0.0064 ; ...
              0.0064; 0.0064; 0.0064; 0.0064; 0.0064; 0.0064;  0.0064; 0.0064; 0.0064;  0.0064; ...
              0.0064; 0.0064; 0.0064; 0.0064; 0.0064; 0.0064;  0.0064; 0.0064; 0.0064;  0.0064 ]';
        dim = 50;     
end
function   PHI = F1(x,Iter,Max_iter)
% B5G

p1=x(1);p2=x(2);p3=x(3);p4=x(4);p5=x(5);
p6=x(6);p7=x(7);p8=x(8);p9=x(9);p10=x(10);

p11=x(11);p12=x(12);p13=x(13);p14=x(14);p15=x(15);
p16=x(16);p17=x(17);p18=x(18);p19=x(19);p20=x(20);
p21=x(21);p22=x(22);p23=x(23);p24=x(24);p25=x(25);
p26=x(26);p27=x(27);p28=x(28);p29=x(29);p30=x(30);
% % 
p31=x(31);p32=x(32);p33=x(33);p34=x(34);p35=x(35);
p36=x(36);p37=x(37);p38=x(38);p39=x(39);p40=x(40);
p41=x(41);p42=x(42);p43=x(43);p44=x(44);p45=x(45);
p46=x(46);p47=x(47);p48=x(48);p49=x(49);p50=x(50);

pm=[ p1 p2  p3  p4  p5  p6  p7  p8  p9  p10 ...
    p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 ...
    p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 ...
    p31 p32  p33  p34  p35  p36  p37  p38  p39  p40 ...
    p41 p42 p43 p44 p45 p46 p47 p48 p49 p50 ];

hm=[500 530 550 610 680 1000 1105 1270 1300 1310  ...
    1370 1500 1530 1650 1780 2000 2080 2200 2500 2710 ...
    2710 2820 2900 2950 2980 3000 3030 3040 3080 3120 ... 
   3210 3250 3340 3350 3450 3490 3610 3740 3940 4020 ...
    4050 4100 4150 4205 4380 5100 5405 5507 5700 5810 ];
A=20;
B=30;
pn=10^(-6);
Pmax=0.32;  % 0.32w
delta=3*sqrt(5)/5;
C=min(A,B-A)/delta;

% original Objective
M = length(pm);
pl = pm;
for i = 1:M
    if i+1 < M
        pl(i)=pm(i+1);  
    end
end
fit = -sum(log2(log2(1+(hm.*pm)./(pn+hm.*sum(sum(pl))))));
% penalty 

 Pf = cumprod((1+sum(sum(pm(:)))*((Iter+2*Max_iter)/Max_iter))/Pmax) * ...
       cumprod((1+sum(sum(sqrt(pm(:))))*((Iter+2*Max_iter)/Max_iter))/C );

% Objective
 PHI = fit + Pf;

end

end
   