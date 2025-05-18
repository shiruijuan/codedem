clc;clear;close all;t0=clock;

global popsize setgen F CR p r index Q gen
global D T ul dl up down Ereal Qload  Psto Prel 
global  Qmax  P_cl lenge Girdmax  Basto Barel Gird
global  Qe wtm stm Ba Bamax GB nac nec Bammax Qmmax
popsize=300; 
setgen=200;  
F=0.2;
CR=0.5;
p=0.9;  
r=rand;
num_states = 1;
num_actions = 5 ;
Q_table = rand(num_states, num_actions);
terminal_state = 0;
alpha = 0.5;
gamma = 0.9;
epsilon = 0.1;
s=1;
D=8;
dim=D;
T=24;
ul=[ 4000 800 300 4000 800 500 230 1000]  ;
dl=[ 2000 100 50 2000 100  50 100 500]; 
up = [  800 800 200 500 800 350 230 1000  ];
down=[  800 800 200 500 800 350 230 1000];
wtm=[1902 1842 1580 1216 848 718 634 966 604 486 1546 1602 1274 830 1512  1702 1764 2104 2254 2360 2342 2290 2170 1984 ]';
stm=[0.00 0.00 0.00	0.00 0.0000 158.4 256 566.4 676.8 1064 1280 1188.2 1156.8 1536 1272 948.8 517.4 320 150.4  0.00	0.00	0.00	0.00	0.00];
Ereal=[3003.2	2573.6	2388.8	2246.4	2406.4	2428.8	3693.6	2868	3036	3543.2	3521.04	4445.16	4524.04	4660.04	4445.16	3963.04	3818.88	3909.32	3684.92	3548.92	3589.6	3545.6	3612.8	3248];
P_cl=[304	354	374	469	473	 502 501 497 570 623	684	756	766	727	688	531	508	502	499 471	504	450 367	401];
Qload=[3083 3083 3052 3083 3098 3120 3298 3482 3236 3413 3520 3744 3713 3575 3635 3620 3551 3329 3398 3320 3498 3604 3389 3083];
Girdmax=5000;  
Gird=zeros(1,T);  
Qmax=1200;   
Bamax=1200;  
Psto=zeros(1,T);
Prel=zeros(1,T);
Basto=zeros(1,T);
Barel=zeros(1,T);
GB=zeros(1,T);
lenge=zeros(1,T);
Q=zeros(1,T+1);
Qe=zeros(1,T+1);
Q(1)=350;
Qmmax=2000;
Ba(1)=400;
Bammax=2000;
nec=0.6;
nac=0.93;

for i=1:popsize
    pop(i)=chushihua(ul,dl);
end
 pop = non_domination_sort_mod(pop);
 denngji_values1 = [pop.denngji];
 AC = pop(denngji_values1 == 1);
 TOP=TOPSIS(AC);
 arraySize = find(TOP == 1);
 BbestX = AC(arraySize);       
 bestX =BbestX.x;
for gen=1:setgen
    Wa=0.9;
    Wb=0.1;
    epsilon=Wa-(Wa-Wb)*(gen/setgen);
    state = randi([1, num_states]); 
    action=randi([1, 5]);
if gen >1
    if rand() < epsilon
         action = randi([1, num_actions]); 
    else
         [~, action] = max(Q_table(state, :)); 
    end
   
end
 
 for i=1:popsize
        n=0.05*exp(-2*(gen/setgen)^2);
        if p<r
            v1(i).x=pop(i).x+n.*(1+sin(r))*pop(i).x;
        else
            v1(i).x=pop(i).x.*(n*(2*rand(1,D)-1)+1);
        end
        index=action;
        m=2*sin(r+pi/2);
        s = randi([1,popsize],1);
        worseness=pop(s);
        ori_value = rand(1,D);
        cauchy_value = tan((ori_value-0.5)*pi);
        AB = [ BbestX ,worseness ];
        Ab=Copy_of_sel_pareto(AB,worseness);
        if Ab==1
            v1(i).x=pop(i).x+cauchy_value(:,D).*( pop(i).x-bestX);
            vll= v1(i).x;
            v1(i).x=mutations(vll,pop,bestX,dl,ul,D,gen,setgen,index,i);
        else
            v1(i).x=pop(i).x+cauchy_value(:,D).* (bestX-m.*pop(i).x);
            vll= v1(i).x;
            v1(i).x=mutations(vll,pop,bestX,dl,ul,D,gen,setgen,index,i);
        end
    u(i)=chushihua(ul,dl);
    u(i).x=v1(i).x;
    u(i)=yueshu(u(i));
     [u(i).obj_1,u(i).obj_2,u(i).obj_3]=evaluation(u(i));
     PopObj(i,:) = [u(i).obj_1, u(i).obj_2,u(i).obj_3];  
 end
Result(gen,1)=MOVI(PopObj,u); 
 if gen>1
     if Result(gen)<Result(gen-1)  
        r = 2;
        new_state =randi([1, 1]);
    else
        r = -1;
        new_state = randi([1, 1]);
     end   
     Q_table(state, action) = Q_table(state, action) + alpha*(r + gamma* max(Q_table(new_state, :)) - Q_table(state, action));
     state = new_state;
end
    pop=u;
    AC=[AC,u];
    AC = Copy_of_non_domination_sort_mod(AC);
    limit = 50;
    denngji_values = [AC.denngji];
    AC_denngji_1 = AC(denngji_values == 1);
    if length(AC_denngji_1) <= limit
        AC = AC_denngji_1;
    else
        z_fifth_column_values = [AC_denngji_1.z_fifth_column];
        [~, sorted_indices] = sort(z_fifth_column_values, 'descend');
        AC = AC_denngji_1(sorted_indices(1:limit));
    end
    BbestX = AC(1);          
    bestX = BbestX.x;    
end

function state = getState(div, con, gen)     
    if div(gen) > div(gen-1) & con(gen) > con(gen-1)
        state = 1;
    elseif div(gen) > div(gen-1) & con(gen) <= con(gen-1)
        state = 2;
    elseif div(gen) <= div(gen-1) & con(gen) > con(gen-1)
        state = 3;
    else
        state = 4;
    end
end
