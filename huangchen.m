%   uunifast single core tasks utilization: 4 tasks; utilization bound 0.55
%   columns are tasks and rows are modes: 4 tasks and 3 modes
u = [0.027717584367889847   0.10954825924790529    0.0014241633805865128        0.41130999300361837;
     0.3277699959119918     0.013513632001002363   0.13101584660539967          0.07770052548160622;
     0.14626938552145202    0.14165934558077242    0.192184073329005            0.0698871955687706];

%   number of tasks
n_tasks = size(u,2);
%   number of modes
n_modes = size(u,1);

%   arbitrary periods (lower prio ------> higher prio)
t = [100,     85,    70,   60
     100,     70,    55,   50
     80,      55,    40,   35];
 
%   implicit deadlines -> deadline = period
d = t;
        
%   calculating the WCET, considering that U = WCET/PERIOD        
c = u.*t;

%   starts the analysis from the lowest priority task
task_k = 1;
%   fixed mode for the FTP, as all modes of a task have the same priority
mode = 3;

%   starts the analysis assuming the system is schedulable
is_sched = true;

while ( (task_k < n_tasks) && is_sched )
    %   allocates an array to store the Ci_max of each task
    betha = zeros(n_tasks - task_k);
    
    %   calculates the betha only for the higher priority tasks than tk
    for i = 1 : n_tasks - task_k
       betha(i) = max( c(:, i + task_k) )/max( u(:, i + task_k) );
    end

    %   allocates an array to store the betha index ordering
    betha_index = zeros(n_tasks - task_k);
    for i = 1 : n_tasks - task_k
        %   finds the position of the biggest value among the bethas
        betha_index(i) =  find(betha==max(betha),1,'first');
        %   attributes a -inf value to the previous biggest value so it is never considered as max again
        betha(betha_index(i)) = -Inf;
    end

    %   theorem 1 given an arbitrary deadline of the task and the mode selected above
    sum_Ci = 0;
    for k = task_k + 1 : n_tasks
       sum_Ci = sum_Ci + max(c(:,k)); 
    end

    %   returns true or false for the inequality
    fprintf('TheoremPt1 i = %i, Ci =  %i \n', task_k, sum_Ci);
    theorem1_pt1 = d(mode, task_k) - sum_Ci - c(mode, task_k) >= 0 ;

    %   calculates the second part of theorem 1 (still needs to be correctly indexed according to betha)
    sum_Ui = 0;
    sum_Cj = 0;
    sum_Ci = 0;  
    for i = 1 : n_tasks - task_k
        for j = i : n_tasks - task_k
           sum_Cj =  sum_Cj + max( c( :, betha_index(j) + task_k ) );
        end
        sum_Ui = sum_Ui + max( u( :, betha_index(i) + task_k ) ) * (d(mode, task_k) - sum_Cj);
        sum_Ci = sum_Ci + max( c( :, betha_index(i) + task_k ) );
    end
    %   returns true or false for the inequality
    fprintf('TheoremPt2 i = %i Ui = %i,  Ci = %i \n', task_k, sum_Ui, sum_Ci);
    theorem1_pt2 = c(mode, task_k) <= d(mode, task_k) - sum_Ui - sum_Ci;
    fprintf('Theorem1 %d %d\n', theorem1_pt1, theorem1_pt2)
    is_sched = theorem1_pt1 && theorem1_pt2;
    %   increments the task counter
    task_k = task_k + 1;
end
