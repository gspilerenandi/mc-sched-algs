%   this script is based on the paper https://ieeexplore.ieee.org/document/9088160

%   tic -> toc: elapsed time to perform all the calculations
tic;
%   number of processors
m = 2;
%   uunifast single core tasks utilization: 4 tasks; utilization bound 0.55
%   columns are tasks and rows are modes: 4 tasks and 3 modes
%   modes are in ascending order of priority, i.e. row 3 has a higher priority than row 1
uti = [0.41130999300361837    0.0014241633805865128   0.10954825924790529     0.027717584367889847;
       0.07770052548160622    0.13101584660539967     0.013513632001002363    0.3277699959119918;
       0.0698871955687706     0.192184073329005       0.14165934558077242     0.14626938552145202];
%   arbitrary periods
p = [60      70     85    100;
     50      55     70    100;
     35      40     55    80];
%   deadlines are equal to periods
d = p;
%   calculating the WCET, considering that U = WCET/PERIOD        
e = uti.*p;
%   calculate the number of tasks and modes (better readability of the code)
n_tasks = size(uti,2);
n_modes = size(uti,1);

%   initial upper-bound slack values of each task 
S = Inf(1,n_tasks);
%   initial slack values of each task in each mode
s = zeros(size(uti));   
%   the first upper bounds are
rkx = zeros(size(uti));
%%rkx(1,:) = (d(1,:) - s(1,:));

%   initializing the verification parameter of the loops
is_sched = true;
%   for every mode transition
for mode = 1 : (n_modes - 1)
    x = mode;
    y = x + 1;
    %   zeroing the current and the upcoming slacks of all tasks between mode and mode + 1 (paper->line 3)
    %   zeroes all the columns of two lines of 's'
    s(x : y, :) = zeros(2,n_tasks);
    %   for every task
    %for task = 1 : n_tasks
    keepgoing = true;
    slack_iter = 0;
    while (is_sched && keepgoing)
        %   stores the slack values of the current and the upcoming mode before updating the slacks
        fprintf('Iteration %i\n', slack_iter)
        old_x_y_s = s(x : y,:);
        %    line 3, unused?
        %   rix = d(x, task) - s(x, task);
        %   riy = d(y, task) - s(y, task);

        %   mode x -> (g)
        for k = 1 : n_tasks
            rkx(x,k) = theorem_1(x, k, d, p, s, e, m, true);
            if rkx(x,k) < d(x, k)
               s(x, k) = min(d(x, k) - rkx(x,k), S(k));
            end
        end
        fprintf('X: %i =>', x)
        disp(s(x,:))
        %   mode y -> (h)
        for k = 1 : n_tasks
            rkx(y,k) = theorem_1(y, k, d, p, s, e, m, false);
            if rkx(y,k) < d(y, k)
               s(y, k) = d(y, k) - rkx(y,k);
            end
        end
        fprintf('Y: %i =>', y)
        disp(s(y,:))
        if isequal(s(x : y,:),old_x_y_s)
            if isequal( le(rkx(x,:), d(x,:)), le(rkx(y,:), d(y,:)) ) 
                S(:) = s(y, :);
                keepgoing = false;
            else
                is_sched = false;
            end 
        end
        %keepgoing = false;
    slack_iter =  slack_iter + 1;
    %fprintf('-------------------------------------\n')
    end
    if ~is_sched
       break 
    end
    %end
    %if ~is_sched
    %   break 
    %end
end
fprintf('-- Slack after: \n')
disp(transpose(s))
%   cheks the amount of time it took to calculate everything
toc;

%   theorem 1
function rku = theorem_1(u, k, d, p, s, e, m, is_it_g)
    rkux = e(u,k);
    %   checks if the given u is either g or h
    if(is_it_g)
        mode = u;
    else
        mode = u - 1;
    end
    %   fixed-point iteration
    while true
        next_rkux = e(u, k) + floor( (1/m)  *  (calcSigma(d, e, p, rkux, mode, k, s, is_it_g)) );
        if(next_rkux == rkux)
            %   if it converged, skip the iteration and return the converged value
            break;
        else
            %   updates the rkux for the next iteration
            rkux = next_rkux;
        end
    end
    %fprintf('r of task %i on mode %i is %f\n', k, u, next_rkux)
    rku = next_rkux;
end

%   deadlines, wcet, periods, interval (rku(x) here), current mode, task not to consider, previous response time, slack, initial mode transition?
function sigma = calcSigma(d, e, p, l, mode, k, s, is_it_g)
    sigma = 0;
    %   from
    g = mode;
    %   to
    h = g + 1;
    %   consider the complete task set of a mode
    for i = 1 : size(e,2)
        %   only considers the intereference of tasks different than task k (i.e, task k cannot interfere itself)
        if (i ~= k)
            Wigh = 0;
            %   verifies if task i has a higher priority than task k (Rate-Monotonic)(NOT SURE IF THIS IS THE CORRECT APPROACH)
            if(p(mode,i) < p(mode,k))
                Wig_l = calcF(l + d(g,i) - s(g,i) - e(g,i), p, e, g, i);
                Wih_l = calcF(l + d(h,i) - s(h,i) - e(h,i), p, e, h, i);
                Wigh1 = 0;
                %   I've used the max function just to make sure that deltag was at least equal to 1
                for deltag = 1 : max(1, floor((l + d(g,i) - s(g,i) - e(g,i))/p(g,i)) )
                    %   equation 9
                    Wigh1 = max( Wigh1, deltag*e(g,i) + calcF(l + (d(g,i) - s(g,i) - e(g,i)) - deltag*p(g,i), p, e, h, i) );
                end
                Wigh2 = 0;
                %   I've used the max function just to make sure that deltah was at least equal to 1
                for deltah = 1 : max(1, floor((l + p(h,i))/p(h,i)))
                    %   equation 10
                    Wigh2 = max( Wigh2, deltah*e(h,i) + calcF( l + p(h,i) - e(h,i) - (p(g,i) - d(g,i) + s(g,i)) - deltah*p(h,i), p, e, g, i) );
                end
                Wigh = max([Wig_l, Wih_l, Wigh1, Wigh2]);
            end
            if(is_it_g)
                sigma = sigma + min(Wigh, l - e(g,k) + 1);
            else
                sigma = sigma + min(Wigh, l - e(h,k) + 1);
            end
        end
    end
end

%   equation 8    
%   interval (rku(x) here), period, wcet, current mode, task position
function F = calcF(l, p, e, mode, i)
    if l > 0
        floor_div = floor(l/p(mode,i));
        F = floor_div * e(mode,i) + min(e(mode,i), l - floor_div*p(mode,i));
    else
       F = 0;
    end
end
