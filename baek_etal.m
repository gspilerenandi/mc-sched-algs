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
p = [60      125     200    125;
     50      100     150    80;
     35      75      50     30];
%   deadlines are equal to periods
d = p;
%   calculating the WCET, considering that U = WCET/PERIOD        
e = uti.*p;
%   calculate the number of tasks and modes (better readability of the code)
n_tasks = size(uti,2);
n_modes = size(uti,1);

%   initial upper-bound slack values of each task 
S = Inf(1,n_tasks);
%   initial slack values of each task in each e
s = zeros(size(uti));   
%   the first upper bounds are
rkx = zeros(size(uti));
rkx(1,:) = (d(1,:) - s(1,:));

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
    for task = 1 : n_tasks
        keepgoing = true;
        while (is_sched && keepgoing)
            %   stores the slack values of the current and the upcoming mode before updating the slacks
            old_x_y_s = s(x : y,:);
            %    line 3
            rix = d(x, task) - s(x, task);
            riy = d(y, task) - s(y, task);
            %   mode x -> (g)
            for k = 1 : n_tasks
                rkx(x,k) = theorem_1(x, k, d, p, s, e, m, true);
                if rkx(x,k) < d(x, k)
                   s(x, k) = min(d(x, k) - rkx(x,k), S(task));
                end
            end
            %   mode y -> (h)
            for k = 1 : n_tasks
                rkx(y,k) = theorem_1(y, k, d, p, s, e, m, false);
                if rkx(y,k) < d(y, k)
                   s(y, k) = d(y, k) - rkx(y,k);
                end
            end
            if s(x : y,:) == old_x_y_s
                for k = 1 : n_tasks
                    if ((rkx(x,k) <= d(x,k)) && (rkx(y,k) <= d(y,k)))
                        S(k) = s(y, k);
                        keepgoing = false;
                    else
                        is_sched = false;
                    end 
                end
            end
        end
        if ~is_sched
           break 
        end
    end
    if ~is_sched
       break 
    end
end
%   cheks the amount of time it took to calculate everything
toc;

%   theorem 1
function rku = theorem_1(u, k, d, p, s, e, m, is_it_g)
    n_iteractions = 0;
    rkux = e(u,k);
    %   checks if the given u is either g or h
    if(is_it_g)
        mode = u;
    else
        mode = u - 1;
    end
    %   fixed-point iteration
    while true
        next_rkux = e(u, k) + floor( (1/m)  *  (calcSigma(d, e, p, rkux, mode, k, s)) );
        if(next_rkux == rkux)
            %   if it converged, skip the iteration and return the converged value
            fprintf('número de iteracoes: %i\n', n_iteractions)
            %disp(n_iteractions);
            break;
        else
            %   updates the rkux for the next iteration
            n_iteractions = n_iteractions + 1;
            rkux = next_rkux;
        end
    end
    rku = next_rkux;
end

%   deadlines, wcet, periods, interval (rku(x) here), current mode, task not to consider, previous response time, slack, initial mode transition?
function sigma = calcSigma(d, e, p, l, mode, k, s)
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
            %   is it really e(mode,k) ?
            sigma = sigma + min(Wigh, l - e(mode,k) + 1);
        end
    end
end

%   equation 8    
%   period, wcet, interval (rku(x) here), current mode, task position
function F = calcF(l, p, e, mode, i)
    if l > 0
        floor_div = floor(l/p(mode,i));
        F = floor_div * e(mode,i) + min([e(mode,i), l - floor_div*p(mode,i)]);
    else
       F = 0;
    end
end
