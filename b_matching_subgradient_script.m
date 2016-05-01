% Author: S M Ferdous
%sferdou@purdue.edu
% 4/1/2016
% b_matching using subgradient methods
% Project of CS520

clc; clear all
%---------------Parameter---------------------
%type=1 if mat format, 2 if dat format 3 if csv format
filetype = 1;
%graphtype =1 if non bipartite 2 otherwise
graphtype = 1;
%0.5 approx or naive algo for bound.1: 0.5 aprox, 2:naive
ubtype = 2;
%do u want to find feasible solution, 1:yes 2:no
feas_need = 1;
%save iteration or not, 1: save 2: dont
save_itn = 2;
%data directory
data_source = 'Data/Input/';
files = {'astro-ph','turon_m','Reuters911','cond-mat-2005','gas_sensor'};
%list of the name of the files
%files = {'toy_graph3'};


% iteration limit to decrease mu
L_it = 40;
% iteration limit to stop
L_max = 120;
% b values
%B_set=[1,3,5,10];
B_set = 1;
%Max iteration Limit
T = 5000;
%-------------------------------------------------------------
for f_ii=1:size(files,2)
	if filetype == 1
		data_url = strcat(data_source,files{f_ii},'.mat');
		load(data_url)
		D = Problem.A;
    elseif filetype == 2
		data_url = strcat(data_source,files{f_ii},'.dat');
		fid = fopen(data_url,'r');
		D = textscan(fid,'%f::%f::%f::%f');
		if graphtype == 2
			d=D{2};
			unik = unique(d);
			for ii=1:size(d,1)
				d(ii) = find(unik==d(ii));
			end
			D{2} =d;
			J = D{1};
			n_A = max(J);
			I = D{2}+n_A;
			S = D{3};
            n_v = max(I);
		end
	else
		data_url = strcat(data_source,files{f_ii},'.csv');
		D = csvread(data_url);
		
		if graphtype == 2
			unik = unique(D(:,2));
			for ii=1:size(D,1)
				D(ii,2) = find(unik==D(ii,2));
			end
			
			J = D(:,1);
			n_A = max(J);
			I = D(:,2)+n_A;
			S = D(:,3);
            n_v = max(I);
		end
	end
	
	% A toy example. with 5 vertices and 5 edges. 
    % the optimal matching using b = 1 for this example is 58. If b>1 then the
    % maximum matching weight is 93 that means to take all the edges.
    %D = [0 29 0 0 2; 29 0 30 0 0; 0 30 0 29 0; 0 0 29 0 3; 2 0 0 3 0];
    
    
	
	if graphtype == 1
        %number of vertices
        n_v = size(D,1);
		[I,J,S] = find(D);
	end
	
    nnzero = length(I);
    E = [J,I,S];
    e = find(J<I);
    %get rid of the duplicate edges.
    E = E(e,:);
    n_e = size(E,1);
    SI = ones(2*n_e,1);
    SJ = zeros(2*n_e,1);
    SW = zeros(2*n_e,1);
    clear  Problem D
    %Find the adjancecy list
    Adj = cell(n_v,1);
    Max_E = zeros(n_v,1);
    f=0;
    for jj=1:n_v
          %find the adjacent edges of a vertex. 
          Adj{jj} = int32([find(E(:,1)==jj)' find(E(:,2)==jj)']);
          if(isempty(Adj{jj})==0) Max_E(jj)= max(E(Adj{jj},3)); end
          l = size(Adj{jj},2);
          SI((f+1):(f+l)) = SI((f+1):(f+l))*jj;
          SJ((f+1):(f+l)) = Adj{jj};
          SW((f+1):(f+l)) = ones(l,1);
          f = f+l;
    end
    %clear Adj
    %sparse incident matrix 
    A = sparse(SI,SJ,SW);


	
    [~,ind_sort] = sort(E(:,3),'descend');
    if ubtype==1
        E_sorted = E(ind_sort,:);
    else
       E_sorted = E;
    end
    degree = zeros(n_v,1);
    %degree of each vertex
    for ii=1:n_v
        degree(ii) = length([find(E(:,1)==ii)' find(E(:,2)==ii)']); 
    end
    
    %storage for each b value the final objective and time 
    Result = zeros(size(B_set,2),6);

	%stroring each iteration output and time
    if save_itn == 1
        G = zeros(T,2*size(B_set,2));
    end
    for b_ii=1:size(B_set,2)
        tic
        %calculate b vector
        b= min(ones(n_v,1)*B_set(b_ii),degree);

        weight = 0;
        btemp = b;
        %0.5 approax algo an upper bound on the dual problem
        %pick the edges satisfying b constraint for each vertex.
        for ii=1:n_e
           u_e = E_sorted(ii,1);
           v_e = E_sorted(ii,2);
           %pick the heaviest one maintaining the constraint
           if btemp(u_e)>0 && btemp(v_e)>0
                weight = weight + E_sorted(ii,3);
                btemp(u_e) = btemp(u_e) - 1;
                btemp(v_e) = btemp(v_e) - 1;
           end

        end

    
        UB = -weight;
        %timer starts
        tic
        %initializing dual variable
        u = zeros(n_v,1);
      
        %initial search direction
        s = zeros(n_v,1);
        max_g = -inf;
        
        
        d = zeros(n_v,1);
        alpha =0;
        mu = 2;
        for ii=1:T

            %calculate x, which minimizes dual function
            x = -E(:,3)+ u(E(:,1)) + u(E(:,2));
            x = (x<=0);

            %calculate objective function of dual
            g = sum(min(-E(:,3)+ u(E(:,1)) + u(E(:,2)),0))- sum(b.*u);
            %bookkeeping
            if save_itn == 1
                G(ii,2*b_ii-1) = -g;
                G(ii,2*b_ii) = toc;
            end

            if g>max_g
                max_g = g;
                it_max = ii;
                x_uopt = x;
                uopt = u;
            end
            %stopping rule
            if(ii-it_max>=L_max)
                break;
            end
            %calculate subgradient vector (slower)
    %         for jj=1:n_v
    %             %find the adjacent edges of a vertex. 
    %             %adj_ind = [find(E(:,1)==jj)' find(E(:,2)==jj)'];
    %             %sub gradient at jj
    %             s(jj) = sum(x(Adj{jj}))-b(jj);
    %         end
    
			%calculate subgradient vector (faster)
            s = A*x-b;
            %updating step length
            if (ii-it_max>=L_it)
                mu = mu/2;
            end
        %     if mu<0.0005
        %         mu=2;
        %     end

            %Traditional update rule
            s = (1-alpha)*s+alpha*d;
            l = mu*(UB-g)/norm(s)^2;
            d = s;
            u = min(u,Max_E);
            u = max(u+l*s,zeros(n_v,1));  
            
            %MFC update rule
            % beta = max(0,-gamma*dot(d,s)/norm(d)^2);
            % s = s+beta*d;
            % l = mu*(UB-g)/norm(s)^2;
            % d = s;
            % u = min(u,Max_E);
            % u = max(u+l*s,zeros(n_v,1));

        end
        %finding feasible solution
        %Ax-b<0
        subgrad_time = toc;
        %if you need a feasible solution. we can run a heuristic for u.
        if feas_need == 1
            tic
            x1=x_uopt;
            e=A*x1-b;
            ind = find(e<0);
            for kk=1:length(ind)

                %[~,s_ind] = sort(E(Adj{ind(kk)},3));
                u = ind(kk);
                %checking still i am under used or not. It might be the case
                %that I have violated the other constraint in the process of
                %finding matched edege
                if e(u)<0
                    adj_ind = Adj{u};
                    [~,s_ind] = sort(E(adj_ind,3),'descend');
                    s_ind = adj_ind(s_ind);
                    cnt = 0;
                    for kk1=1:length(s_ind)
                        e1 = s_ind(kk1);
                        if(E(e1,1)==u)
                            v = E(e1,2);
                        else
                            v = E(e1,1);
                        end
                        if(x1(e1)==0 && E(e1,3)>=0 && e(v)<0)
                            x1(e1) = 1;
                            cnt = cnt+1;
                        end
                        if cnt == abs(e(u))
                            break;
                        end
                    end
                end
                e=A*x1-b;
            end
            %Ax-b>0

            e=A*x1-b;
            ind = find(e>0);
            for kk=1:length(ind)
                %it might be the case that fixing the previous constraint fixed
                %this constraint too. so checking whether this constraint is
                %still violated
                if(e(ind(kk))>0)
                    %[~,s_ind] = sort(E(Adj{ind(kk)},3));
                    adj_ind = Adj{ind(kk)};
                    [~,s_ind] = sort(E(adj_ind,3));
                    s_ind = adj_ind(s_ind);
                    cnt = 0;
                    for kk1=1:length(s_ind)
                        if(x1(s_ind(kk1))==1)
                            x1(s_ind(kk1)) = 0;
                            cnt = cnt+1;
                        end
                        if cnt == e(ind(kk))
                            break;
                        end
                    end
                end
                e=A*x1-b;
            end
            feas_time = toc;
        end
        
        
        if feas_need == 1
            Result(b_ii,:)=[B_set(b_ii), -max_g, it_max, subgrad_time,dot(E(:,3),x1),feas_time];
        else
            Result(b_ii,:)=[B_set(b_ii), -max_g, it_max, subgrad_time,dot(E(:,3),x1),0];
        end
    end
    col_header = {'B','obj','iteration','time','feasible'};
    if(ubtype==1)
        output_file = strcat(files{f_ii},'_result_subgradient','_ubhalf');
    else
        output_file = strcat(files{f_ii},'_result_subgradient','_ubnaive');
    end
    xlswrite(output_file,Result);
    %xlswrite(output_file,col_header,'sheet1','B1');
    if save_itn == 1
        if(ubtype==1)
            output_file = strcat(files{f_ii},'_iteration_subgradient','_ubhalf');
        else
            output_file = strcat(files{f_ii},'_iteration_subgradient','_ubnaive');
        end
        xlswrite(output_file,G);
    end
    %col_header = {'obj','time','obj','time','obj','time','obj','time'};
    %xlswrite(output_file,col_header,'sheet1','B1');
    Result
    %clearvars -except data_source files f_ii L_it L_max B_set filetype graphtype T ubtype
end










