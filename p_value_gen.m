% p_value_gen: a function taking as inputs 
% (a) the three time series X,Y,Z (row vectors) that each of their elements take values (0, 1, ..., X_alph_size) 
%       (0, 1, ..., Y_alph_size) and (0, 1, ..., Z_alph_size) respectively
% (b) the alphabet sizes of series X, Y, Z (X_alph_size, Y_alph_size, Z_alph_size)
% (c) the Markov depth D

% Note: If we wish to check the causal influence of X to Y without
% conditioning on Z, set Z as an all-zero vector of the same length as X and Y, and set Z_alph_size=1







function p_val= p_value_gen(X_alph_size, Y_alph_size, Z_alph_size, X, Y, Z, D )

%Initialise context (with a large amount of zero rows)
a_context=zeros(length(X),3*D+3);
count=zeros(length(X),1);
% Initialise complementary variables
row_num=1;
cell_non_zeros={[]};
r=0;

%In the end of the "for" loop below, want: 
%(a) "a_context" to list all X,Y,Z contexts (the first row will be the all zero context)
% (Note that by context we mean each vector [X((i-D):i), Y((i-D):i), Z((i-D):i)]   )
%(b) Each element in "count" counts the number of times each context appears in the respective row of a_context (the first element is the number of times the all-zero context appears)
%(c) "row_num" will be the number of different contexts in our data
%(d) Each element of "cell_non_zeros" lists the indices of non-zeros in the context in the respective row of "a_context"



for i=(D+1):length(X)
    
    % "new_row" is the current context
    new_row= [X((i-D):i), Y((i-D):i), Z((i-D):i)];
    % "non_zero" finds indices of non-zero elements in new_row
    non_zero=find(new_row);
    
    
    % Increase the count of the all-zero vector by 1
    % if there are no ones in the vector above (useful step computationally when data is sparse)
    if isempty(non_zero)==1
        count(1)= count(1)+1;
    
    else
        
        % This for loop identifies if the current context in "new_row" has already appeared in the data
        % If the current context has already appeared, the auxiliary variable r takes the
        % value of the index of this context
        for j=1:row_num
            if length(cell_non_zeros{j})==length(non_zero)              
                if range(a_context(j,:)-new_row)==0                                
                    r=j;
                end
            end
        end
        
        % If the current context has already appeared, increase its count by 1 
        if r~=0
            count(r)=count(r)+1;
            r=0;
          
            
        % If the current context has not yet appeared, 
        % add a new entry to variables "cell_non_zeros", "a_context" and "count"
        % and increase "row_num" by 1    
        else          
            cell_non_zeros{row_num+1}=non_zero;
            a_context(row_num+1, :)= new_row;            
            count(row_num+1)=1;
            row_num=row_num+1;            
        end    
    end                
end



% Remove any unused rows in "a_context" and the unused elements in "count"
count((row_num+1):end)=[];
a_context((row_num+1):end,:)=[];
if count(1)==0
    count(1)=[];
    a_context(1,:)=[];
end



















%Initialise "a_cxt_yz_draft" as the Y and Z elements of "a_context"  
a_cxt_yz_draft=a_context(:, (D+2):end);
%Initialise "a_context_yz" and "count_yz"
a_context_yz=zeros(1,2*D+2);
count_yz=[0];
% In the end, we want "a_context_yz" to list all Y,Z contexts
% And each element in "cout_yz" to be the number of times the respective context appears


% Now, some rows in a_cxt_yz draft are the same, 
% need to merge them in "a_context_yz" using the "for" loop below
for i=1:size(a_cxt_yz_draft,1)
    new_row_yz= a_cxt_yz_draft(i,:);
    if new_row_yz==zeros(1,2*D+2)
        count_yz(1)= count_yz(1)+count(i);
    else
    check=ismember(a_context_yz, new_row_yz, 'rows');   
    idx=find(check==1);
        if isempty(idx)==1
            a_context_yz= [a_context_yz ; new_row_yz];
            count_yz=[count_yz; count(i)];
        else
            count_yz(idx(1))=count_yz(idx(1))+count(i);
        end
    end
end

% Remove the all zero context vector if it does not appear
if count_yz(1)==0
    count_yz(1)=[];
    a_context_yz(1,:)=[];
end










% Initialise I1, I2, which will be used to calculate the Statistic of the hypothesis test
I1=zeros(length(a_context),1);
I2=zeros(length(a_context_yz),1);

%Calculate I1
for i=1:size(a_context,1)
    %Initialise sum_I1
    sum_I1=count(i);
    
    %Calculate sum_I1 using the "for" loop below
    for j=setdiff([1:Y_alph_size], (a_context(i, 2*D+2)+1) ) 
        
        new_row = [a_context(i,1:(2*D+1)),    j-1, a_context(i,(2*D+3):(3*D+3))  ];
        log_1 =ismember(a_context, new_row, 'rows');
        idx=find(log_1==1);
            if isempty(idx)==0
                sum_I1=sum_I1+count(idx);
            end 
    end
    
    %Calculate I1(i)
    I1(i)= count(i) * log(count(i)/sum_I1);
end



%Calculate I2
for i=1:size(a_context_yz,1)
    %Initialise sum_I2
    sum_I2=count_yz(i);
    %%Calculate sum_I2 using the "for" loop below
    for j=setdiff([1:Y_alph_size], (a_context_yz(i, D+1)+1) ) 
        
        new_row_2 = [a_context_yz(i,1:D),    j-1, a_context_yz(i,(D+2):(2*D+2))  ];
        log_2 =ismember(a_context_yz, new_row_2, 'rows');
        idx_2=find(log_2==1);
            if isempty(idx_2)==0
                sum_I2=sum_I2+count_yz(idx_2);
            end    
    end
    
    %Calculate I2(i)
    I2(i) = count_yz(i) * log(count_yz(i)/sum_I2);
end


% Calculate the statistic
Statistic=2*(sum(I1)-sum(I2));

% Find the degrees of freedom of the hypothesis test
dof=  Y_alph_size^D * Z_alph_size^(D+1) *(X_alph_size^(D+1) -1)*(Y_alph_size-1);

% Calculate the p-value of the hypothesis test
p_val=1-chi2cdf(Statistic,dof);



end