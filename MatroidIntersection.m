function I=MatroidIntersection(E,S1,S2,option) % find the matroid intersection between matroids fromed by columns of S1 and S2; E is the base set;
%   optin=[1,2]% S1 is a numerical matrix, and S2 is a structured matrix
%   option=[1,1]% S1 and S2 are both numerical matrices
%   I is the intersection of the two matriods 
%   The algorithm is due to the document: Graph Algorithms and
%   Combinatorial Optimization Dr. R. Chandrasekaran  Presenters: Benjamin
%   Ferrell and K. Alex Mills%  http://www.utdallas.edu/~kam093020/papers/matroid-intersection.pdf
%%%%%%%%%
%%%%%%%%% We just consider option=[1,2] by default; the case option=[1,1]
%%%%%%%%% can be easily obtained. 
%  Example:  E=[1 2 3 4 5];S1=[1 0 1 0 1;1 0 1 1 0]; S2=[1 0 1 0 0;0 1 0 1
%  0;0 0 0 0 1]; option=[1,2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin==3 % option=[1,2]
I=[];E_I=E; % Step 1: Initilize   I: left set; E_I: right set E_I=E-I;
MinRow=min(size(S1,1),size(S2,1));
n=numel(E); % cardinality of the base set
%%%% Step 2: 


while 1

for z=E_I;   % step 2
   Iaddz=[I,z];
   leftrank=rank(S1(:,Iaddz));
   rightrank=sprank(S2(:,Iaddz)); % if option=[1,1], change sprank to be rank;
    if leftrank==numel(Iaddz) && rightrank==numel(Iaddz) % add z to I
       I=[I,z]; 
    end   
    if leftrank==MinRow
           break
    end
end
E_I=setdiff(E,I);% update E_I;

%  step 3  construct the Krogdahl graph
D=zeros(n,n);  % D_{i,j}=1 means a edge from i to j
for i=I
    for j=E_I
   I_i=setdiff(I,i); I_iAddj=[I_i,j];
   if rank(S1(:,I_iAddj))==numel(I_iAddj);
       D(i,j)=1;
   end
    if sprank(S2(:,I_iAddj))==numel(I_iAddj)
       D(j,i)=1;   
    end
    end
end
% step 3 continued construct x1 and x2
X1=[];X2=[];

if numel(I)<MinRow  % only if this case, X1 or X2 will not be empty
for z=E_I
    I_z=[I,z];
    if rank(S1(:,I_z))==numel(I_z) 
        X1=[X1,z];
    end
    if sprank(S2(:,I_z))==numel(I_z) 
        X2=[X2,z];
    end   
end
end

% step 4  find the short path from X1 to X2
  if numel(X1)>0 && numel(X2)>0 % either wise stop 
    KroGraph=zeros(n+2,n+2); % start node: n+1;  end node: n+2
    KroGraph(1:n,1:n)=D;
    KroGraph(n+1,X1)=1;
    KroGraph(X2,n+2)=1;
    
    SparseKroGraph=sparse(KroGraph);
    [dist,path,pred]=graphshortestpath(SparseKroGraph,n+1,n+2); % find the shortest path from X1 to X2
    
    if dist>2 && dist~=Inf % the shorest path exists
        Qj=[];Qi=[];
        for vecsort=1:dist/2-1
            Qj=[Qj,path(vecsort*2)];
            Qi=[Qi,path(vecsort*2+1)];
        end
        Qj=[Qj,path(dist)];
        I=setdiff(I,Qi);
        I=union(I,Qj);  %augment I using Q
        I_E=setdiff(E,I);      
    else      
        break
     %   stop=1;
    end
  else
      break    
  end
  
  if numel(I)==MinRow
      break
  end
  
% step 4 find the 
end  % end the procedure 
else 
   % I=null;
    I=[];
    disp('please do not input the "option"');
    disp('If S1 and S2 are both numerical matrices, you just need to replace all the "sprank" by "rank"');
end