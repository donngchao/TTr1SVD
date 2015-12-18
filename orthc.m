function O=orthc(A,varargin)
% O=orthc(A) or O=orthc(A,tol)
% ----------------------------
% Computes the outer product vectors that span the orthogonal complement
% tensors.The result that has been calculated is to iterate along the tree
% data structure.
%
% O         =   cell, each column of O contains outer product vectors that form
%               a tensor orthogonal to A,
%
% A         =   array, d-way array,
%
% tol       =   scalar, optional tolerance to decide whether computed singular
%               values are numerically zero. Default value = max(size(A))*eps(sigmas(1)).
%
% Reference
% ---------
%
% A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms
% http://arxiv.org/abs/1407.1593
%
% 2015, Kim Batselier, Haotian Liu, Ngai Wong

n=size(A);
r=zeros(1,length(n)-1);%The tree has a deep degree of length(n)-1.
for i=1:length(n)-1
    r(i) = min(n(i),prod(n(i+1:end)));%For every level,there should be some nodes.
end
totalsvd=1;
svdsperlevel=zeros(1,length(r));   % add length 1 for the first level
svdsperlevel(1)=1;  % The first level should just contains only 1 svd node,we view it as a node here.
for i=2:length(r)   %Here explains other level contain how many nodes are contained in the structure.
    svdsperlevel(i)=prod(r(1:i-1));%I see r plays a very important role here in deciding how many svds on this level.But why?
    totalsvd=totalsvd+svdsperlevel(i);%Total number of svds in the whole levels of the tree. 
end
nleaf=prod(r);%This will tell us the tree has how many leaves.After all it is still a tree.

U=cell(1,totalsvd);
S=cell(1,totalsvd);
V=cell(1,totalsvd);
%First reshape the tensor to become a new matrix,and then using the svd operation on this matrix.
%Precisely,it's better to understand this process as decrease the dimension of the tensor to only 2,
%and then calculate the svd as usual. 
[Ut St Vt]=svd(reshape(A,[n(1),prod(n(2:end))]));
U{1}=Ut;%U{1} stores the first U matrix after calculating the svd for the first time after the reshaping of the tensor.
S{1}=diag(St);%S{1} stores the eigen values of the first svd result.
V{1}=Vt;%V{1} stores the firsr V matrix after calculating the svd for the first time after the reshaping of the tensor.
counter=2; % this counter keeps track during the iterations which V{i} we are computing. This is a linear counter that counts breadth-first
whichvcounter=1;    % this counter keeps track during the iterations of which V we are taking the svd
%In fact this is how the algorithm will retrieve the whole tree.
for i=1:length(r)-1           % outer loop over the levels
    for j=1:prod(r(1:i))      % inner loop over the number of svds for this level 
        if rem(j,r(i)) == 0
            col=r(i);
        else
            col=rem(j,r(i));
        end
        [Ut St Vt]=svd(reshape(V{whichvcounter}(:,col),[n(i+1),prod(n(i+2:end))]));
        U{counter}=Ut;
        S{counter}=diag(St);
        V{counter}=Vt;
        counter=counter+1;
        if rem(j,length(S{whichvcounter}))==0
%             V{whichvcounter}=[];
            whichvcounter =  whichvcounter+1;
        end
    end
%     whichvcounter = whichvcounter+1;
end

Slevel=cell(1,length(r));   % cat each level singular values into 1 vector
counter=1;
for i=1:length(r),
    for j=1:svdsperlevel(i),
        Slevel{i}=[Slevel{i}; S{counter}];
        counter=counter+1;
    end
end

for i=1:length(r),             % make all singular value vectors the same size (number of leaves)
    Slevel{i}=kron(Slevel{i}, ones(nleaf/length(Slevel{i}),1));
end

sigmas=ones(nleaf,1);         % output singular values at each leaf
for i=1:length(r),
    sigmas=sigmas.*Slevel{i};
end
%This MATLAB function returns the number of input arguments passed in the call to the currently executing function.
if nargin==2                  %returns the number of function input arguments.
    tol=varargin{1};
else
    tol= max(n)*eps(sigmas(1));
end

counter=1;
for i=1:svdsperlevel(end)
    % for each svd of the last level
    
    % how many leaves per svd?
    leavespersvd=prod(r)/svdsperlevel(end);
    indices=leave2ind((i-1)*leavespersvd+1,n);    
    numberofUs=size(U{indices(1)},2);
    
    % check for numerically zero sigmas
    I=find(sigmas((i-1)*leavespersvd+1:i*leavespersvd)<tol);
    if ~isempty(I)
        for k=1:length(I)
            for j=1:length(r)-1
                O{j,counter}=U{indices( (length(r)-j)*2+1 )}(:,indices((length(r)-j)*2+2));
            end
            O{length(r),counter}=U{indices(1)}(:,I(k));
            O{length(r)+1,counter}=V{indices(1)}(:,I(k));
            counter =counter+1;
        end        
    end
    
    for k=1:numberofUs
        Vindices=1:n(end);
        if k <=n(end)
            Vindices(k)=[];
        end
        for l=1:length(Vindices)        % first we go down the tree up until the last level            
            for j=1:length(r)-1
                O{j,counter}=U{indices( (length(r)-j)*2+1 )}(:,indices((length(r)-j)*2+2));
            end
            O{length(r),counter}=U{indices(1)}(:,k);
            O{length(r)+1,counter}=V{indices(1)}(:,Vindices(l));
            counter =counter+1;
        end
    end    
end

end
