%% small demo to illustrate the use of ttr1svd.m and getAtilde.m

clear all;
clc;

% example tensor from "A Constructive Algorithm for Decomposing a Tensor into a Finite Sum of Orthonormal Rank-1 Terms"
A(:,:,1)=[1 4 7 10;2 5 8 11;3 6 9 12];
A(:,:,2)=[ 13 16 19 22; 14 17 20 23; 15 18 21 24];

% store dimensions of A in n
n=size(A);

% compute the TTr1 decomposition of A
[U,S,V,sigmas]=ttr1svd(A);

% inspect sigmas
sigmas

% convert the TTr1 format to the Tucker (HOSVD) format
[S,Q]=ttr12tucker(U,sigmas,V,size(A));

% check orthogonality of mode vectors
norm(Q{1}'* Q{1}-eye(3),'fro')
norm(Q{2}'* Q{2}-eye(4),'fro')
norm(Q{3}'* Q{3}-eye(2),'fro')

% sparsity of the tensor core S in percent
100*numel(find(S))/numel(S)

% try a rank-2 approximation
Atilde=getAtilde(U,sigmas(1:2),V,1:2,n);

% compare our approximation error with bound from singular values
[norm(reshape(A-Atilde,[1 prod(n)])) norm(sigmas(3:end))]

% convert the rank-2 approximation to the Tucker (HOSVD) format
% we can keep the original U,V but need to make a new 'sigmas' vector that
% contains only the nonzero singular values we are interested in
sigmas2=zeros(length(sigmas),1);sigmas2(1:2)=sigmas(1:2);
[S2,Q2]=ttr12tucker(U,sigmas2,V,size(A));


% try a different rank-2 approximation
Atilde2=getAtilde(U,sigmas([2 4]),V,[2 4],n);

% compare our approximation error with bound from singular values
[norm(reshape(A-Atilde2,[1 prod(n)])) norm(sigmas([1 3 5 6]))]

% determine orthonormal nullspace terms: needs a sigma of value 1
nullspace1=getAtilde(U,1,V,5,n);
nullspace2=getAtilde(U,1,V,6,n);

% check whether nullspace tensors are orthogonal to original A
A(:)'*nullspace1(:)
A(:)'*nullspace2(:)

% demonstrate use of leave2ind.m to determine which U,V vectors we need to reconstruct leaves
indices=leave2ind([1 3 4],n)

% odd columns of indices contain cell indices i of U{i},V{i},
% even columns of indices contain column indices j of U{i}(:,j),V{i}(:,j)

% reconstruct the first, 3rd and 4th leaf
firstTerm=sigmas(1)*mkron(V{indices(1,1)}(:,indices(1,2)),U{indices(1,1)}(:,indices(1,2)),U{indices(1,3)}(:,indices(1,4)));
thirdTerm=sigmas(3)*mkron(V{indices(2,1)}(:,indices(2,2)),U{indices(2,1)}(:,indices(2,2)),U{indices(2,3)}(:,indices(2,4)));
fourthTerm=sigmas(4)*mkron(V{indices(3,1)}(:,indices(3,2)),U{indices(3,1)}(:,indices(3,2)),U{indices(3,3)}(:,indices(3,4)));

% compare to getAtilde
Atilde3=getAtilde(U,sigmas([1 3 4]),V,[1 3 4],n);
norm(firstTerm+thirdTerm+fourthTerm-Atilde3(:))
