% Submitted with the manuscript:
% A TECHNICAL NOTE ON IMPROVING COMPUTATIONAL EFFICIENCY IN LARGE LINEAR 
% INVERSE PROBLEMS: AN EXAMPLE FROM CARBON DIOXIDE FLUX ESTIMATION
% Vineet Yadav1 and Anna M. Michalak2
% [1]{Department of Global Ecology, Carnegie Institution for Science, 
% Stanford, California, USA 94305}
% THIS IS MATLAB CODE THAT CAN BE RUN EITHER IN MATLAB OR OCTAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose:To Demonstrate Algorithm Proposed in Section 2 of the manuscript
% ALL VARIABLES AND INDICES ARE AS DEFINED IN THE MANUSCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code demonstrates matrix multiplication between an arbitrary matrix
% H and an arbitrary matrix Q that can be expressed as a Kronecker product
% of two matrices.
% The code can also be used to compute H*Q*H' provided that Q is expressed
% as Kronecker product of two square matrices by changing the "type"
% parameter in the script
% Note: For small matrices direct method would be faster however for larger 
% matrices the indirect method would be faster than the direct method
% Note: Indirect method is also amenable for large scale out of core parallel
% matrix multiplication (See manuscript for details) 
% Code Written By : Vineet Yadav
% Date: 9/06/2012
%%%%%%%%%%%%%%%%%%%%%%%CODE BEGINS%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc % Clear Console
clear all % Clear All Variables in the Workspace
%% MULTIPLICATION OF H and Q WHERE Q IS A KRONECKER PRODUCT: DIRECT METHOD
% Some Parameters
% Q can be decomposed as a kronecker product of two small matrices D and E
% where D has dimensions (p*q)and E has dimensions (r*t) thus Q has
% dimensions (pr*qt).
% H has dimensions n * m; where n is arbitrary and m should be equal
% to number of rows in Q.
%----------------Change these parameters to experiment with different size
% of matrices--------------------------------------------------------------
display('Warning: Demonstration Code Do Not Use Large Values For Dimensions of Matrices')
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
p = 35;
q = 35;
r = 70;
t = 70;
n = 500;
% NOTE:: IMP!! type parameter has two options 'HQ' or 'HQHT' and gives the
% user option to compute HQ or HQHT from the method described in Sect 2 of 
% the manuscript . However H*Q*H' is only defined when p=q and r=t
type='HQ'; % type has two options as HQ or HQHT 
% ----------------Do not change anything after this line-------------------
m = p*r;
% Draw normal random numbers with mean 1 and std 2
D = 0 + 1.*randn(p,q);
E = 0 + 1.*randn(r,t);
Q = kron(D,E);
H = rand(n,m);
%----------------Time of direct multiplication-----------------------------
a=tic;
if strcmpi (type,'HQ')
    HQ_DIRECT = H*Q;
    display(['Total Time Taken in "DIRECT" Multiplication of HQ is: ' num2str(toc(a)) ' Seconds'])
elseif strcmpi (type,'HQHT')
    HQHT_DIRECT = H*Q*H';
    display(['Total Time Taken in "DIRECT" Multiplication of HQHT is: ' num2str(toc(a)) ' Seconds'])
end
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
clear a % clear variables not required
%% COMPUTE H*Q or H*Q*H' WHERE Q IS A KRONECKER PRODUCT: INDIRECT METHOD
% Divide H in equal vertical strips; This operation can be also performed
% first by dividing H in strips and saving them on harddrive making
% multiplication of HQ amenable for distributed parallel computing
Hstrips=mat2cell(H,n,repmat(r,p,1)); % p columns divided into n*r matrices % SEE EQUATION 6 or EQ 8
% in dimensions
% H created earlier is not required as it has already been divided in
% strips
clear H
%--------------- Assign Space For Storing Output of HQ multiplication------
if strcmpi (type,'HQ')
    HQ_INDIRECT=zeros(n,q*t);
elseif strcmpi (type,'HQHT')
    HQHT_INDIRECT=zeros(n,n);
end
counter = t-1; % a counter variable
a=tic;
for i=1:q
    temp = i*t; % Just for indices computation
    HQsum = zeros (n,r);% Assign space for for storing block sum
    for j=1:p
        if D(j,i) ~= 0 && D(j,i) ~=1 % To avoid multiplying H column blocks by one or zero
            HQsum=HQsum+Hstrips{j}.*D(j,i);
        elseif D(j,i) == 1 % Only add if it is 1
            HQsum=HQsum+Hstrips{j};
        end
    end
    %---------------- Put matrix multiplication output at correct place----
    if strcmpi (type,'HQ')
        HQ_INDIRECT(:,temp-counter:temp)=HQsum*E; % EQUATION 7 or 9 IN THE MANUSCRIPT
    elseif strcmpi (type,'HQHT')
        HQHT_INDIRECT=HQHT_INDIRECT+HQsum*E*(Hstrips{i}'); % EQUATION 16 IN THE MANUSCRIPT
    end
end
if strcmpi (type,'HQ')
    display(['Total Time Taken in "INDIRECT" Multiplication of HQ is: ' num2str(toc(a)) ' Seconds'])
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display('All Variables are as defined in the manuscript Check HQ_direct and HQ_indirect for comparing answers')
elseif strcmpi (type,'HQHT')
    display(['Total Time Taken in "INDIRECT" Multiplication of HQHT is: ' num2str(toc(a)) ' Seconds'])
    display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
    display('All Variables are as defined in the manuscript Check HQHT_direct and HQHT_indirect for comparing answers')
end
clear a counter HQsum Hstrips i j m n p q r t temp % clear variables
