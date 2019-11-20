% Submitted with the manuscript:
% A TECHNICAL NOTE ON IMPROVING COMPUTATIONAL EFFICIENCY IN LARGE LINEAR 
% INVERSE PROBLEMS: AN EXAMPLE FROM CARBON DIOXIDE FLUX ESTIMATION
% Vineet Yadav[1] and Anna M. Michalak[1]
% [1]{Department of Global Ecology, Carnegie Institution for Science, 
% Stanford, California, USA 94305}
% THIS IS MATLAB CODE THAT CAN BE RUN EITHER IN MATLAB OR OCTAVE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Purpose:To Demonstrate Algorithm Proposed in Section 3 of the manuscript
% ALL VARIABLES AND INDICES ARE AS DEFINED IN THE MANUSCRIPT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This code demonstrates fast computation and extraction of a posteriori 
% covariance for bayesian inversions through arbitrary random matrices. 
% Vshat=Q-(H*Q)'inv(H*Q*H'+R)*(H*Q))
% The validity of the proposed method can be checked by comparing answers 
% from the direct and indirect method i.e.the method proposed in the manuscript.
% In this script we have only shown how these computations are performed
% not the accuracy or resonableness of the values used in different
% matrices.
% Code Written By : Vineet Yadav
% Date: 09/06/2012
%% Clear Console and All the Variables From the Memory
clc
clear all
%% ----------------Define Some Parameters----------------------------------
%----------------Change these parameters to experiment with different size
% of matrices--------------------------------------------------------------
% NOTE: Always p=q and r=t as Q is symmetric a priori covariance matrix 
% see submitted manuscript for details
display('Warning: Demonstration Code Do Not Use Large Values For Dimensions of Matrices')
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
p = 35; % the dimension of the temporal covariance
q = 35; % the dimension of the temporal covariance
% NOTE: max time periods cannot be greater than p or q (p=q) or lower than
% 1
r = 70; % the dimension of the spatial covariance
t = 70; % the dimension of the spatial covariance
n = 500; % no of observations
% Define time periods for which you want to compute uncertainity
Begin_Time_Period = 1;
End_Time_Period = 10; % should not be greater than q
%% --------------Do not change anything after this line--------------------
m = p*r;
ms = r;
mt = p;
% Draw normal random numbers with mean 1 and std 2
D = 0 + 1.*randn(p,q); % Temporal Covariance
E = 0 + 1.*randn(r,t); % Spatial Covariance
Q = kron(D,E); 
H = rand(n,m);
R = diag(randn(n,1));
%% These costs are same in the direct and indirect method so we compute
% before a posteriori covariance
HQ = H*Q;
HQHTR=HQ*H'+R;
IHQHTR=inv(HQHTR);
%% -------------Directly Compute the Requested Uncertainity----------------
% After computing full uncertainity
a=tic;
Vshat = Q-HQ'*IHQHTR*HQ; %#ok<MINV>
% A posteriori covariance for an Inversion over requested time periods
% sum of all blocks and then divison of all the summed blocks
for i = Begin_Time_Period:End_Time_Period
    for j=Begin_Time_Period:End_Time_Period
        if i == Begin_Time_Period && j == Begin_Time_Period
            Vshat_Direct = Vshat(((j-1)*ms)+1:ms*j,...
                ((i-1)*ms)+1:ms*i);
        else
            Vshat_Direct = Vshat_Direct+Vshat(((j-1)*ms)+1:ms*j,...
                ((i-1)*ms)+1:ms*i);
        end
    end
end
Vshat_Direct = Vshat_Direct./(mt^2);
display(['Total Time Taken for "DIRECTLY" computing a posteriori covariance is: ' num2str(toc(a)) ' Seconds'])
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
%% ---------------Indirectly Compute A Posteriori Covariance---------------
% This method indirectly computes a posteriori covariance by reducing both the 
% memory and computational costs and even allows performing inversions on a
% desktop. 
% Sum of HQ
a=tic;
% EQUATION 17 IN THE MANUSCRIPT
Dsum=D(Begin_Time_Period:End_Time_Period,Begin_Time_Period:End_Time_Period);
Dsum=sum(sum(Dsum));
Q_Sum=Dsum*E;
%%% EQUATION 18 IN THE MANUSCRIPT
for i = Begin_Time_Period : End_Time_Period
    if i == Begin_Time_Period 
        HQ_Sum = HQ(:,((i-1)*ms)+1:ms*i);
    else
        HQ_Sum = HQ_Sum+HQ(:,((i-1)*ms)+1:ms*i);
    end
end
%%% EQUATION 19 IN THE MANUSCRIPT
Vshat_Indirect = Q_Sum-(HQ_Sum)'*IHQHTR*HQ_Sum; %#ok<MINV>
Vshat_Indirect = Vshat_Indirect./(mt^2);
%%% DONE
display(['Time Taken for "INDIRECTLY" computing a posteriori covariance is: ' num2str(toc(a)) ' Seconds'])
clear a i j m ms mt n p q r t  Begin_Time_Period End_Time_Period Dsum Q_Sum...
    HQHTR HQ_Sum IHQHTR 
display('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
display('All Variables are as defined in the manuscript Check Vshat_direct and Vshat_indirect for comparing answers')
%----------------------------END-------------------------------------------