%% Initialize parameters
%%short  -- 100copies_2

clear all;
set(0,'DefaultFigureWindowStyle','docked')
%amp = 0.0261;
numRuns = 100;
numCopies = 100;
iOrder = 1;
%% Read strain data and get true singular values of augmented strain matrix
% Read strain matrices
load('Eps1.mat');
load('Eps2.mat');
load('Eps6.mat');

% Compose row-augmented matrix and compute SVD
A = [E1Matrix; E2Matrix; E6Matrix]; % row augmentation
[m,n]=size(A);
if (m < n) 
    A = A';
    iflip = 1;
end
[m,n]=size(A);
[U, S, V] = svd(A);

% Plot the true singular values 
maxRank = min(m,n);
trueSingVal = diag(S(1:maxRank, 1:maxRank));
% figure; loglog(trueSingVal);

%% Partition U, S, V
m_1 = 1;  % number of SVs of interest
m_2 = n - m_1;
% U
U1 = U(:, 1:m_1);
U2 = U(:, m_1+1:n);
U3 = U(:, n+1:m);
% S
S1 = S(1:m_1, 1:m_1);
S2 = S(m_1+1:n , m_1+1:n);
%
V1 = V(:,1:m_1);
V2 = V(:,m_1+1:n);
%
u1 = U(:,1);
v1 = V(:,1);
%% Compute 0th order perturbation
% initialize $x_i$, $y_i$ and $\gamma_i$
sigma1 = S1(1:1);
% generate Error matrix
% rng('default');
% normE = sqrt(m) + sqrt(n);
%amp = norm(A)/normE;

sigNoiseRatio = [200 180 160 140 120 100 80 60 40 20 0];
Doc_rightEV = zeros(size(U,2),size(sigNoiseRatio,2),numCopies);
Doc_rightEV_true = zeros(size(U,2),size(sigNoiseRatio,2),numCopies);
% S3 = RandStream.create('mt19937ar','NumStreams',1,'StreamIndices',1);
% RandStream.setGlobalStream(S3);
for counter_copies = 1 : 1 :numCopies
	
	counter=1;
	for counter=1:size(sigNoiseRatio,2)
		[counter_copies,counter]
		 
		[Atilde,E]=add_wgn(A,sigNoiseRatio(counter));
		 
		[Utrue,Strue,Vtrue] = svd(Atilde);  
		% initialize iterates
		x0 = v1;
		gamma0 = sqrt(norm(sigma1 + u1'*E*v1)^2 + norm(U3'*E*v1)^2);
		% gamma1 = gamma0;
		y0 = u1*gamma0;

		%% get higher order perturbation estimates 
		u1 = U(:,1);
		v1 = V(:,1);
		xvec = zeros(n,iOrder+1);
		yvec = zeros(m,iOrder+1);
		uvec = zeros(m,iOrder+1);
		gammavec = zeros(iOrder+1,1);
	   
		%
		xvec(:,1) = x0;
		yvec(:,1) = y0;
		gammavec(1) = gamma0;
		uvec(:,1)=y0/gamma0;
		%
		% figure; 
		%
		for i = 1:iOrder
			%
			j = i + 1;
			x = xvec(:,i);% right
			y = yvec(:,i);% left
			gamma = gammavec(i);
			%
			temp = gamma^2*eye(size(S2)) - S2*S2;
			temp1 = temp\(V2'*E'*y + S2*U2'*E*x);
			VHx = [sqrt(1 - norm(V2'*x)^2); temp1];
			xvec(:,j) = V*VHx;
			%
			temp1 = sqrt(1 - norm(V2'*x)^2)*sigma1 + u1'*E*x;
			temp2 = temp\(S2*V2'*E'*y + gamma^2*U2'*E*x);
			temp3 = U3'*E*x;
			UHy = [temp1; temp2; temp3];
			yvec(:,j) = U*UHy;
			%
			gammavec(j) = norm(Atilde*xvec(:,j));   
			uvec(:,j) = yvec(:,j)/gammavec(j);

		end
		
		
		

		
		Doc_rightEV(:,counter,counter_copies) = uvec(:,2); %Since A=A' we are using uvec instead of v.
        Doc_rightEV_true(:,counter,counter_copies) = Utrue(:,1); %First column of U here is the first right eigen vector of the untransposed row-augmented matrix
		
	end
end

save('100copies_11NLevels_2.mat','U','sigNoiseRatio','Doc_rightEV','Doc_rightEV_true');




    
    
    

    


