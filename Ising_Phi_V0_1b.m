%% Avg. Phi vs. Temperature in the Ising Model
%
% Monte Carlo Methods for the evaluation of the Phi and thermodynamics
% properties of the one-dimensional Generalized Ising Model.
%
clc; clear all; close all;
tic
%% ------------------------------------------------------------------------
% Loop Parameters

timeSteps = 250;

T_start = 0.01; % Initial temparature
T_end = 2;   % Final temparature
dT = 0.01; % Temperature steps

numBin = 50; % # of bin
% ------------------------------------------------------------------------
%--------------------------- Network Properties --------------------------
% 
% CONVENTION: 
% This code uses the convention that spin up is represented by
% 0 and spin down is represented by 1. The matrices used in calculating the
% TPM the dE of the system and other network related matrices follow the
% Lowest Order Lowest Index (LOLI) convention.
N = 3;

% ------------------------- Connectivity Matrix  -------------------------

spars = 0; % 0 means full connected, 1 means no connections
antiCorSpars = 1; % 0 means all anticorrelations, 1 means no anticorrelations
J = (rand(N) > spars);
J = J.*((rand(N) > 1- antiCorSpars)*2 -1); % Includes anti correlations



% Make sure either spars and antiCorsSpars are not both 1. Generates
% non-trivial Connectivity Matrix for small N.

while (sum(sum(J)) == 0)
    J = (rand(N) > spars);
    J = J.*((rand(N) > 1- antiCorSpars)*2 -1); % Includes anti correlations
    

    
    J = triu(J) + triu(J,1)';
    J(1:N+1:end) = 0;
end
% The code below ensures: that connections have real values and are not
% simply binary; that the matrix is symmetric about the diagonal, and its
% diagonal is 0.
J = J.*rand(N);
J = ones(N);
J = triu(J) + triu(J,1)'; % make sure the matrix is symmetric
J(1:N+1:end) = 0; % make sure the diagonals are 0
% ------------------------------------------------------------------------
% --------------------------- Matrices for Phi ---------------------------
% Complete Repertoire
M = ((dec2bin(0:(2^N)-1)=='0') - 0)*2 - 1;
PF = zeros(length(M),1);

% dE and E of complete repertoire
M = fliplr(((dec2bin(0:(2^N)-1)=='0') - 0)*2 - 1);
dE_sys = M*J.*M*2;

Jtemp = J;
Jtemp(logical(triu(ones(N)))) = 0;
ener_sys = sum(M*(-Jtemp).*M,2);
clear Jtemp

% Matrices used in calculating the TPM
M_TPM = fliplr((dec2bin(0:(2^N)-1)) == '1');

% Generate Matrix Indicies
PartIndex = 1:N;
% Generates all possible combinations of partitions from size 2:(MAX-1)
for i = 1:(size(PartIndex,2)-1)
    combinationsSet{i} = nchoosek(PartIndex,i);
    % NOTE: nchoosek function valid when length(PartIndex) < 15.
end
clear PartIndex

numPartSet = size(combinationsSet,2);
% -------------------------------------------------------------------------

% Initializations of Matrices and Loop Variables
% -------------------------------------------------------------------------
temp = T_start:dT:T_end;
temp_len = length(temp);
% B = B_start:dB:B_end;
% beta = 1./B;

phi_mean = zeros(numBin,1);
ener_mean = zeros(numBin,1);
mag_mean = zeros(numBin,1);
mag2_mean = zeros(numBin,1);
ener_sqr_mean = zeros(numBin,1);
mag_sqr_mean = zeros(numBin,1);
specheat_mean = zeros(numBin,1);
sus_mean = zeros(numBin,1);

Ener = zeros(1, temp_len);
Mag = zeros(1, temp_len);
% Spec_Heat= zeros(1, temp_len);
Sus = zeros(1, temp_len);
Phi = zeros(1, temp_len);
TPM = zeros(2^N);

count = 0;
totalCount = 0;

for T = temp;
    
    % Generate a random starting state
    count = count + 1;

    
    clear FPMtemp iTPM jTPM lgclP
    
    % For each bin , we perform no_flip number of flips
    for b = 1 :numBin;
        totalCount/(numBin*length(temp))*100
        totalCount = totalCount + 1;
        
        stateVec = (rand(N,1) > 0.5)*2 - 1;
        currentState = bi2de(logical(stateVec' == -1));
        
        mag = sum(stateVec);
        
        stateDist = zeros(1,N^2);
        close all
        
        % This generates the Flip Probablity Matrix (FPM) for the system. It is
        % written in the LOLI (Low-Order bits correspond to Low-Index nodes)
        % convention as described in the IIT 3.0 python documentation
        % (https://pythonhosted.org/pyphi/examples/2014paper.html). This matrix
        % is then used to generate the Transition Probability Matrix (TPM) for
        % each temperature. This is essentially the entire Metropolis Algorithm
        % in this line.
        detFlip = (dE_sys <= 0);
        FPM = dE_sys; FPM(detFlip) = 1; FPM(~detFlip) = exp(-FPM(~detFlip)/T);
        
        for iTPM = 1:2^N
            for jTPM = 1:2^N
                
                lgclP = logical(M(iTPM,:) - M(jTPM,:));
                FPMtemp = FPM(iTPM,:);
                FPMtemp(~lgclP) = 1 - FPMtemp(~lgclP);
                TPM(iTPM,jTPM) = prod(FPMtemp);
                if sum(sum(isnan(TPM))) > 0
                    
                    pause()
                end
                
            end
        end
        
        ener_sum = 0;
        phi_sum = 0;
        phi_sum2 = 0;
        ener_sqr_sum = 0;
        mag_sum = 0;
        mag_sqr_sum = 0;
        
        %%% NEST THIS LOOP IN ANOTHER TIMESTEP LOOP AND MAKE THIS LOOP RUN
        %%% IN SUCH A WAY THAT THE NODES DON'T CHANGE INSTANTANEOUSLY IN A
        %%% TIMESTEP.
        
        %%% QUESTION: What is the least amount of time until the system
        %%% equilibrates? Need to setup some tests to solve this one.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for ttStep = 1:timeSteps
            stateIndex = currentState + 1;
            flipVec = ones(1,N);
            
            % Go through the system and check for flips one by one. The
            % decision to flip is not effective until the entire loop is
            % complete.
            TESTRUN = randi([1 2],1,100);
            for nodeIndex = 1:N;
                % Compute the change in energy
                dE_flipEl = dE_sys(stateIndex,nodeIndex);
                
                if (dE_flipEl <= 0);
                    flipVec(nodeIndex) = -1;
                    %                 mag = mag + 2*spin_vec(Flip);
                    %  elseif (rand <= exp(-dE/T));
                    %     elseif (rand <= exp(-dE/(p(Flip,1)*T)));
                elseif (rand <= exp(-(dE_flipEl)/T));
                    flipVec(nodeIndex) = -1;
                    %                 spin_vec(flipEl) = - spin_vec(flipEl); mag = mag +
                    %                 2*spin_vec(Flip);
                end
            end
            
            % Update the state of the system according to the flips that
            % occured in this timestep.
%             stateVec
%             flipVec
            
            stateVec = stateVec.*flipVec';
%             stateVec
            currentState = bi2de(logical(stateVec' == -1));
%             stateDist(currentState + 1) = stateDist(currentState + 1) + 1;
%             bar(stateDist)
%             drawnow;
%             
            % ------------------------- PHI -----------------------------
            pCount = 0;
            % Runs through all possible bi-partition pairs of network to
            % calculate Phi.
            for i_part = 1:ceil(numPartSet/2)
                partLength = length(combinationsSet{i_part});
                if i_part == ceil(numPartSet/2) && mod(numPartSet,2) ~= 0
                    jLoop = partLength/2;
                else
                    jLoop = partLength;
                end
                
                % Generate Phi for this bi-partition pair
                for j_part = 1:jLoop
                    pCount = pCount + 1;
                    part1 = combinationsSet{i_part}(j_part,:);
                    part2 = combinationsSet{numPartSet - (i_part - 1)}(partLength - (j_part - 1),:);
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%% TPM ALGORITHM %%%%%%%%%%%%%%%
                    % For each row in the whole system's cause repertoire,
                    % generates indices representing to which index in the
                    % partitioned cause repertoire they fall under.
                    % REPRESENTS HOW THE ROWS OF THE TPM SHOULD BE AVERAGED
                    indP1 = bi2de(M_TPM(:,part1));
                    indP2 = bi2de(M_TPM(:,part2));
                    
                    % Logical array representing all whole system states
                    % that are equal to the partitioned current state.
                    % REPRESENTS COLUMNS OF THE TPM The Columns of the TPM
                    % that correspond to the partitioned current state
                    indTPV1 = bi2de(M_TPM(:,part1)) == bi2de(stateVec(part1)' == -1);
                    indTPV2 = bi2de(M_TPM(:,part2)) == bi2de(stateVec(part2)' == -1);
                    
                    % mean changed to sum
                    TPM1 = mean(TPM(:,indTPV1),2);
                    %                     TPM1 = TPM1./max(TPM1)
                    TPM2 = mean(TPM(:,indTPV2),2);
                    %                     TPM2 = TPM2./max(TPM2);
                    
                    % Transition probabily vector (TPV): Cause repertoire
                    % of first and second partition
                    for iTPV1 = 1:2^length(part1)
                        TPV1(iTPV1) = mean(TPM1((indP1 == (iTPV1 - 1)),:));
                    end
                    TPV1 = TPV1./(sum(TPV1));
                    for iTPV2 = 1:2^length(part2)
                        TPV2(iTPV2) = mean(TPM2((indP2 == (iTPV2 - 1)),:));
                    end
                    TPV2 = TPV2./(sum(TPV2));
                    %                     sum(TPV1) sum(TPV2)
                    
                    for iFinal = 1:2^N
                        TPV(iFinal) = TPV1(indP1(iFinal)+1) * TPV2(indP2(iFinal)+1);
                    end
                    
                    
                    
                    %------------------------ Kullback-Leibler Divergance
                    %---------------------
                    %--------------------------------------------------------------------------
                    % % TO DO: IIT 3.0 uses EMD (Earth Mover's Distance).
                    % Might solve Phi > N % problem when N = 2? Check this
                    % out.
                    %--------------------------------------------------------------------------
                    
                    PTPM = TPM(:,currentState+1)./sum(TPM(:,currentState+1));
                    
                    H2 = PTPM.*log2(PTPM./TPV');
                    H2(isnan(H2)) = 0;
                    H2 = sum(H2);
                    effInfo2(pCount) = H2;
                    H2 = H2/(min([length(part1) length(part2)]));
                    effInfoNorm2(pCount) = H2;
                end
            end
            % ----------------------- END PHI ----------------------------
            ener = ener_sys(stateIndex);
            mag = sum(stateVec);
            
            phi_sum = phi_sum + min(effInfoNorm2);
            ener_sum = ener_sum + ener;
            ener_sqr_sum = ener_sqr_sum + ener^2;
            mag_sum = mag_sum + mag;
            mag_sqr_sum = mag_sqr_sum + mag^2;
        end
        phi_mean(b) = phi_sum/timeSteps;
        
        ener_mean(b) = ener_sum/timeSteps; 
        mag_mean(b) = mag_sum/timeSteps;
        ener_sqr_mean(b) = ener_sqr_sum/timeSteps; 
        mag_sqr_mean(b) = mag_sqr_sum/timeSteps; 
        % specheat_mean(b) = (ener_sqr_mean(b)-(ener_mean(b))^2)/N/T^2;
        sus_mean(b) = (mag_sqr_mean(b)-(mag_mean(b))^2)/N/T;
    end
    
    
    Phi(count) = mean(phi_mean);
    Ener(count) = mean(ener_mean)/N ; Mag(count) = mean(mag_mean)/N ;
    Sus(count)= mean(sus_mean)/N ;
    mm = 0;
    
%     figure
%     imagesc(TPM)
    
end
clearvars -except Phi Ener Sus Mag temp

figure
plot(temp,Sus,'ro')
title('Susceptability vs. Temperature');
xlabel('T');
% Plot Phi ------------------
figure
plot(temp,Phi,'ro')
title('Phi2 vs. Temperature');
xlabel('T');

figure
plot(temp,Ener,'ro')
title('Energy vs.Temperature');
xlabel('T'); ylabel('E');

toc
