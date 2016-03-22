% Monte Carlo Methods for the evaluation of the thermodynamics properties
% of the one-dimensional Ising model without taking the nearest-neighbour
% coupling
%
clc; clear all; close all;
tic 
% ------------------------------------------------------------------------
% Simulation Parameters
N = 6;

T_start = 0.01; % Initial temparature
T_end = 5;   % Final temparature
dT = 0.01; % Temperature steps

numBin = 1; % # of bin
timeCorr = 201;
% ------------------------------------------------------------------------
% ------------------------- Generating Networks --------------------------

% ------------------------- Connectivity Matrix  -------------------------
spars = 0.3; % 0 means full connected, 1 means no connections
antiCorSpars = 1; % 0 means all anticorrelations, 1 means no anticorrelations

J = (rand(N) > spars);
J = J.*((rand(N) > 1- antiCorSpars)*2 -1); % Includes anti correlations
J = triu(J) + triu(J,1)'; % make sure the matrix is symmetric
J(1:N+1:end) = 0; % make sure the diagonals are 0

% Make sure either spars and antiCorsSpars are not both 1.
while (sum(sum(J)) == 0) % Generate non-trivial Connectivity Matrix for small N
    J = (rand(N) > spars);
    J = J.*((rand(N) > 1- antiCorSpars)*2 -1); % Includes anti correlations

    %     J = J.*rand(N);

    J = triu(J) + triu(J,1)';
    J(1:N+1:end) = 0;
end

J = J.*rand(N); % This line ensures the values in J are not just 0/1 but
% a number between -1 and 1. Comment if you need integer values.

% fileDir = 'Masters Simulations/Brain/';
% fileLoad = 'J_Brain.mat';
% fileSave = 'Ising_Phi_Brain.mat';
% 
% J = load([fileDir,fileLoad],'J_brain');
% J = J.J_brain;

% ------------------------------------------------------------------------
% Matrices for Phi

M = fliplr(((dec2bin(0:(2^N)-1)=='0') - 0)*2 - 1);
Mtest = fliplr((dec2bin(0:(2^N)-1)) == '1');
dE_sys = M*J.*M*2;
MD = M.*(((dE_sys <= 0 ).^2 ==1).*(-1));

TPM = zeros(2^N);
PF = zeros(length(M),1);

% ------------------------------------------------------------------------
% Generate Partition Indicies

PartIndex = 1:N;
% Generates all possible combinations of partitions from size 2:(MAX-1)
for i = 1:(size(PartIndex,2)-1)
    combinationsSet{i} = nchoosek(PartIndex,i);
    % NOTE: nchoosek function valid when length(PartIndex) < 15.
end
clear PartIndex

numPartSet = size(combinationsSet,2);

% Jtemp = J;
% Jtemp(logical(triu(ones(N)))) = 0;
% ener_sys = sum(M*(-Jtemp).*M,2);
% clear Jtemp

% ------------------------------------------------------------------------
% -------------------- Initialization of Matrices-- ----------------------
phi_mean = zeros(numBin,1);
ener_mean = zeros(numBin,1);
mag_mean = zeros(numBin,1);
mag2_mean = zeros(numBin,1);
ener_sqr_mean = zeros(numBin,1);
mag_sqr_mean = zeros(numBin,1);
specheat_mean = zeros(numBin,1);
sus_mean = zeros(numBin,1);
phi_data = [];

t_eq=N^2*100;
no_flip = t_eq; % No of flip at each temperature
temp = T_start:dT:T_end;
temp_len = length(temp);

Phi = zeros(1, temp_len);
Ener = zeros(1, temp_len);
Mag = zeros(1, temp_len);
Spec_Heat= zeros(1, temp_len);
Sus = zeros(1, temp_len);
% ------------------------------------------------------------------------

count = 0;
totalCount = 0;
randcoord = randperm(N);

for T = temp;
    %     T
    spin_vec = (rand(N,1) > 0.5)*2 - 1;
    count = count + 1;
    mag = sum(spin_vec);

    stateDist = zeros(1,N^2);
    close all

    % This generates the Flip Probablity Matrix (FPM) for the system. It is
    % written in the LOLI (Low-Order bits correspond to Low-Index nodes)
    % convention as described in the IIT 3.0 python documentation
    % (https://pythonhosted.org/pyphi/examples/2014paper.html). This matrix
    % is then used to generate the Transition Probability Matrix (TPM) for
    % each temperature.
    detFlip = (dE_sys <= 0);
    FPM = dE_sys; FPM(detFlip) = 1; FPM(~detFlip) = exp(-FPM(~detFlip)/T);
%     FPM = 1-M.*tanh(M*J'/T); FPM = FPM./(max(max(FPM))); % NEW TEST


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

    clear FPMtemp iTPM jTPM lgclP

    % For each bin , we perform no_flip number of flips
    for b = 1 :numBin;
        
        spin_vec = (rand(N,1) > 0.5)*2 - 1;
        
        totalCount/(numBin*length(temp))*100
        totalCount = totalCount + 1;

        ener2_sum = 0;

        ener_sum = 0;
        phi_sum = 0;
        phi_sum2 = 0;
        ener_sqr_sum = 0;
        mag_sum = 0;
        mag_sqr_sum = 0;

        nn = 0;

        %%% QUESTION: What is the least amount of time until the system
        %%% equilibrates? Need to setup some tests to solve this one.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        for i = 1:no_flip;
            % Run through the N points of the coordinate matrix once, and
            % when the counter is larger than N, it goes back to the first
            % point, and the coordinate is rerandomized before it is called
            % upon. This ensures the random but unique choice of index when
            % checking to flip.
            nn = nn + 1;
            if (nn > N);
                nn = nn-N;
                randcoord = randperm(N);
            end
            Flip = randcoord(nn);

            % Compute the change in energy

            dE=0;
            for j=1:N
                if (j~= Flip)
                    dE = dE + J(Flip,j)*spin_vec(j);
                end
            end
            dE=2*dE*spin_vec(Flip);
            
            if (dE <= 0);
                spin_vec(Flip) = - spin_vec(Flip);
                mag = mag + 2*spin_vec(Flip);
            elseif (rand <= exp(-dE/T));
                spin_vec(Flip) = - spin_vec(Flip);
                mag = mag + 2*spin_vec(Flip);
            end

            currentState = bi2de(logical(spin_vec' == -1));
            
            ener = 0;
            for ii = 1:N;
                for jj=ii+1:N
                    ener = ener - J(ii,jj)*spin_vec(ii)*spin_vec(jj);
                end
            end
            
            %           -------------- PHI ------------
            pCount = 0;

            % Generates the simplified TPM (transition probability matrix).
            % This matrix determines the nature of a transition from 3
            % groups for the Ising model: Deterministic flips towards
            % current state, deterministic flips away from current state,
            % and probabilistic flips in either direction. This matrix is
            % used for calculating the transition probability matrix in the
            % transP_Par2 function.

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
                    indP1 = bi2de(Mtest(:,part1));
                    indP2 = bi2de(Mtest(:,part2));

                    % Logical array representing all whole system states
                    % that are equal to the partitioned current state.
                    % REPRESENTS COLUMNS OF THE TPM
                    % The Columns of the TPM that correspond to the
                    % partitioned current state
                    indTPV1 = bi2de(Mtest(:,part1)) == bi2de(spin_vec(part1)' == -1);
                    indTPV2 = bi2de(Mtest(:,part2)) == bi2de(spin_vec(part2)' == -1);

                    TPM1 = mean(TPM(:,indTPV1),2);
                    TPM2 = mean(TPM(:,indTPV2),2);

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

                    for iFinal = 1:2^N
                        TPV(iFinal) = TPV1(indP1(iFinal)+1) * TPV2(indP2(iFinal)+1);
                    end



                    % Kullback-Leibler Divergance
                    %----------------------------------------------------
                    % % TO DO: IIT 3.0 uses EMD (Earth Mover's Distance).
                    % Might solve Phi > N % problem when N = 2? Check this
                    % out.
                    %-----------------------------------------------------

                    PTPM = TPM(:,currentState+1)./sum(TPM(:,currentState+1));

                    H2 = PTPM.*log2(PTPM./TPV');
                    H2(isnan(H2)) = 0;
                    H2 = sum(H2);
                    effInfo2(pCount) = H2;
                    H2 = H2/(min([length(part1) length(part2)]));
                    effInfoNorm2(pCount) = H2;


                    % TO DO: This normalization doesn't exist in IIT 3.0.

                end
            end
            phi_sum = phi_sum + min(effInfoNorm2);
            %            -------------- PHI ------------


            ener_sum = ener_sum + ener;
            mag_sum = mag_sum + mag;
            ener_sqr_sum = ener_sqr_sum + ener^2;
            mag_sqr_sum = mag_sqr_sum + mag^2;
            %             phi_sum = phi_sum + phi;
            
        end
       
        phi_mean(b) = phi_sum/no_flip;
        ener_mean(b) = ener_sum/no_flip;
        mag_mean(b) = mag_sum/no_flip;
        ener_sqr_mean(b) = ener_sqr_sum/no_flip;
        mag_sqr_mean(b) = mag_sqr_sum/no_flip;
        specheat_mean(b) = (ener_sqr_mean(b)-(ener_mean(b))^2)/N/T^2;
        sus_mean(b) = (mag_sqr_mean(b)-(mag_mean(b))^2)/N/T;
    end

    Phi(count) = mean(phi_mean);

    Ener(count) = mean(ener_mean)/N ;
    Mag(count) = mean(mag_mean)/N ;
    Spec_Heat(count)= mean(specheat_mean)/N ;
    Sus(count)= mean(sus_mean)/N ;
    
    mm = 0;

end
clearvars -except Phi Ener Mag Spec_Heat Sus S temp Corr_DTI J

subplot(3,2,1)
scatter(temp,Ener,'.')
title('Energy vs. Temp')
subplot(3,2,2)
scatter(temp,Mag,'.')
title('Magnetization vs. Temp')
subplot(3,2,3)
scatter(temp,Spec_Heat,'.')
title('Specific Heat vs. Temp')
ylim([-0.05 1])
subplot(3,2,4)
scatter(temp,Sus,'.')
title('Susceptibility vs. Temp')
ylim([-0.05 1])
subplot(3,2,[5 6])
scatter(temp,Phi,'.')
title('Phi vs. Temp')
ylim([-0.05 max(Phi)+1])