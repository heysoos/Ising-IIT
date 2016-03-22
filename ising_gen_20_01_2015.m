%=========================================================================
%
% Monte Carlo Methods for the evaluation of the thermodynamics properties
% of the one-dimensional Ising model without taking the nearest-neighbour
% coupling
% ========================================================================

clear all;

% load Corr_FMRI_DMN_6x6
% load DMN_con
% J = DMN2;
% J = J./max(max(J));
% clear DMN2
tic

N = 20; % Dimension of the connectivity matrix
temp = .0001:0.0001:3;
% B = 0.1:0.01:7;
% temp = 0.7:0.001:3;
% temp = 1./B;

time_corr= 0; %t_eq=m1*m2*10;
no_flip = 100*N^2; % No of flip at each temperature
no_bin = 1; % # of bin

% Connectivity Matrix
spars = 0.9; % 0 means full connected, 1 means no connections
antiCorSpars = 1; % 0 means all anticorrelations, 1 means no anticorrelations
J = (rand(N) > spars);
J = J.*((rand(N) > 1- antiCorSpars)*2 -1); % Includes anti correlations
J = J.*rand(N); % This line ensures the values in J are not just 0/1 but a number between -1 and 1.
% J = triu(J) + triu(J,1)'; % make sure the matrix is symmetric
J(1:N+1:end) = 0; % make sure the diagonals are 0

while (sum(sum(J)) == 0)
    J = (rand(N) > spars);
    J = J.*((rand(N) > 1- antiCorSpars)*2 -1); % Includes anti correlations
    J = J.*rand(N);
end
J = triu(J) + triu(J,1)';

% r2 = random(makedist('Weibull'),6);
% r2 = r2./max(max(r2));

% J = ones(N);
% J(1:N+1:end) = 0;
% 
% fileDir = 'Masters Simulations/Brain/';
% fileLoad = 'J_brain';
% 
% J = load([fileDir,fileLoad, '.mat'],fileLoad);
% J = J.(fileLoad);

J = load('J_7_cluster');
J = J.J;

[m1, m2]=size(J);

% A coordinate matrix is generated from random permutation of N length
randcoord = randperm(N);
count = 0;
temp_len = length(temp);
Ener = zeros(1,temp_len);
Mag = zeros(1, temp_len);
Spec_Heat= zeros(1, temp_len);
Sus = zeros(1, temp_len);
Corr_DTI = zeros(N,N, temp_len);
S = zeros(N,time_corr,temp_len);

totalCount = 0;

for T = temp;
    %     tic

    count = count + 1;

    for b = 1 :no_bin;
        spin_vec = (rand(N,1) > 0.5)*2 - 1;
        mag = sum(spin_vec);
        
        totalCount/(no_bin*length(temp))*100
        totalCount = totalCount + 1;
        
        ener_sum = 0;
        ener_sqr_sum = 0;
        mag_sum = 0;
        mag_sqr_sum = 0;
        nn = 0;
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
            ener = 0;
            for ii = 1:N;
                for jj=ii+1:N
                    ener = ener - J(ii,jj)*spin_vec(ii)*spin_vec(jj);
                end
            end
            ener_sum = ener_sum + ener;
            mag_sum = mag_sum + mag;
            ener_sqr_sum = ener_sqr_sum + ener^2;
            mag_sqr_sum = mag_sqr_sum + mag^2;
        end
        ener_mean(b) = ener_sum/no_flip;
        mag_mean(b) = mag_sum/no_flip;
        ener_sqr_mean(b) = ener_sqr_sum/no_flip;
        mag_sqr_mean(b) = mag_sqr_sum/no_flip;
        specheat_mean(b) = (ener_sqr_mean(b)-(ener_mean(b))^2)/N/T^2;
        sus_mean(b) = (mag_sqr_mean(b)-(mag_mean(b))^2)/N/T;
    end
    Ener(count) = mean(ener_mean)/N ;
    Mag(count) = mean(mag_mean)/N ;
    Spec_Heat(count)= mean(specheat_mean)/N ;
    Sus(count)= mean(sus_mean)/N ;
    %
    % Generating 198 time points for each of 1015 spinsites
    %
    mm = 0;
    for jj = 1:time_corr;
        for i = 1:N;
            % Compute the change in energy
            mm = mm + 1;
            if (mm > N);
                mm = mm-N;
                randcoord = randperm(N);
            end
            Flip = randcoord(mm);
            dE = 0;
            for j=1:N
                if (j~=Flip)
                    dE = dE + J(Flip,j)*spin_vec(j);
                end
            end
            dE=2*dE*spin_vec(Flip);
            % Decide whether the change is possible
            if (dE <= 0);
                spin_vec(Flip) = - spin_vec(Flip);
            elseif (rand <= exp(-dE/T));
                spin_vec(Flip) = - spin_vec(Flip);
            end
        end
        S(:,jj,count)=spin_vec;  % temperature x timepoints x spinsites
    end
    % Correlation matrix
    Corr_DTI(:,:, count)=corrcoef(squeeze(S(:,:, count))');
    Sus(count)= mean(sus_mean)/N ;
    display([num2str(count/length(temp)*100)])
    %toc
end

% Corr_DTI(isnan(Corr_DTI)) = 0;
% 
% for i = 1:size(Corr_DTI,3)
% corr_dist(i) = sum(sum(pdist2(squeeze(Corr_DTI(:,:,i)),Corr_FMRI2)));
% end
% 
% for i = 1:size(Corr_DTI,3);
% corr = Corr_DTI(:,:,i);
% corr1 = reshape(corr,[N*N,1]);
% corr2 = reshape(Corr_FMRI2,[N*N,1]);
% [~,~,ks2stat(i)] = kstest2(corr1,corr2); % Obtaining the test statistic for the 10 realizations as a function of temperature
% end

clear Flip J_max MJ_all T b count dE ener ener_mean ener_sqr_mean ener_sqr_sum ener_sum i ii j jj m1 m2 mag mag_mean mag_sqr_mean mag_sqr_sum mag_sum mm nn no_bin no_flip randcoord specheat_mean spin_vec sus_mean t_eq temp_len time_corr

% figure
% plot(temp,Ener,'ro')
% title('Energy vs. Temperature');
% xlabel('T');
% ylabel('E');
% % % Plot Magnetization
% % % ------------------
% figure
% plot(temp,Mag,'ro')
% title('Magnetization vs. Temperature');
% xlabel('T');
% ylabel('M');
% % % Plot Specific Heat
% % % ------------------
% figure
% plot(temp,Spec_Heat,'ro')
% title('Specific Heat vs. Temperature');
% xlabel('T');
% Plot Susseptiblity
% ------------------
% figure
% plot(1./temp,Sus,'ro')
% title('Susceptability vs. Temperature');
% xlabel('1/T');
% axis square
% figure
% plot(temp,Sus,'ro')
% title('Susceptability vs. Temperature');
% xlabel('T');
% axis square

% % Plot pdist2 Distance
% % ------------------
% figure
% plot(temp,corr_dist,'ro')
% title('Correlation Distance by pdist2');
% xlabel('T');
% axis square
% % % Plot kstest2 Distance
% % % ------------------
% figure
% plot(temp,ks2stat,'ro')
% title('Correlation Distance by kstest2');
% xlabel('T');
% axis square

% for i = 1:size(Corr_DTI,3)*10
%     figure(1)
%     imagesc(Corr_DTI(:,:,i))
%     xlabel(['T = ', num2str(temp(i))])
%     colormap(linspecer)
%     colorbar
%     caxis([-1 1])
%     figure(2)
%     hist(Corr_DTI(:,:,i)')
%     figure(3)
%     plot(S(:,:,i)')
%     ylim([-3 3])
%     pause()    
% end

subplot(2,2,1)
scatter(temp,Ener,'.')
axis tight
title('Energy vs. Temp')
subplot(2,2,2)
scatter(temp,Mag,'.')
axis tight
title('Magnetization vs. Temp')
subplot(2,2,3)
scatter(temp,Spec_Heat,'.')
axis tight
title('Specific Heat vs. Temp')
ylim([-0.05 1])
subplot(2,2,4)
scatter(temp,Sus,'.')
axis tight
title('Susceptibility vs. Temp')
ylim([-0.05 1.5])

toc
