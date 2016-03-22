% Given its cause repertoire (M), its partition cause repertoir (Mpart),
% its temperature (T), and the partition indicies (A), and its (simplified)
% transition probability matrix, output the full transition probability
% matrix for the (partitioned) system to enter the state defined by MD2. 

function trans_P = transP_Par2(M,Mpart,MD2,dE,T,A)

k = 1;
P = ones(size(M,1),1);

% The matrix MD2 has three possible values:

% 0: represents nodes where it will deterministically flip towards final
% state.

% 1: represents nodes where it will probabilistically flip in either
% direction (this will later be distinguished between favorable and
% unfavorable flips within the loop below when it is multiplied by .*M.)

% 2: represents nodes where it will deterministically flip away from final
% state (therefore all rows with 2 entries have 0% probability of
% transition into final state)


% MD2 = MD - ones(length(M),1)*F;


for jj = 1:length(M)
    % deterministic impossibilities
    if or(sum(MD2(jj,A) == 2),sum(MD2(jj,A) == -2 )) > 0
%         MD2(jj,:) = 0;
        P(jj) = 0;
    else
        MD2(jj,:) = MD2(jj,:).*M(jj,:);
        
        for ii = A
            if MD2(jj,ii) > 0 % favorable flips
                P(jj) = P(jj)*exp(-dE(jj,ii)/(k*T));
            elseif MD2(jj,ii) < 0 % unfavorable flips
                P(jj) = P(jj)*(1 - exp(-dE(jj,ii)/(k*T)));
%             else % (MD == 3)
%                 P(jj) = P(jj)*(0.5);
            end            
        end
        
    end
end
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       
trans_P = zeros(length(Mpart),1);

for i = 1:length(Mpart)
%     partIndex = find(ismember(M(:,A),Mpart(i,:),'rows'));
    partIndex = sum(M(:,A) - ones(length(M),1)*Mpart(i,:) == 0,2) == length(A);
    trans_P(i) = sum(P(partIndex));
    
end
    


trans_P = trans_P./sum(trans_P);

