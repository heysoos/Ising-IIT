clc; clear all; close all

load indice_map.mat
network = {'Brain','Auditory', 'DMN', 'ECN', 'Salience', 'Sensorimotor', 'Visual'};

% Change 'work_dir' according to the present work directory
folder = 'Brain Data/Propofol_timeseries/';
folderList = dir(folder);

brainStates = {'S1','S2','W1','W2'};
% w1 awake
% s1 mild sedation
% s2 unconsciousness
% w2 recovery

Subject = 'Subject';
BrainState = 'BrainState';
Network = 'Network';
TimeSeries = 'TimeSeries';
Correlation = 'Correlation';
ogTime = 'ogTime';

dataEntryCount = 0;

% cycle through all folders
for i = 3:length(cat(1,folderList.isdir))
    subName = folderList(i).name;
    subFolder = strcat([folder,subName,'/']);
    % cycle through 4 brain states
    for j = 1:length(brainStates)
        % generate filename
        fileDir = strcat([subFolder,subName,char(brainStates(j)),'.mat']);
        
        % if the file exists
        if numel(dir(fileDir)) ~= 0
            
            temp = load(fileDir,'func_mean');
            timeSeries = temp.func_mean;
            
            % using the indice maps average the timeseries and generate a
            % correlation matrix. Need to cycle through the 6 networks and
            % then the whole brain.
            % Outputs 7 correlations.
            
            % FUNCTION: indice_map = INDICE_MAP_GENERATOR(timeSeriesCorr,MJ_all_mn,region);
            
            for ii = 1:(size(indice_map,1)+1) % + 1 for whole system
                dataEntryCount = dataEntryCount + 1;
                tempSeries = [];
                close all
                if ii == 1
                    % ii = 1 corresponds to the full brain network.
                    % Generates timeseries by averaging through the
                    % 6 networks as defined by the indice map.
                    
                    for iind = 1:size(indice_map,1)
                        tempSeries(iind,:) = mean(timeSeries(indice_map(iind,:),:));
                        
                    end
                else
                    
                    tempSeries = timeSeries(indice_map(ii-1,:),:);
                end
                
                %----------------------- CORRELATION GENERATION --------------------------%
                corrFinal = corrcoef(tempSeries');
                corrFinal(1:(length(indice_map)+1):end) = 0;
                
                data(dataEntryCount).(Subject) = subName;
                data(dataEntryCount).(BrainState) = brainStates{j};
                data(dataEntryCount).(Network) = network{ii};
                data(dataEntryCount).(TimeSeries) = tempSeries;
                data(dataEntryCount).(Correlation) = corrFinal;
                data(dataEntryCount).(ogTime) = timeSeries;

%                 display([subName, brainStates{j},'-', network{ii}])
                
            end
            
            
        end
    end
    
end
clearvars -except data

% [DATA, FUNC(DATA)] = [STATE, NETWORK, DATA]
state = 'S2';
network = 'Brain';
DATA = 'ogTime';

[extract, source, ind] = structTimeCorr(data, 'BrainState', state, 'Network', network, DATA);

LINSPEC = linspecer;
for i = 1:size(extract,3)
    figure(1)
    imagesc(extract(:,:,i))
    colormap(linspecer)
    colorbar
    caxis([-1 1])
    figure(2)
    plot(source(:,:,i)')
    axis tight
    figure(3)
    hist(extract(:,:,i)')
    pause()
end


