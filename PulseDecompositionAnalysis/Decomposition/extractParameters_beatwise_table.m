clear all
close all
clc

%% Paths
% add path to functions that are required (derivation, adding noise,...)
addpath('..\NeededFunctions');
% add path to decomposition functions and cell with algorithm names
addpath('..\Algorithms');
% add path to parameter functions and cell with parameter names
addpath('..\Parameters');
% load file containing algorithms here
load('algorithmsSORELLI.mat','algorithms');
% load file containing parameters here
load('parameters.mat','parameters');
% check that there is an algorithm for every entry in the list
for actualAlgorithm = 1:size(algorithms,1)
    if(~(isfile(['..\Algorithms\' 'decompose' algorithms{actualAlgorithm} '.m'])))
        errordlg('Unknown decomposition algorithm selected','Input Error','modal');
        return
    end
end
% specify source, data and results folder
dataset ='CPT';
sourceFolder=['Data\' dataset '\realData\'];
resultsFolder=['Data\' dataset '\decompositionParameters_beatwise\'];
dataFolder=['Data\' dataset '\results_beatwise\'];
% load data information
load([sourceFolder 'physiologicalMeasuresTable.mat']);
load([sourceFolder 'epochs.mat']);
patients=physiologicalMeasuresTable.SubjectID; % loads list with patients
% initialize struct in which data is stored
parameterTables = struct;
parameterTable = table;
parameterTable.SubjectID = patients;

%% Extract parameters
for actualAlgorithm = 1:size(algorithms,1)
    for actualPatientNumber=1:size(patients,1)
        for currentInterval=1:size(epochs,1)
            currentFile = [dataFolder patients{actualPatientNumber} '\' epochs{currentInterval} '\' algorithms{actualAlgorithm}];
            % try loading current file
            try
                load(currentFile)
            catch
                for actualParameter = 1:size(parameters,1)
                    parameterTable.([epochs{currentInterval} '_' parameters{actualParameter,1}])(actualPatientNumber) = NaN;
                end
                continue
            end
            currentParameterList = zeros(size(beatIndicesEnsembleBeat,1),1);
            %% calculate parameters
            for actualParameter = 1:size(parameters,1)
                for beatNumber = 1:size(beatIndicesEnsembleBeat,1)
                    try
                        currentParameterList(beatNumber) = feval(['calculate_' parameters{actualParameter,1}], ...
                            decompositionResults(beatNumber).signal_mod, ...
                            decompositionResults(beatNumber).y, ...
                            decompositionResults(beatNumber).opt_params, ...
                            algorithms{actualAlgorithm}, ...
                            samplingFreq);
                    catch
                        currentParameterList(beatNumber) = NaN;
                    end
                end
                parameterTable.([epochs{currentInterval} '_' parameters{actualParameter,1}])(actualPatientNumber) = mean(currentParameterList,'omitnan');
            end
        end
    end
    parameterTables.(algorithms{actualAlgorithm}) = parameterTable;
end
%% save parameters
if(exist(resultsFolder,'dir')~=7)
    mkdir(resultsFolder)
end
save([resultsFolder 'parameterTables.mat'],...
    'parameterTables');