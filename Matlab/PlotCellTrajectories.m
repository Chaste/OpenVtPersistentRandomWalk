clear, close all

addpath ../../Chaste/anim/matlab/
addpath ../../OpenVT/Matlab

% Dirs to load

BaseType = 'OpenVT/Test01PersistentRandomWalk';

Models = {'Model004'}




for ModelIndex = 1:length(Models)
    Model = Models{ModelIndex};


    dir = [BaseType, '/', Model,]

    data = load(['../../../testoutput/',dir,'/results_from_time_0/results.viznodes']);



    for i=1:length(data(1,2:2:end))
        plot(data(:,2*i),data(:,2*i+1))
        hold on
    end
    
    



    SaveAsPngEpsAndFig(-1,['Figs/', Model], 12, 7/5, 12);
end