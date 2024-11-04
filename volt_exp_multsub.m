function volt_exp_multsub(task,datapath,IDs)
%%volt_exp_multsub
% calculates a Volterra expansion of the motion parameters (from the
% realingment procedure) 

% scale MP?
scale_MP = 0;


%% Loop to load data from folders and run the job
for i = 1:numel(IDs)

    MPFile = cfg_getfile('FPListRec',fullfile(datapath,num2str(IDs(i)),[task{1} '_EPI']),'^rp_.*.txt');

    if size(MPFile,1)>1
        MPFile(2,:) = [];
    end

    if isempty(MPFile)
        warning('No rp-file for data set: %d\n',IDs(i))
    end

    MP = load(cell2mat(MPFile));

    if scale_MP
        % scale MP (mean correct and divide by standard deviation)
        for CurrCol = 1:size(MP,2)
            MP(:,CurrCol) = MP(:,CurrCol) - mean(MP(:,CurrCol));
            MP(:,CurrCol) = MP(:,CurrCol) ./ max(eps,std(MP(:,CurrCol)));
        end
    end

    MP1stderiv = [zeros(1,6); diff(MP)];
    %MP2ndderiv = [zeros(1,6); diff(MP1stderiv)];
    VoltxMP = [MP MP1stderiv MP.^2 MP1stderiv.^2];
    %VoltxMP2 = [MP MP1stderiv MP.^2 MP1stderiv.^2 MP2ndderiv MP2ndderiv.^2];

    [P,N,E] = fileparts(cell2mat(MPFile));
    save(fullfile(P,['Vx',N,E]),'VoltxMP','-ascii');
    %save(fullfile(P,['Vx2',N,E]),'VoltxMP2','-ascii');

    %% Clean up
    clear MP MP1stderiv VoltxMP % MP2ndderiv VoltxMP2

end

