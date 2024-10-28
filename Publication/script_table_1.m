for iAnimal = 1:length(ai)
    name = string(ai(iAnimal).name);
    if contains(name, "daisy")
        gender = 'F';
    end
    if contains(name, "desmond")
        gender = 'M';
    end
    id = str2double(regexp(name, '\d*', 'Match'));
    ai(iAnimal).displayName = sprintf('%i (%s)', iAnimal, gender);

    switch ai(iAnimal).probe
        case 'tetrode'
            ai(iAnimal).probeDisplayName = 'Chronic NiCr 32-wire bundle';
        case 'double-bundle'
            ai(iAnimal).probeDisplayName = 'Chronic NiCr 32-wire bundle';
        case 'bundle'
            ai(iAnimal).probeDisplayName = 'Chronic NiCr 32-wire bundle';
        case '4shank-neuronexus'
%             ai(iAnimal).probeDisplayName = 'Chronic NeuroNexus A4x8-7mm-100-200-177-CM32';
            ai(iAnimal).probeDisplayName = 'Chronic NeuroNexus A4x8';
        case '4shank-acute-wide'
            ai(iAnimal).probeDisplayName = 'Acute UCLA 128D';
        case '4shank-acute'
            ai(iAnimal).probeDisplayName = 'Acute UCLA 128DN';
    end

    ai(iAnimal).nUnits = nnz(strcmpi(eu.getAnimalName, name));
    ai(iAnimal).nSessions = length(unique({eu(strcmpi(eu.getAnimalName, name)).ExpName}));
    ai(iAnimal).nUnitsHasPress = nnz(strcmpi(eu.getAnimalName, name) & c.hasPress);
    ai(iAnimal).nSessionsPress = length(unique({eu(strcmpi(eu.getAnimalName, name) & c.hasPress).ExpName}));
    ai(iAnimal).nUnitsHasLickAndPress = nnz(strcmpi(eu.getAnimalName, name) & c.hasPress & c.hasLick);
    ai(iAnimal).nSessionsLickAndPress = length(unique({eu(strcmpi(eu.getAnimalName, name) & c.hasPress & c.hasLick).ExpName}));
    ai(iAnimal).nUnitsHasPos = nnz(strcmpi(eu.getAnimalName, name) & c.hasPress & c.hasLick & c.hasPos);
    ai(iAnimal).nSessionsHasPos = length(unique({eu(strcmpi(eu.getAnimalName, name) & c.hasPress & c.hasLick & c.hasPos).ExpName}));
end


%% Grouped table        
%                           Press   Press vs. Lick  2Tgt(spont)     4Tgt(Cued)  Spontaneous
% Group A: press only
% Group B: press vs lick
% Group C: Multi-target
% Group D: Spontaneous
clear aiGrps

N2 = horzcat(trajCombined2tgt.eta.N);
selUnits2 = all(N2 >= 4, 2);

N4 = arrayfun(@(eta) eta.N, trajCombined.eta, 'UniformOutput', false);
selUnits4 = N4{1, 1} >= 4 & N4{2, 1} >= 4 & N4{3, 1} >= 4 & N4{4, 3} >= 4;

SEL = {c.hasPress; c.hasPress & c.hasLick; []; cSpontaneous.hasPress};

aiGrps(2).animals = unique(eu(SEL{2}).getAnimalName);
aiGrps(1).animals = unique(eu(SEL{1}).getAnimalName);
aiGrps(1).animals = aiGrps(1).animals(~ismember(aiGrps(1).animals, aiGrps(2).animals));

aiGrps(1).eu = eu(SEL{1} & ismember(eu.getAnimalName, aiGrps(1).animals));
aiGrps(2).eu = eu(SEL{2} & ismember(eu.getAnimalName, aiGrps(2).animals));

aiGrps(3).eu = euReachDir2Tgt(selUnits2);
aiGrps(3).eu = aiGrps(3).eu(:);
aiGrps(3).animals = unique(aiGrps(3).eu.getAnimalName);

aiGrps(4).eu = euReachDir4Tgt(selUnits4);
aiGrps(4).eu = aiGrps(4).eu(:);
aiGrps(4).animals = unique(aiGrps(4).eu.getAnimalName);

aiGrps(5).animals = unique(euSpontaneous(SEL{4}).getAnimalName);
aiGrps(5).eu = euSpontaneous(:);

for i = 1:5
    aiGrps(i).eu = aiGrps(i).eu(:);
    aiGrps(i).sessions = unique({aiGrps(i).eu.ExpName});
    aiGrps(i).animals = aiGrps(i).animals(:);
    aiGrps(i).nByGender(1) = nnz(contains(aiGrps(i).animals, 'desmond'));
    aiGrps(i).nByGender(2) = nnz(contains(aiGrps(i).animals, 'daisy'));
end

for i = 1:2
    aiGrps(i).selfTimedReach = [length(aiGrps(i).eu), length(aiGrps(i).sessions)];
end

for i = 2
    aiGrps(i).selfTimedReachVsLick = [length(aiGrps(i).eu), length(aiGrps(i).sessions)];
end

for i = 3
    aiGrps(i).reach2tgt = [length(aiGrps(i).eu), length(aiGrps(i).sessions)];
end

for i = 4
    aiGrps(i).reach4tgt = [length(aiGrps(i).eu), length(aiGrps(i).sessions)];
end

for i = 5
    aiGrps(i).spontaneousReach = [length(aiGrps(i).eu), length(aiGrps(i).sessions)];
end
