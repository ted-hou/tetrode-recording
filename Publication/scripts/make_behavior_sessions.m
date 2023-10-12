%% Make behavior session objects and save to server (only run once)
dirs = dir('C:\SERVER'); 
dirs = cellfun(@(x, y) sprintf('%s\\%s', x, y), {dirs.folder}, {dirs.name}, 'UniformOutput', false); 
dirs = dirs(isfolder(dirs) & contains(dirs, {'daisy', 'desmond'}) & ~contains(dirs, 'unfinished'));

for iAnimal = 1:length(dirs)
    animalName = strsplit(dirs{iAnimal}, '\');
    animalName = animalName{end};
    f = dir(sprintf('C:\\SERVER\\%s\\%s_*\\%s_*.mat', animalName, animalName, animalName));
    startIndex = regexpi({f.name}, '(daisy|desmond).*[0-9]{8}\.mat$');
    isMatch = ~cellfun(@isempty, startIndex);
    f = f(isMatch);
    files = cellfun(@(x, y) sprintf('%s\\%s', x, y), {f.folder}, {f.name}, 'UniformOutput', false);
    bs{iAnimal} = BehaviorSession(files);
end

mkdir('C:\SERVER_PRIVATE\data')
save('C:\SERVER_PRIVATE\data\bs.mat', 'bs');