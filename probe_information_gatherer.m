filePaths = fileread('celFilePaths.txt');
filePaths = strsplit(filePaths);
filePaths = filePaths(1:(end-1));
noOfFiles = size(filePaths);
noOfFiles = noOfFiles(2);

filenames = strings(noOfFiles,1);
for i = 1:1:noOfFiles
    filename = filePaths{i};
    splitFilename = strsplit(filename(1:(end-4)), '/');
    filename = splitFilename{2};
    filename = strrep(filename, '-', '_');
    filenames(i) = filename;
end



CELStructure = struct(); 
for i = 1:1:noOfFiles
    CELStructure.(filenames{i}) = affyread(filePaths{i});
end



cdfStruct = affyread('Sp20b_M_v04.CDF'); 

SumOfProbeSetSizes = 0;
for i = 1:1:size(cdfStruct.ProbeSets)
    CurrentProbeSetSize = size(cdfStruct.ProbeSets(i).ProbePairs);
    CurrentProbeSetSize = CurrentProbeSetSize(1);
    SumOfProbeSetSizes = SumOfProbeSetSizes + CurrentProbeSetSize;
end



ProbeDataColNames = strings((2*noOfFiles)+6,1);
ProbeDataColNames(1) = 'ProbeSetName';
ProbeDataColNames(2) = 'Direction';
ProbeDataColNames(3) = 'PMX';
ProbeDataColNames(4) = 'PMY';
ProbeDataColNames(5) = 'MMX';
ProbeDataColNames(6) = 'MMY';


PMColNames = strcat('PM', filenames);
MMColNames = strcat('MM', filenames);

for i = 1:1:noOfFiles
    PMName = PMColNames(i);
    MMName = MMColNames(i);
    PMIndex = (i * 2) - 1 + 6;
    MMIndex = i * 2 + 6;
    ProbeDataColNames(PMIndex) = PMName;
    ProbeDataColNames(MMIndex) = MMName;
end



ProbeSetColNames = cdfStruct.ProbeSetColumnNames;
ProbeSets = cdfStruct.ProbeSets;

NoOfVariables = size(ProbeDataColNames);
NoOfVariables = NoOfVariables(1);

ProbeData = cell(SumOfProbeSetSizes,NoOfVariables);
counter = 1;
for i = 1:1:size(ProbeSets)
    disp(i);
    ProbeSet = ProbeSets(i);
    ProbeSetPairs = ProbeSet.ProbePairs;
    ProbeSetName = ProbeSet.Name;
    NoOfProbeSetPairs = size(ProbeSetPairs);
    NoOfProbeSetPairs = NoOfProbeSetPairs(1);
    for j = 1:1:NoOfProbeSetPairs
        ProbeRow = ProbeSetPairs(j,:);
        Direction = ProbeRow(2);
        PMX = ProbeRow(3);
        PMY = ProbeRow(4);
        MMX = ProbeRow(5);
        MMY = ProbeRow(6);
        ProbeData{counter,1} = ProbeSetName;
        ProbeData{counter,2} = Direction;
        ProbeData{counter,3} = PMX;
        ProbeData{counter,4} = PMY;
        ProbeData{counter,5} = MMX;
        ProbeData{counter,6} = MMY;
        for k = 1:1:noOfFiles
            filename = filenames{k};
            celStruct = CELStructure.(filename);
            PMDataColIndex = (k * 2) - 1 + 6;
            MMDataColIndex = (k * 2) + 6;
            PMProbeColName = ProbeDataColNames(PMDataColIndex);
            MMProbeColName = ProbeDataColNames(MMDataColIndex);
            PMCELRowIndex = find((celStruct.Probes(:,1) == PMX) & (celStruct.Probes(:,2) == PMY));
            PMIntensity = celStruct.Probes(PMCELRowIndex,3);
            MMCELRowIndex = find((celStruct.Probes(:,1) == MMX) & (celStruct.Probes(:,2) == MMY));
            MMIntensity = celStruct.Probes(MMCELRowIndex,3);
            ProbeData{counter,PMDataColIndex} = PMIntensity;
            ProbeData{counter,MMDataColIndex} = MMIntensity;
        end
        counter = counter + 1;
    end
end

ProbeDataTable = cell2table(ProbeData,'VariableNames',ProbeDataColNames);
writetable(ProbeDataTable,'ProbeDataTable.csv')




