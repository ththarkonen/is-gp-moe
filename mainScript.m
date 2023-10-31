
clear
close all
    
includeFolders = genpath('include');
addpath( includeFolders );

fileNames = ["stationary", "nonstationary"];
alphas = [ 0.1, 1, 3.5];

nFiles = length( fileNames );
nAlphas = length( alphas );

for fileName = fileNames
    for concentrationParameter = alphas

        filePath = "./data/" + fileName + ".mat";
        data = load( filePath );
        data = preprocess( data );

        [ data, settings] = preprocess( data );
        settings.prior.psi.concentrationParameter = concentrationParameter;

        ismoeObject = isgpmoe( data, settings);
        
        timeStamp = datetime('now', 'Format', 'yyyy-mm-dd-HH-MM-SS');
        timeStamp = string( timeStamp );
        
        resultPath = "./results/" + timeStamp + "_" + fileName;
        resultPath = resultPath  + "_ismoe_" + num2str( concentrationParameter ) + ".mat"
        
        save( resultPath, 'ismoeObject')
    end
end