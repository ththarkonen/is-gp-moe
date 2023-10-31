function [ highDensityMatrix, medianEstimate] = computeHighDensityArea( ismoeObject, nX, nY)

    resultCDF = cumsum( ismoeObject.yDataGrid );
    highDensityMatrix = NaN( nY, nX);
    medianEstimate = zeros( nX, 1);

    for ii = 1:nX

        C_ii = result.yDataGrid(:,ii);
        resultCDF_ii = resultCDF(:,ii);

        [cOpt, ~] = fminsearch( @(c) ( sum(C_ii( C_ii <= c )) - 0.10 )^2, 0.000);
        
        nonUniques = resultCDF_ii > 0.99 | resultCDF_ii < 1e-3;
        nNonUniques = sum( nonUniques );

        highDensityMatrix( C_ii > cOpt, ii) = 1;

        resultCDF_ii( nonUniques ) =  resultCDF_ii( nonUniques ) + 0.001 * randn( nNonUniques, 1);
        resultCDF_ii = resultCDF_ii + 0.000000001 * rand( nY, 1);
        medianEstimate(ii) = interp1( resultCDF_ii, yStar, 0.50);
    end
end