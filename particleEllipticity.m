function particleEllipticity(dataX, dataY,fileName,output,saveResults)
    % saveResults = if set to 1, will save results to Results folder
    % output = folder specifying where output should be saved
    % fileName = name of file being analysed (for listing in the results file)
    % data X and Y : X and Y axis data (calibrated waveforms)
    
    % Define Output folders
    if ispc
        pathResults = [output '\Results\Ellipticity.csv'];
    else
        pathResults = [output '/Results/Ellipticity.csv'];
    end

    n = pow2(nextpow2(min(length(dataX),length(dataY)))); % calculate next power of 2, for padding fft length
    
    % Maximum amplitude value (uPa, nm/s, or um/s^2)
    maxHist = max([dataX;dataY]);

    % Define bin boudaries for histogram (200x200 bins)
    edges = {-maxHist:maxHist*2/99:maxHist; -maxHist:maxHist*2/99:maxHist};

    % Plot histogram
    histogram = hist3([dataX dataY],'Edges',edges);
    Max = max(max(histogram));
    thresh = find(histogram>(Max/4)); % thresh = Index of all values over 25% of max amplitude
    y = mod(thresh,100) + 1;
    x = floor(thresh/100) + 1;
    cvh = convhull(x,y);
    vectors = [x(cvh) y(cvh)];

    %longest distance between 2 convex hull points
    MajorAxis = 0;
    for u=1:length(vectors)
        for o=1:length(vectors)
            tempHyp = sqrt((vectors(u,1)-vectors(o,1))^2 + (vectors(u,2)-vectors(o,2))^2);
            if  tempHyp > MajorAxis;
                MajorAxis = tempHyp;
                I1 = u;
                I2 = o;
            end
        end
    end

    %Get slope of major axis
    slope = (vectors(I1,2) - vectors(I2,2)) / (vectors(I1,1) - vectors(I2,1));
    %Rotate slope 90 degrees
    slopex = -(vectors(I1,1) - vectors(I2,1)) / (vectors(I1,2) - vectors(I2,2));
    % define end points of line
    % y = mx+b
    intercept = 50 - slopex*50;
    point1 = slopex*100 + intercept;
    point2 =  intercept;

    intercept2 = 50 - slope*50;
    point1x = slope*100 + intercept2;
    point2x =  intercept2;

    [xi,yi] = polyxpoly([0 100],[point2 point1],vectors(:,1),vectors(:,2));

    MinorAxis = sqrt((xi(1) - xi(2))^2 + (yi(1) - yi(2))^2);

    ellipticity = atan(MinorAxis/MajorAxis);

    bearing = atan2(slope,1);

    %% Plot proofing graph                                  
    figure(25);
    imagesc(hist3([dataX dataY],'Edges',edges));
    hold on;
    plot(x(cvh),y(cvh),'r');
    plot([0 100],[point2 point1],'g');
    plot([0 100],[point2x point1x],'k');
    set(gca,'YDir','normal');
    axis equal;
    xlabel('X axis');
    ylabel('Y axis');
    title('lines reperesent the Major axis, and the Adjacent axis');

    disp(['Ellipticity (degrees): ' num2str(ellipticity*(180/pi))]);
    disp(['Bearing (degrees): ' num2str(bearing*(180/pi))]);
    %%%%%%%%%%%%%%%%%%%%%%%%%

    if saveResults == 1;
        fw = fopen(pathResults,'a');
        fprintf(fw, '%s',[fileName ',' num2str(ellipticity*(180/pi)) ',' num2str(bearing*(180/pi))]);
        fprintf(fw, '\n');
        fclose(fw);
    end
