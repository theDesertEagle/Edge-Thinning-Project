% Author: the.desert.eagle
% Program: Edge Detection, Thining and Gradients' Calculations

function a = imageProcessing(fileName, thresholdSobelGradMag, thresholdSobelGradDir, thresholdRobertsGradMag, thresholdRobertsGradDir)
%%%Default Parameters for rose.jpg
% fileName = 'rose.jpg';
% thresholdSobelGradMag = '80';
% thresholdSobelGradDir = '1.5';
% thresholdRobertsGradMag = '15';
% thresholdRobertsGradDir = '1.5';
%%%Default Parameters for flower.jpg
% fileName = 'flower.jpg';
% thresholdSobelGradMag = '150';
% thresholdSobelGradDir = '1.5';
% thresholdRobertsGradMag = '33';
% thresholdRobertsGradDir = '1.5';


%Reading image from graphics file
inputImage = imread(fileName);

%Display Input Image
figure(1)
imshow(inputImage); title('Input Image');

[initNumberOfRows, initNumberOfColumns] = size(inputImage);

%Image Padding Using Replicaiton - to account image borders during kernel convolution
paddedImage = [inputImage(:, 1) inputImage(:,:) inputImage(:, end)];
paddedImage = [paddedImage(1, :); paddedImage(:,:); paddedImage(end, :)];
paddedImage = double(paddedImage);

%Gradient Measurement Matrices Initialization
gX = zeros(initNumberOfRows, initNumberOfColumns, 'double');
gY = zeros(initNumberOfRows, initNumberOfColumns, 'double');
gX2 = zeros(initNumberOfRows, initNumberOfColumns, 'double');
gY2 = zeros(initNumberOfRows, initNumberOfColumns, 'double');
gXR = zeros(initNumberOfRows, initNumberOfColumns, 'double');
gYR = zeros(initNumberOfRows, initNumberOfColumns, 'double');
gX2R = zeros(initNumberOfRows, initNumberOfColumns, 'double');
gY2R = zeros(initNumberOfRows, initNumberOfColumns, 'double');
[padNumberOfRows, padNumberOfColumns] = size(paddedImage);

%Edge Detector Kernel Initialization - Sobel Operator
delta1 = [-1, 0, 1; -2, 0, 2; -1, 0, 1];
delta1 = double(delta1);
delta2 = [1, 2, 1; 0, 0, 0; -1, -2, -1];
delta2 = double(delta2);

%Edge Detector Kernel Initialization - Roberts Operator
deltaR1 = [0, 1; -1, 0];
deltaR1 = double(deltaR1);
deltaR2 = [1, 0; 0, -1];
deltaR2 = double(deltaR2);

%Application of Sobel kernels to all pixels of the image
for delta = 1:2
    for i = 2:padNumberOfRows - 1
        for j = 2:padNumberOfColumns - 1

            %Kernel Convolution for the pixel in consideration
            conResultOfPixel = 0;
            ii = 1;
            for p = i-1:i+1
                jj = 1;
                for q = j-1:j+1
                    if delta == 1
                        conResultOfPixel = conResultOfPixel + paddedImage(p, q) * delta1(ii, jj);
                    else
                        conResultOfPixel = conResultOfPixel + paddedImage(p, q) * delta2(ii, jj);
                    end
                    jj = jj + 1;
                end
                ii = ii + 1;
            end
            
            if delta == 1
                gX(i-1, j-1) = conResultOfPixel;
                gX2(i-1, j-1) = conResultOfPixel.^2;
            else
                gY(i-1, j-1) = conResultOfPixel;
                gY2(i-1, j-1) = conResultOfPixel.^2;
            end
        end
    end
end

%Application of Roberts kernels to all pixels of the image
for delta = 1:2
    for i = 2:padNumberOfRows-1
        for j = 2:padNumberOfColumns-1

            %Kernel Convolution for the pixel in consideration
            conResultOfPixel = 0;
            ii = 1;
            for p = i:i+1
                jj = 1;
                for q = j:j+1
                    if delta == 1
                        conResultOfPixel = conResultOfPixel + paddedImage(p, q) * deltaR1(ii, jj);
                    else
                        conResultOfPixel = conResultOfPixel + paddedImage(p, q) * deltaR2(ii, jj);
                    end
                    jj = jj + 1;
                end
                ii = ii + 1;
            end
            
            if delta == 1
                gXR(i-1, j-1) = conResultOfPixel;
                gX2R(i-1, j-1) = conResultOfPixel.^2;
            else
                gYR(i-1, j-1) = conResultOfPixel;
                gY2R(i-1, j-1) = conResultOfPixel.^2;
            end
        end
    end
end

%Gradient Magnitude Calculation - Sobel Operator
gX2gY2 = gX2 + gY2;
gMagnitude = sqrt(gX2gY2);
gMagView = zeros(initNumberOfRows, initNumberOfColumns);
%Normalizationo of Gradient Magnitude - Sobel Operator
sum = 0;
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        sum = sum + gMagnitude(i, j);
    end
end
meanVal = sum/(initNumberOfRows * initNumberOfColumns);
sum = 0;
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        sum = sum + (gMagnitude(i, j)-meanVal).^2;
    end
end
standardDeviationVal = sqrt(sum/((initNumberOfRows * initNumberOfColumns)-1));
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        gMagView(i, j) = (gMagnitude(i, j) - meanVal)/standardDeviationVal;
    end
end

%Display Gradient Magnitude = Sobel Operator
figure(2)
subplot(2,2,1)
imshow(gMagView); title('Gradient Magnitude - Sobel');


%Gradient Magnitude Calculation - Roberts Operator
gX2gY2R = gX2R + gY2R;
gMagnitudeR = sqrt(gX2gY2R);
gMagViewR = zeros(initNumberOfRows, initNumberOfColumns);
%Normalization of Gradient Magnitude - Roberts Operator
sum = 0;
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        sum = sum + gMagnitudeR(i, j);
    end
end
meanVal = sum/(initNumberOfRows * initNumberOfColumns);
sum = 0;
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        sum = sum + (gMagnitudeR(i, j)-meanVal).^2;
    end
end
standardDeviationVal = sqrt(sum/((initNumberOfRows * initNumberOfColumns)-1));
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        gMagViewR(i, j) = (gMagnitudeR(i, j) - meanVal)/standardDeviationVal;
    end
end

%Display Gradient Magnitude - Roberts Operator
figure(2)
subplot(2,2,2)
imshow(gMagViewR); title('Gradient Magnitude - Roberts');


%Gradient Direction Calculation - Sobel Operator
gDirection = atan2(gY, gX);
%Diplay Gradient Direction - Sobel Operator
figure(2)
subplot(2,2,3)
imshow(gDirection); title('Gradient Direction - Sobel');


%Gradient Direction Calculation - Roberts Operator
gDirectionR = atan2(gYR, gXR);
%Display Gradient Direction - Roberts Operator
figure(2)
subplot(2,2,4)
imshow(gDirectionR); title('Gradient Direction - Roberts');

%Thresholding the Sobel Operator Gradient Matrices
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        if gMagnitude(i, j) >= str2num(thresholdSobelGradMag)
            gMagnitude(i, j) = 255;
        else
            gMagnitude(i, j) = 0;
        end
    end
end
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        if gDirection(i, j) >= str2num(thresholdSobelGradDir)
            gDirection(i, j) = 255;
        else
            gDirection(i, j) = 0;
        end
    end
end

%Display the Gradient Matrices subjected to thresholding - Sobel Operator
figure(3)
subplot(2,2,1)
imshow(gMagnitude); title('Gradient Magnitude (Threshold) - Sobel');
figure(3)
subplot(2,2,3)
imshow(gDirection); title('Gradient Direction (Threshold) - Sobel');


%Thresholding the Roberts Operator Gradient Matrices
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        if gMagnitudeR(i, j) >= str2num(thresholdRobertsGradMag)
            gMagnitudeR(i, j) = 255;
        else
            gMagnitudeR(i, j) = 0;
        end
    end
end
for i = 1:initNumberOfRows
    for j = 1:initNumberOfColumns
        if gDirectionR(i, j) >= str2num(thresholdRobertsGradDir)
            gDirectionR(i, j) = 255;
        else
            gDirectionR(i, j) = 0;
        end
    end
end


%Display the Gradient Matrices subjected to thresholding - Roberts Operator
figure(3)
subplot(2,2,2)
imshow(gMagnitudeR); title('Gradient Magnitude (Threshold) - Roberts');
figure(3)
subplot(2,2,4)
imshow(gDirectionR); title('Gradient Direction (Threshold) - Roberts');


%Padding the Gradient Matrices for Line Thinning Preparation 
gMagnitude = [gMagnitude(:, 1) gMagnitude(:, :) gMagnitude(:, end)];
gMagnitude = [gMagnitude(1, :); gMagnitude(:, :); gMagnitude(end, :)];
gDirection = [gDirection(:, 1) gDirection(:, :) gDirection(:, end)];
gDirection = [gDirection(1, :); gDirection(:, :); gDirection(end, :)];
gMagnitudeR = [gMagnitudeR(:, 1) gMagnitudeR(:, :) gMagnitudeR(:, end)];
gMagnitudeR = [gMagnitudeR(1, :); gMagnitudeR(:, :); gMagnitudeR(end, :)];
gDirectionR = [gDirectionR(:, 1) gDirectionR(:, :) gDirectionR(:, end)];
gDirectionR = [gDirectionR(1, :); gDirectionR(:, :); gDirectionR(end, :)];

%Application of Zhan Seun Thinning Algorithm
%For matrix in each iteration, apply the algorithm iteratively
for iteration = 1:4
    if iteration == 1
        image = gMagnitude;
    elseif iteration == 2
        image = gMagnitudeR;
    elseif iteration == 3
        image = gDirection;
    elseif iteration == 4
        image = gDirectionR;
    end
    grid = zeros(9);
    %List of marked points for removal
    removeList1 = [-1];
    removeList2 = [-1];
    %Algorithm Loop
    while ~isempty(removeList1) || ~isempty(removeList2)
        removeList1 = [];
        removeList2 = [];
        %For each pass, carry out the following steps
        for pass = 1:2
            for i = 2:padNumberOfRows-1
                for j = 2:padNumberOfColumns-1
                    condition1 = 0;
                    condition2 = 0;
                    condition3 = 0;
                    condition4 = 0;
                    numberOfWhitePixels = 0;
                    numberOfBlackWhitePixelTransitions = 0;
                    if pass == 1
                        %Counting the number of White Pixels
                        if image(i, j) == 0
                            continue;
                        end
                        ii = 1;
                        for p = i-1:i+1
                            jj = 1;
                            for q = j-1:j+1
                                    grid(ii, jj) = image(p, q);
                                    jj = jj +1;
                                if p == i && q == j
                                    continue;
                                else
                                    if image(p, q) == 255
                                    numberOfWhitePixels = numberOfWhitePixels + 1;
                                    end
                                end
                            end
                            ii = ii + 1;
                        end
                        if numberOfWhitePixels >= 2 && numberOfWhitePixels <= 6
                            condition1 = 1;
                        end
                        %Black to White Transition Checking 
                        if grid(1) == 255 && grid(4) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(2) == 255 && grid(1) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(3) == 255 && grid(2) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(4) == 255 && grid(7) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(6) == 255 && grid(3) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(7) == 255 && grid(8) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(8) == 255 && grid(9) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(9) == 255 && grid(6) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if numberOfBlackWhitePixelTransitions == 1
                            condition2 = 1;
                        end
                        %Atleast One Black Neighbor from Three Specfic Neighbors
                        %Checking 
                        if grid(2) == 0 || grid(4) == 0 || grid(6) == 0
                            condition3 = 1;
                        end
                        if grid(6) == 0 || grid(4) == 0 || grid(8) == 0
                            condition4 = 1;
                        end
                        %Mark Pixel for removal if all conditions satisfied                        
                        if condition1 == 1 && condition2 == 1 && condition3 == 1 && condition4 == 1
                            removeList1 = [removeList1 [i j]];
                        end
                    else
                        %Counting the number of White Pixels
                        if image(i, j) == 0 
                            continue;
                        end                    
                        ii = 1;
                        for p = i-1:i+1
                            jj = 1;
                            for q = j-1:j+1
                                    grid(ii, jj) = image(p, q);
                                    jj = jj +1;
                                if p == i && q == j
                                    continue;
                                else
                                    if image(p, q) == 255
                                    numberOfWhitePixels = numberOfWhitePixels + 1;
                                    end
                                end
                            end
                            ii = ii + 1;
                        end
                        if numberOfWhitePixels >= 2 && numberOfWhitePixels <= 6
                            condition1 = 1;
                        end
                        %Black to White Transition Checking
                        if grid(1) == 255 && grid(4) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(2) == 255 && grid(1) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(3) == 255 && grid(2) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(4) == 255 && grid(7) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(6) == 255 && grid(3) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(7) == 255 && grid(8) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(8) == 255 && grid(9) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if grid(9) == 255 && grid(6) == 0
                            numberOfBlackWhitePixelTransitions = numberOfBlackWhitePixelTransitions + 1;
                        end
                        if numberOfBlackWhitePixelTransitions == 1
                            condition2 = 1;
                        end
                        %Atleast One Black Neighbour from Three Specfic
                        %Neighbors Checking
                        if grid(2) == 0 || grid(4) == 0 || grid(8) == 0
                            condition3 = 1;
                        end
                        if grid(2) == 0 || grid(6) == 0 || grid(8) == 0
                            condition4 = 1;
                        end
                        if condition1 == 1 && condition2 == 1 && condition3 == 1 && condition4 == 1
                            removeList2 = [removeList2 [i j]];
                        end
                    end
                end
            end
            %Pixel Removal which are marked in the list
            if pass == 1 && ~isempty(removeList1)
                [s1,s2] = size(removeList1);
                for k = 1:s1
                    a = removeList1(k);
                    b = removeList1(k+1);
                    image(a, b) = 0;
                    k = k + 1;
                end
            elseif pass == 2 && ~isempty(removeList2)
                [s1, s2] = size(removeList2);
                for k = 1:s1
                    a = removeList2(k);
                    b = removeList2(k+1);
                    image(a, b) = 0;
                    k = k + 1;
                end
            end
        end
    end
    %Display relevant diagram on each iteration
    if iteration == 1
        figure(4)
        subplot(2,2,1)
        imshow(image);
        title('Zhang Suen Thinning On Gradient Magnitude - Sobel');
    elseif iteration == 3
        figure(4)
        subplot(2,2,3)
        imshow(image);
        title('Zhang Suen Thinning On Gradient Direction - Sobel');
    elseif iteration == 2
        figure(4)
        subplot(2,2,2)
        imshow(image);
        title('Zhang Suen Thinning On Gradient Magnitude - Roberts');
    elseif iteration == 4
        figure(4)
        subplot(2,2,4)
        imshow(image);
        title('Zhang Suen Thinning On Gradient Direction - Roberts');
    end
end