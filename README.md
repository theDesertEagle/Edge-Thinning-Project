# Outlining Potential Objects of Interest Using Zhan Suen Thinning   
A Traditional Computer-Vision program to compare the results of Sobel and Roberts Edge Detectors by computing their gradient magnitudes and gradient directions respectively. The edges are then manually thresholded, segmented and thinned using the Zhan Suen Thinning algorithm, to estimate an outline of the object of interest.     

## Instructions
Run the `imageProcessing.m` script in the MATLAB console terminal, using the following command: </br> </br>
`path$: imageProcessing(fileName, thresholdSobelGradMag, thresholdSobelGradDir, thresholdRobertsGradMag, thresholdRobertsGradDir)` </br> </br>
where, </br>
- **fileName**                 = String depicting the graphics file which is to be used for image processing </br>
- **thresholdSobelGradMag**    = Thresholding value that is to be applied on the Sobel Gradient-Magnitude Matrix
- **thresholdSobelGradDir**    = Thresholding value that is to be applied on the Sobel Gradient-Direction Matrix
- **thresholdRobertsGradMag**  = Thresholding value that is to be applied on the Roberts Gradient-Magnitude Matrix
- **thresholdRobertsGradDir**  = Thresholding value that is to be applied on the Roberts Gradient-Direction Matrix

The individual thresholds, mentioned above, are user-driven manual intensity-values for improved object-segmentation from noise.

The recommended set of threshold values for the given set of sample images are given below:
- Script Parameters for rose.jpg
  - thresholdSobelGradMag = 80
  - thresholdSobelGradDir = 1.5
  - thresholdRobertsGradMag = 15
  - thresholdRobertsGradDir = 1.5
- Script Parameters for flower.jpg
  - thresholdSobelGradMag = 150
  - thresholdSobelGradDir = 1.5
  - thresholdRobertsGradMag = 33
  - thresholdRobertsGradDir = 1.5

## Observations
It was noticed that Sobel edge-detection was slightly slower than Roberts, although the smoother edges obtained by Sobel indicated that Roberts was more sensitive to noise. Moreover, thresholding and thinning gradient magnitudes produced meaningful object-outline results, unlike in the case of gradient directions which produced uninterpretable visual results.
### Sample Image: flower.jpg
#### Input Image
![flower image](/flower.jpg)
#### Edge Detection Results
![flower edge detection image](/docImg/flEdgeDetect.PNG)
#### Thrershold Results
![flower edge thresholded image](/docImg/flThresh.PNG)
#### Edge Thinning Results
![flower edge thinning image](/docImg/flEdgeThin.PNG)
### Sample Image: rose.jpg
#### Input Image
![rose image](/rose.jpg)
#### Edge Detection Results
![rose edge detection image](/docImg/roEdgeDetect.PNG)
#### Thrershold Results
![rose edge thresholded image](/docImg/roThresh.PNG)
#### Edge Thinning Results
![rose edge thinning image](/docImg/roEdgeThin.PNG)

## Known Bugs/Issues
The script may run slow for considerably larger image sizes due to the underlying thinning process. This might in turn delay the output of already calculated-results for tasks prior to edge-thinning. It is aimed to improve the speed of the edge-thinning process in future releases.

## @the.desert.eagle
