# Image Analysis Hub

This MATLAB was developed for use by Dr. Paula Barrett's lab in the UVa Department of Pharmacology. This app allows users to quickly analyze calcium oscillatory events and view fluorescent traces.

## System requirements
This software will run on any computer that has MATLAB version 2018a or later installed. It is compatible with Windows, Mac, and Linux operating systems. This app has been tested on Windows, Mac, and Linux operating systems and with MATLAB versions 2018a, 2018b, and 2019a.

## Installation Instructions
No true installation is required. Simply download and extract the ImageAnalysisHub folder on a machine that has MATLAB version 2018a or later. Download and extraction should take approximately 5 minutes.

# Demo Instructions (Expected Run Time: 5-10 minutes)

## Getting Started

1.	Navigate to the ImageAnalysisHub folder on your computer.

2.	Launch the app by double-clicking ImageAnalysisHub.mlapp.

## Load Event Data

### Choose Data
1.	Click the “Load Event Data” button to open the data selection menu.

2.	Click the “Select” button in the “Event Detection Algorithm” box. (minEASE should be selected by default)

3.	Click “x3nMAng_data30” in the “Directory Contents” box to select that dataset and click “Add Selected”.

4.	Click the “Next” button.

### Parameters
1.	Enter “60” in the “Bin Time” field. This will divide the data into 60 second bins.

2.	Enter “90” in the “Start Time” field. In this experiment the Angiotensin II was added at 90 seconds.

3.	Enter “600” in the “End Time” field.

4.	Click “Load Data”. This will create a MAT-file containing your chosen data and parameters. The MAT-file will be saved to the “1_Loaded_Datasets” folder and will have the timestamp in its name (e.g. eventData_dd_MMM_YYYY_HH.mm.mat).

5.	Close the “Load Event Data” menu.

## Event Analysis

### Parameters

1.	Click the “Event Analysis” button to open the Event Analysis menu.

2.	Click the “Browse” button and then click the MAT-file you created in the “Load Event Data” section.

3.	Click “Open”. (This sometimes brings the main MATLAB window to the forefront. You may need to click on the “Event Analysis” window to return to that menu.)

4.	Click the “Load” button. This will show you the data and parameters selected in the “Load Event Data” section.

5.	Click the “Next” button.

### Threshold Algorithm

1.	Select “Gaussian + Exponential (Intersection)”.

2.	Click the “Next” button.

### Measurements
1.	Select the measurements you would like in the output Excel file. (By default, everything is selected).

2.	Click the “Run Analysis” button. This will create an Excel file with the dataset and timestamp in its name. The Excel file will be written to the “2_Event_Analysis” folder.

3.	Close the “Event Analysis” menu.

## Burst Analysis

### Parameters

1.	Click the “Event Analysis” button to open the Event Analysis menu.

2.	Click the “Browse” button and then click the MAT-file you created in the “Load Event Data” section.

3.	Click “Open”. (This sometimes brings the main MATLAB window to the forefront. You may need to click on the “Event Analysis” window to return to that menu.)

4.	Click the “Load” button. This will show you the data and parameters selected in the “Load Event Data” section.

5.	Click the “Next” button.

### Threshold Algorithm

1.	Select “Gaussian + Exponential (Intersection)”.

2.	Click the “Next” button.

### Measurements
1.	Select the measurements you would like in the output Excel file. (By default, everything is selected).
2.	Click the “Run Analysis” button. This will create an Excel file with the dataset and timestamp in its name. The Excel file will be written to the “3_Burst_Analysis” folder.
3.	Close the “Burst Analysis” menu.

