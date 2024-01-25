# SimulatedEB1TipTracking
Script to simulated EB1 tip tracking along a dynamic microtubule. Used for https://www.biorxiv.org/content/10.1101/2022.06.07.495114v1.abstract. 


## About the project

This code is used to simulated EB1 tip tracking a dynamic microtubule. The user can adjust certain parameters in the top level script (labelled with EB1DynamicsBatching...). Of note, values is the matrix that contains all parameters. In the script, they are labelled as to what they control. Of note, the on/off rates for EB1 are first, followed by the number of steps (labelled as iterations), then tubulin characteristics such as concentration, hydrolysis rate, pibreak, etc. Of note, the last three variables (values 15-17) allow for slight adjustments to the script. If the user wants the simulation to remove binding at edge/lattice sites partway through the simulation, they can adjust values(15). Of note, at the moment, this removal is hard coded to occur at iteration 100,000 and 350,000 with it going back to its original value at 150,000 and 400,000 respectively. This can be adjusted in the rules script if wanted. Values(16) controls whether or not the simulation allows for splaying (which is achieved by simply saying all EB1 binding sites ahead of the highest lateral bond are considered and edge site). Of note, a value of 1 corresponds to that simple rule change (anything in front of the highest lateral bond is an edge site) while a value of 2 does the same rule as one and leads to the splaying being much greater (by lowering the rate at which lateral bonds form so that the highest lateral is further away from the end of the microtubule). Values(17) is used to determine if you want to record which GTP site an EB1 is bound to (edge versus lattice). Of note, the code will output four types of files. Two of them, proteinBindingAndRemovals... and protofilamentLengths... are both used downstream to generate the simulated images from this data which is the most useful output of the code. The proteinBins... is legacy code that I have left in at the moment; it is not very useful and is out of date. Finally, parameterTestingResultsSet will have some useful top level information about the simulation and will record the main parameters (as metadata so to speak). These file names should be changed with each run so that data is written over. 

The general way the code works is the top level script just feeds the user variables to the rules script. This is where the work is done. The way this underlying script works is for every iteration, it will calculate the time required for every type of event to occur (lateral bond form/break, EB1 bind/dissociate, tubulin dimer bind/dissociate). It will then perform whichever event takes the least amount of time to do and update the matrices that store the tubulin and EB1 information. Next, the time required for hydrolsyis to occur is calculated and if it is faster than the fastest other event, it will also occur at a tubulin dimer. Finally, this process is repeated many times with every cortime number of events being printed to the associated movie and recorded for the output excel files (by default, this occurs every 1000 steps in the simulation which ends up being around 1-2 seconds in the simulaton). 


### Built With
MATLAB by Mathworks

## Getting Started

To use the code, download the scripts. Next, make any necessary adjustments you want to parameters and adjust file names. Finally, running the script should start to generate simulaton data that will eventually be saved as excel files, two of which are used for downstream simulated image generation. Of note, there are no files that you need to input for this script. 


## Prerequisites

Ensure Matlab is operational on your device. 

## Installation

This script requires some add-ons including: 
- Image Processing Toolbox


## Usage

We have provided example output data for reference. The code should work out of the box and therefore the user can adjust parameters of interest and examine how the EB1 tip tracking is altered. 


Of note, this code was used for some publications from the Gardner lab including: 

https://www.biorxiv.org/content/10.1101/2022.06.07.495114v1.abstract 




## Contact

Samuel Gonzalez-(https://www.linkedin.com/in/samuel-gonzalez-081504163/) - samueljgonzalez@hotmail.com
