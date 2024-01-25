%%% Copyright Gardner lab 2023.
%%%This code is the batching script to run the simualtion for EB1 tip
%%%tracking. Essentially, in this script, the user can manually decide the
%%%on/off rates for EB1, the number of iterations per simulation run, the
%%%number of simulations to run, the tubulin concentration and hydrolysis
%%%rate, the tubulin on rate, the pibreak, whether or not to generate a
%%%movie output (of note, the movie output is at the dimer level and not a
%%%simulated image from an experimental TIRF assay), and the acceptable
%%%taper length. The user can also choose the file name. The script will
%%%first generate the output excel files and populate them with headers.
%%%Then, it will run the underlying simulation script, take its outputs,
%%%and store them into excel files. Of note, the relevant excel files would
%%%be the parameter testing results file, the protein binding and removals
%%%file, and the protofilament lenghts file (the later two are used later
%%%to generate the simulated experimental images while the first one will
%%%contain some metadata and other legacy characteristics). This script
%%%will also keep track of how long it takes to run each loop so that the
%%%user can plan accordingly (in the past, depending on the iteration
%%%number and how long the MT gets (based on tubulin concentration and
%%%therefore growth rate), the simulation can take a long time to run as it
%%%takes longer to run as the MT gets longer.

close all;
clear all;
rng('shuffle');%create random number generator seed for catastrophe...
%Initialize values used for the simulation. 
 values(1)=0.000023;%%EB1 on rate at lattice sites. nM-1 site-1 s-1
 values(2)=0.0016;%% EB1 on rate at edge sites nM-1 site-1 s-1
 values(3)=25;%%EB1 off rate at GDP Edge sites s-1
 values(4)=1.7;%%EB1 off rate at GDP Lattice sites s-1
 values(5)=2.9;%%EB1 off rate at GTP Edge sites s-1
 values(6)=0.29;%%EB1 off rate at GTP Lattice Sites, s-1 
 values(7)=160000;%iterations number, essentially the number of steps in the simulation. 
 values(8)=12;%tubulin concentration uM
 values(9)=0.55/13;%%Tubulin Hydrolyis rate per protofilament (0.55 is per MT, s-1)
 values(10)=0.65*13;%%Tubulin on rate per microtubule (0.65 is per PF, uM-1 s-1 pf-1)
 values(11)=10;%pibreak%%1000 so its only really growing...no large tubulin removal occurs at this pibreak "orignally" 5.4...For all for paper, use 10...
 values(12)=0.200;%EB1 concentration normally 0.250, changed to 200 nM since thats what I use in the experiments...
 values(13)=1;%whether or not a movie is generated. 0 will not generate a movie while 1 will generate a movie (slows code but easier for seeing effect real time). 
 values(14)=75;%acceptable taper lenght in Dimers. 
 values(15)=0;%If you want edge binding or lattice binding to be removed partway through the simulation. 0 means no, 1 for edge, 2 for lattice. 
 values(16)=1;%whether or not there is splaying. 0 is no, 1 only having in front of highest lateral bond as splaying, 2 is 1 plus it makes the lateral bonds further back (ie making it extreme splaying).  
 values(17)=1;%whether or not you want to record if EB1 is bound at GTP is at an edge or a lattice
 number=num2str(6);%File number used to quickly change the name of each file (for when you are running many files worth). 
 filename=strcat("parameterTestingResultsSet",number,"Condition1.xlsx");%establish file name for parameter results, stores parameters and some top level measurements which may be out of date. 
 filename2=strcat("proteinBinsResultsSet",number,"Condition1.xlsx");%establish file name for proteins binned results, legacy code; not useful. 
 filename3=strcat("proteinBindingAndRemovals",number,"Condition1.xlsx");%Establish file for protein tracking and chosen events lists
 filename4=strcat("protofilamentLengths",number,"Condition1.xlsx");%makes files to hold PF  lengths
if values(17)==1
    filename5=strcat("proteinBindingAndRemovalsWithEvL",number,"Condition1.xlsx");%Makes a file that keeps protein position AND where it originally bound. 1 is edge, 2 is lat...
end
%Decide how many iterations you want the simulation to be run. Time per
%simulation varies with iteration number. Of note, the code gets
%considerably slower as the microtubule gets longer so keep that in mind
%when deciding the iteration number. 
for i=1
%keeping track of time
     tic
 index=1;%initializing index for file names.
 i %%printing out which loop is being run for book keeping in command window. 
                        %generate file titles 
                        resultsTitles=["parameterDetails","length","taper length","total proteins","# proteins in taper","# proteins 30 dimers before taper","#GTP edge proteins","#GTP lattice proteins","%GTP proteins at edge","#GDP edge proteins","#GDP lattice proteins","%GDP proteins at edge","absolute time","decay factor","fittedMTLength"];
                        MTbinValues=[0:10:400];
                        %store current parameter sets for later use (as
                        %metadata so to speak). 
                        parameters=strcat("tub conc ",num2str(values(8)),"um, khydr=",num2str(values(9)),", kprotlatt=",num2str(values(1)),", kprotEdge=,",num2str(values(2)),", kpoffedgeGDP=",num2str(values(3)),", kpofflattgdp=",num2str(values(4)),", kpoffedgegtp=",num2str(values(5)),", kpofflattGTP=",num2str(values(6)),", kplusMT=",num2str(values(10)),",nIteraations=",num2str(values(7)),", imageFrequency=every1000Iterations, FunctionV10","pibreak=",num2str(values(11)), "EB1Concentration=",num2str(values(12)),"TaperLengthForRule=",num2str(values(14))," Removing Binding=",num2str(values(15)), " Splaying=",num2str(values(16)), " GTPEvL=",num2str(values(17)));%change later...
                        %titles for protofilament length file. 
                        protofilamentNumbers=["PF1","PF2","PF3","PF4","PF5","PF6","PF7","PF8","PF9","PF10","PF11","PF12","PF13"];
                        %write titles and initial data (like parameter
                        %sets) to each file ahead of time. Can take some
                        %time to note. 
                             writematrix(resultsTitles,filename,'Sheet',i, 'AutoFitWidth',false,'UseExcel',true);
                             writematrix(parameters,filename,'Sheet',i,'Range','A2', 'AutoFitWidth',false,'UseExcel',true);
                             writematrix(MTbinValues,filename2,'Sheet',i,'Range','B1', 'AutoFitWidth',false,'UseExcel',true);
                             writematrix("MT position bin",filename2,'Sheet',i,'Range','A1', 'AutoFitWidth',false,'UseExcel',true);
                             writematrix("number of proteins",filename2,'Sheet',i,'Range','A2', 'AutoFitWidth',false,'UseExcel',true);
                             writematrix(parameters,filename2,'Sheet',i,'Range','A3', 'AutoFitWidth',false,'UseExcel',true);
                             writematrix(protofilamentNumbers,filename4,'Sheet',i, 'AutoFitWidth',false,'UseExcel',true);
                            %Run the underlying simulation code 
                            [binVals,outputdatas,proteinTracking,chosenEventList,proteinPlaces,proteinPlacesWithEvL,protoLengths,fittedMTLength]=MTV15_newRules_RESET_transfer_to_samSJGUpdateCommentedBetter(values);%run the function 
                            %write the current output from this simulation
                            %run to file. 
                            writematrix(outputdatas,filename,'Sheet',i,'Range','B2', 'AutoFitWidth',false,'UseExcel',true);
                            writematrix(fittedMTLength,filename,'Sheet',i,'Range','O2', 'AutoFitWidth',false,'UseExcel',true);
                            writematrix(binVals,filename2,'Sheet',i,'Range','B2', 'AutoFitWidth',false,'UseExcel',true);
                            proteinPlaces=proteinPlaces';
                            overallProteinPlaces=NaN(20000,(values(7)/1000)*2);%this is because do an image every 1000 iterations and want x and y values for each protein at each time point. 
                            %work to setup the output file that has the
                            %protein (EB1) positions in dimers. 
                            for k=1:size(proteinPlaces,2)
                            proteinPlacesMatrix=cell2mat(proteinPlaces(k));
                            indexNumber=k*2-1;
                            overallProteinPlaces(1:size(proteinPlacesMatrix,1),k*2-1:k*2)=proteinPlacesMatrix;
                            indexNumber=xlscol(indexNumber);
                            %indexNumberString=num2str(indexNumber);
                            startingSpot=strcat(indexNumber,"1");
                            
                            end
                            if values(17)==1
                                for k=1:size(proteinPlaces,2)
                                proteinPlacesMatrixEvL=cell2mat(proteinPlacesWithEvL(k));
                                indexNumber=(k-1)*3+1;
                                overallProteinPlacesEvL(1:size(proteinPlacesMatrixEvL,1),(k-1)*3+1:(k-1)*3+3)=proteinPlacesMatrixEvL;
                                indexNumber=xlscol(indexNumber);
                                startingSpot=strcat(indexNumber,"1");

                                end
                            end
                            %Okay, all the extra zeros in protein place
                            %matrix cause an issue when making the
                            %simulated videos, so need to delete them. To
                            %do that, do: This will essentially go column
                            %by column and once you hit a zero, you set
                            %that row value to the end of the matrix as
                            %emtpy "" so that it won't affect the simulated
                            %MT growth at the end. 
                            for increment2=1:size(overallProteinPlaces,2)
                                for increment1=1:size(overallProteinPlaces,1)
                                
                                    if overallProteinPlaces(increment1,increment2)==0
                                        overallProteinPlaces(increment1:size(overallProteinPlaces,1))="";
                                        break
                                    end
                                end
                            end
                            writematrix(overallProteinPlaces,filename3,'Sheet',i,'AutoFitWidth',false,'UseExcel',true);
                            writematrix(protoLengths,filename4,'Sheet',i,'Range','A2', 'AutoFitWidth',false,'UseExcel',true);
                            if values(17)==1
                                for increment4=1:size(overallProteinPlacesEvL,2)
                                    for increment3=1:size(overallProteinPlacesEvL,1)
                                        if overallProteinPlacesEvL(increment3,increment4)==0
                                            overallProteinPlacesEvL(increment3:size(overallProteinPlacesEvL,1),increment4)="";
                                            break
                                        end
                                    end
                                end
                                writematrix(overallProteinPlacesEvL,filename5,'Sheet',i,'AutoFitWidth',false,'UseExcel',true);
                            end
                            binVals=[];
                            outputdatas=[];
                            proteinPositions=[];
                            events=[];
%%If you want to play with parameters in a systematic way, using a for loop
%%here where you increment every so many loops will allow for that. For
%%example, you can use the mod operator to every 5 loops say double the
%%value for the lattice on rate. This was helpful in initial parameter
%%testing. 
 
                        index=index+1;

toc %end tracking time per run (including the portion required to write as well since that can take some time). 
end

