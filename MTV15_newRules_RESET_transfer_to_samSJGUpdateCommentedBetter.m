%%
%%%The overall goal of this function is to run the simulation of EB1 tip
%%%tracking. In general, the data from the batching script will be loaded
%%%into variables in the file. Then, there are two matrices that keep track
%%%of all the information for EB1 and the tubulin dimers labelled as
%%%"EB1BindingPositions" and "tubAProps". TubAProps has 5 layers (legacy
%%%code) but layer 1 is important (it keeps track if a tubulin dimer is GTP
%%%(0) or GDP (1). As for EB1BindingPositions, this keeps track of all EB1
%%%binding sites. The first layer keeps track of if the binding site is an
%%%edge vs lattice site. The second layer keeps track of the hydrolysis
%%%state of the pocket (that depends on the tubulin dimers on the minus end
%%%of the pocket), and layer 3 includes whether or not an EB1 is bound. In
%%%short, the way the simulation works is that it will calculate the rate
%%%for a multitude of things occuring (lateral bond break/form, tubulin
%%%dimer add/remove, EB1 bind/remove). Then, whichever event occurs fastest
%%%(by taking the least amount of time) occurs. Next, the time needed for
%%%hydrolysis to occur is calculate and if it is faster than the fastest
%%%other event, hydrolysis will also occur at a tubulin dimer. Then, all
%%%positions are updated in both TubAProps and EB1BindingPositions and the
%%%process is repeated. Every cortime (number of steps that you want to
%%%take snapshots of), the simulation will record data to a file that can
%%%then be output to the batching file and subsequently saved to file. Of
%%%note, the speed of the code can likely be improved (the main crux at the
%%%moment is the calculation of the on-rate which depends on determining
%%%the total number of sites for each type, edge vs lat). The current
%%%iteration of the code essentially sums the number of edge vs lattice
%%%sites and then can calculate the time required for EB1 to bind to either
%%%site to occur, but because the microtubule lenght increases over time,
%%%the code slows down over time. There is likely a cleaner way to keep
%%%track of the total number of sites that are edge or lattice, but at the
%%%time of writing I wasn't able to figure one out and therefore the code
%%%is as is.
function[binValues,outputData,proteinTracking,chosenEventList,proteinPlaces,proteinPlacesWithEvL,protosLengths,fittedMTLength]=MTV15_newRules_RESET_transfer_to_samSJGUpdateCommentedBetter(values)%kProtLattice,kProtEdge,,kPOffEdgeGDP,kPOffEdgeGTP,kPOffLatticeGDP,kPOffLatticeGTP,videoOutputName, (order of indices for values)
tic

%for debugging
%likely redundant here with the paths, but whatever.
addpath('R:\cbs_lab_klei0091\Sam\Taylor code for MT dynamics with tpx2\Microtubule\Include')
addpath('Z:\Groups\LAB-klei0091\Sam\Taylor code for MT dynamics with tpx2\MicroTubule\Include')
addpath('Y:\cbs_lab_klei0091\TaylorReid\MatlabCode\MicroTubule\Include')

%% Set Simulation Perameters
if(exist('values','var'))
    kProtLattice=values(1); % on rate for EB1 at lattice sites.
    kProtEdge=values(2);% on rate for EB1 at edge sites.
    kPOffEdgeGDP=values(3);% off rate for EB1 at GDP edge sites.
    kPOffLatticeGDP=values(4); % off rate for EB1 at GDP lattice sites.
    kPOffEdgeGTP=values(5);% off rate for EB1 at GTP edge sites.
    kPOffLatticeGTP=values(6); % off rate for EB1 at GTP lattice sites.
    nIterations=values(7); %number of steps for the simulation to run
    tubGtpC=values(8); %tubulin concentration in uM
    kHyd=values(9); %hydrolysis rate per protofilament
    kPlusMT=values(10); %tubulin on rate
    Pibreak=values(11); %correction factor for the lateral bond breakage rate if all neighbor PFs have lateral bonds. Essentially, controls how tapered the MT is with a higher value leading to reduced taper and a lower value leading to a higher taper.
    LateralBreakingFoldFactor=1;%Will divide the off rate for all lateral bonds breaking by this factor, initialized and always set as 1. Reset if conditions are met.

    protConcUM=values(12);%EB1 concentration
    movie_on=values(13);%0 means dont, 1 means yes...
    filN=13;%number of protofilaments for the microtubule. Assumed to be 13 throughout.
    LongestTaperLength=values(14);%#of dimers that the longest PF can be greater than smallest two PFs before lateral bonds break quickly. For simulations was 75 but was increased to allow for larger microtuble splaying for split comets.
    removeEdgeOrLatBinding=values(15);%if you want edge or lattice binding removed partway through the simulation.
    splayingOccurs=values(16);%if you want to have splaying or not.
    GTPEvL=values(17);%If we want to track EB1 bound to GTP edge versus GTP lattice. 
    % videoOutputName=values(7);
    videoOutputName='mymovie.avi';%for now...can change later...
else
    videoOutputName='mymovie.avi';%Give the name for the output video
    kProtLattice=0.000023;%%EB1 on rate at lattice sites. nM-1 site-1 s-1
    kProtEdge=0.0016;%% EB1 on rate at edge sites nM-1 site-1 s-1
    kPOffEdgeGDP=25;%%EB1 off rate at GDP Edge sites s-1
    kPOffLatticeGDP=1.7;%%EB1 off rate at GDP Lattice sites s-1
    kPOffEdgeGTP=2.9;%%EB1 off rate at GTP Edge sites s-1
    kPOffLatticeGTP=0.29;%%EB1 off rate at GTP Lattice Sites, s-1
    nIterations=450000;%%Iteration number
    tubGtpC=12;%Tubulin concentration in uM
    filN=13;%Filament number
    kHyd=0.55/filN;%Hydrolysis rate s-1 per PF
    kPlusMT=0.65*filN;%Tubulin on rate per microtubule uM-1 s-1
    protConcUM=0.200;%EB1 concentration in uM
    LongestTaperLength=75;%Longest allowed taper length before tubulin off rate increased.
    removeEdgeOrLatBinding=0;%don't removed edge or lattice binding.
    splayingOccurs=0;%if you want to have splaying or not.
    GTPEvL=0;%If we want to track EB1 bound to GTP edge versus GTP lattice. 

    movie_on=1;%whether or not a movie is generated.
    Pibreak = 10; % correction for when both lat neighboring bonds exist.
    LateralBreakingFoldFactor=1;%Will divide the off rate for all lateral bonds breaking by this factor, initialized and always set as 1. Reset if conditions are met.

end
rng('shuffle');%create random number generator seed for catastrophe...
cortime = 1000; %This determines how often a "snapshot" of the simulation is capture and stored for downstream use. 1000 is equivalent to around 1-2 seconds of real simulation time average. Can increase if you want the simulation to run faster and don't need as much temporal resolution.
%for debugging:
%DistanceFromLatAndEnds=zeros(nIterations,3);%NOT SURE, RETURN SJG. LIkely can be tossed without any repercussion.
% % if GTPEvL==1 1-25-24SJG return to once I find proper script. 
% % GTPEdge=zeros(nIterations/cortime)-1;%get this value every cortime steps, added 5-10-23
% % GTPLat=zeros(nIterations/cortime)-1;%Get this value every cortime steps, added 5-10-23
% % end

if movie_on == 1;%deciding if the movie will be generated or not...
    aviobj = VideoWriter(videoOutputName, 'Uncompressed AVI');
    aviobj.FrameRate = 10;  %frame rate in frames/sec
    open(aviobj);  %now opening the avi object
end

percDec1N=1; % multiplication factor on-rate if one neighboring subunits are longer%was0.9; for paper, used 1 %as of 10-10-23 was 5, likely for splaying
percDec2N=1; % multiplication factor on-rate if both neighboring subunits are longer%Was0.5; for paper, used 1% as of 10-10-23 was 25 likely for splaying
klatbond=100; %lateral bond formation rate (1/sec) %was 100, set back to 100 on 8-9-22; was 40 on 10-10-23, likelly for splaying..
if splayingOccurs==2 %to make more sites ahead of the highest lateral bond. This is for the extreme splaying. 
    percDec1N=4; % multiplication factor on-rate if one neighboring subunits are longer%was0.9; for paper, used 1 %as of 10-10-23 was 5, likely for splaying
    percDec2N=16; % multiplication factor on-rate if both neighboring subunits are longer%Was0.5; for paper, used 1% as of 10-10-23 was 25 likely for splaying
    klatbond=40; %lateral bond formation rate (1/sec) %was 100, set back to 100 on 8-9-22; was 40 on 10-10-23, likelly for splaying..

end

klatbondseam=100; %lateral bond formation rate at MT seam (1/sec)
klatbreak_TT=70; %lateral bond breakage rate (GTP/GTP) (1/sec)
klatbreak_TD=90; %lateral bond breakage rate (GTP/GDP) (1/sec)
klatbreak_DD=400; %lateral bond breakage rate (GDP/GDP) (1/sec)
klatbreakseam_TT=140; %lateral bond breakage rate at MT seam (GTP/GTP) (1/sec)
klatbreakseam_TD=180; %lateral bond breakage rate at MT seam (GTP/GDP) (1/sec)
klatbreakseam_DD=1200; %lateral bond breakage rate at MT seam (GDP/GDP) (1/sec)
testNum=0;
checkTagsBool=0; %UNSURE, double check
onPenaltyBool=1;%Used to determine if tubulin on-rate should be altered by presence of neighbor tubulins.
useBelow=0;  % method of calculating kon penalty for steric hinderance, legacy..
useHydro=1; %If you want hydrolylsis in the simulation. Should always be 1..but can change if you want. Was only used with hydrolysis =1.

protOn=1; %EB1 can arrive at MT's
fixed=0; % Fixed: percentage of subunits with tpx2 at start of simulation. Legacy code but for safety kept (likely can be removed if associated portions are also removed).
plotBool=1;%If you wnat to plot the output of each cortime loops; should keep but can set to zero and should increase speed.


protC=protConcUM*(10^(-6)); % [protein] in M,
kPlus=kPlusMT/filN; %tubulin on per protofilament uM-1 s-1 pf-1

fittedMTLength=zeros(nIterations/cortime,1);%Fitted microtubule length (average protofilament lenght)


kOn=kPlus*tubGtpC;%on-rate constant s-1, per protofilament
kOnProtLattice=kProtLattice*protC*(10^(9));%Eb1 on rate at lattice cites in sites-1 s-1
kOnProtEdge=kProtEdge*protC*(10^(9));%Eb1 on rate at edge cites in sites-1 s-1
kOff=1;%initialize off-rate constant, calculated later

kshorten_TBelow = 0.2; %1/sec Tubulin off rate when GTP dimer below
kshorten_DBelow = 200; %1/sec Tubulin off rate when GDP dimer below

%more initializations
startL=25; %starting seed length (will be all GTP behind it).
fixedSeedSize=25;%similar to startL. Could likely replace one with the other in a future code addition.
startHydro=0;%If you want the MT to be GDP to start rather than GTP
onesCol=ones(filN,1);%Don't remember, come back to.
highestFullHydro=startL;%To keep track of the highest layer in the microtubule that has hydrolyzed (for bookkeeping of GDP vs GTP tubulin dimers)
hydP=Inf;%initialize the time to hydroyslsis, will be calculated per event later.
jx=0;%variable to keep track of X position along the microtubule (PF #)
jy=0;%variable to keep track of y position along the microtubule (height position)
jy2=0;%variable to keep track of y position along the microtubule when looking at hydrolysis of the tubulin dimer (height position)
disP=0;%Initialize variable for calculating time for tubulin dissociation.
startArrayLen=300;%Number of items to initialize (y positions) in the underlying matrices that keep track of the MT and EB conditions. Can be increased if you want, but started arbitrarily at 300.
numProt=0; %keeps track of the number of protein that have bound. Likely legacy and can be removed, but kept for now.

% creates array of microtubule elements, 2 dimensional matrix of elements
%3rd dimension: layer 1: tubulin hydrolyzed (0/1) 0=GTPTub, 1=GDPTub; at
%end, for fig4, label non-existant tub subunits in matrix with 3/4
%               layer 2: tubulin tagged (0/1)
%               layer 3: protein bound (0/1)
%               layer 4: NOT USED!  lat bond candidate (0/1) -- 0 if tub below not yet  Right-side laterally bonded or if contains R lat. bonds
%               layer 5: tubulin has RIGHT lateral bond;  (0/1)0=No Bond
tubAProps=zeros(filN,startArrayLen,5);
matrixSize=size(tubAProps);%Generated so that the tubAProps can be increased in size as needed.
% tubAProps(:,fixedSeedSize+1,4)=1; %first row after the seed is candidate for lateral bond.
TrackSmallest=zeros(nIterations,3);%A variable used to track the fastest event and where it occured.
highestFullLatMade = startL*ones(1,filN); %For each protofilament, this is the index of the highest laterally bonded (Right side) subunit

if(startHydro)  %not used right now, starting length is GDP if on
    tubAProps(:,1:startL,1)=1;%makes them all hydrolyzed at the start if you want; never used for the simulation.
end
%%For EB1 protein dynamics, added by Sam below
%Added by Sam. This is to keep track of EB1 binding positions
%This matrix will keep track of the positions states.
%Dimension 3 layer 1: binding site Lattice/Edge single/ edge double
%lateral/edge double horizontal (and triple too) (0/1/2/3). I know that
%case 1 is probably overcoding but can allow for differences to be added
%later so seemed worth setting up.
%            layer 2: dimer hydrolysis state:
%            all GTP, split GTP/GDP, all GDP
%            (0/1/2) so lateral edges are either 0 or 2. (never 1)
%            layer 3: EB1 bound (0/1)
%It will have an element for each area across the MT that keeps track if the element is an
%edge binding spot (less than 4 dimers at the binding site) or a binding
%spot in the lattice (and therefore 4 dimers around the binding site). For
%layer 2, currently only three possibilities but can be changed at a later
%time to include GDP-Pi simply by adding more possibilities.
EB1BindingPositions=zeros(filN+1,startArrayLen,3);
if(startHydro)
    EB1BindingPositions(:,1:startL,2)=2; %This makes it so every binding site has GDP dimers associated with its base rather than
    %GTP since all tubulin dimers at the start are initialized to GDP
    %rather than GTP.
end

%These values are for the binning of proteins to better understand the
%dynamics.
binSize=10;%determining bin size, part of legacy code but left behind.
binValues=zeros(nIterations/cortime,30);%part of legacy code but left behind.
binValues(:,:)=-1; %initializing so you can see where no values are versus where zero proteins are bound part of legacy code but left behind.
rowSums=zeros(1,30);%variable to place summed bins to then store into binValues. part of legacy code but left behind.
outputData=zeros(nIterations/cortime,11);%This is to keep track of MT length, pf difference (taper), proteins bound, and proteins bound in tapered portion, and protiens bound 30 dimers before taper, %GTP edges bound, %GTP lattice bound, %GDP edges bound, %GDP lattice bound. part of legacy code but left behind.
proteinBindingArea=0; %Initializing this which will state which region of the MT should have the protein bind.
proteinRemovalArea=0; %Initializing this which will state which region of the MT should have a protein fall off of if this is the fastest action.
bindcount=1;%initializing value to keep track of proteins
proteinTracking=zeros(nIterations/2,10);%initialize matrix to hold where proteins bind and fell off from; keeps track of all the proteins for the downstream applications.
protein=0;%for debugging
falloff=0;%fore debugging
%removals=0;%for debugging
proteinPlaces=cell(nIterations/cortime,1);
proteinPlacesWithEvL=cell(nIterations/cortime,1);
protosLengths=zeros(nIterations/cortime,13);
presentProteins=[];%zeros(nIterations/cortime,4);%Only have X and Y, all other data kept in the other matrix. works with proteinTracking.
%presentProteins(:,:)=-1;%This is so everything is not a useful value to start (since zero is useful for some of the later portions).
tubAProps(:,1:startL,4)=1; %This is defining each tubulin dimer as present so that they can be accounted for later for EB1 binding

%%For protein binning to see if it reaches steady state.

%Initializing the EB1BindingPositions.
for jx=1:filN+1%This setup only works for Taylor's code because it starts flat at the plus end...this will likely need to be different for other code if they start not blunt.
    if jx==1 %This is to account for the special case of the seam upfront on the left side of PF1
        %EB1BindingPositions(jx,1,1)=1; %makes the one binding spot not at
        %the - end for PF1 that is an edge due to PF13 and Pf1 offset
        %actually an edge. Adjust to 1:2 if you move EB1bindingpositions
        %down one to account for complete - end. %no binding at seed so
        %removed this.
        for l=fixedSeedSize:startL; %go across the whole area to account for the seam as well..
            if tubAProps(filN,l-1,4)==0 && tubAProps(jx,l+1,1)==1% not sure if adding here or elsewhere... && tubAProps(jx,l+1,4)==1; %added second and so only binds if two spots
                EB1BindingPositions(jx,l,1)=2; %currently same but can change as needed if you adjust effects for 3, 2, or 1 dimers at edgebinding site. Set so binding only occurs at extnesion of 13 or 1 but never when equal.
            elseif tubAProps(filN,l-1,4)==0 && tubAProps(jx,l+1,1)==0
                EB1BindingPositions(jx,l,1)=1; %While this won't have binding, I coded it in in case I wanted to change it later
            else
                EB1BindingPositions(jx,l,1)=0;%redundant so don't need but don't like leaving an else without something...
            end
        end
    elseif jx>1 && jx<=filN %this is for the middle PFs
        for l=startL:-1:startL-24 %could assign summed value to a local variable so only needs to be checked once if slow..%%%%CHANGE THISSSS
            if tubAProps(jx,l,4)+tubAProps(jx,l+1,4)+tubAProps(jx-1,l,4)+tubAProps(jx-1,l+1,4)==4
                EB1BindingPositions(jx,l,1)=0;
            elseif tubAProps(jx,l,4)+tubAProps(jx,l+1,4)+tubAProps(jx-1,l,4)+tubAProps(jx-1,l+1,4)==3
                EB1BindingPositions(jx,l,1)=3;%three dimers currently considered horizontal edge, but all treated the same
            elseif tubAProps(jx,l,4)+tubAProps(jx,l+1,4)==2 && tubAProps(jx-1,l,4)==0%after and probably redundant but not sure so kept
                EB1BindingPositions(jx,l,1)=2; %so have a lateral pocket
            elseif tubAProps(jx,l,4)+tubAProps(jx-1,l,4)==2
                EB1BindingPositions(jx,l,1)=3; %so horizontal pocket
            else
                EB1BindingPositions(jx,l,1)=1; %so when you have a single tub dimer at pocket....just for bookkeeping to add into analysis later if wanted.
            end
        end
    else %This is for the last PF so jx=filN+1
        for l=fixedSeedSize:startL
            if tubAProps(1,l+2,4)==0 && tubAProps(filN,l+1,4)==1
                EB1BindingPositions(jx,l,1)=2;%lat dimers pair
            elseif tubAProps(1,l+2,4)==0 && tubAProps(filN,l+1,4)==0 %single tub dimer for pocket...ignored now but overcoded just in case
                EB1BindingPositions(jx,l,1)=1;
            else
                EB1BindingPositions(jx,l,1)=0;
            end
        end
    end
end

%%Above is for EB1 protein dynamics added by Sam.8-8-22

%keep track of proto-filament length and chance for association
pfLen=zeros(1,filN);
pfLenI=zeros(nIterations,filN);
% RLatBondI5=zeros(nIterations,100);
pfAssoc=zeros(1,filN);

%Fill first n rows
pfLen(:)=startL; % MT starts as 25 subunits long
pfLenI(1,:)=startL; % MT starts s 25 subunits long

%more initialization
shortestFilLength=zeros(nIterations,1);
longestFilLength=zeros(nIterations,1);
avgFilLength=zeros(nIterations,1);
stdev=zeros(nIterations,1);
chosenEventsList=zeros(nIterations,1);
timeX=zeros(nIterations,1);
highestFullH=zeros(nIterations,1);
numProtI=zeros(nIterations,1);


%initial conditions
timeX(1)=0;%sum up total time elapsed
longestFilLength(1)=startL;
avgFilLength(1)=startL;
stdev(1)=std(pfLen(1,:));
highestFullH(1)=startL;
numProtI(1)=0;
numlatbonds = 0;
if splayingOccurs==1
    %The goal of this code is to assign sites past the highest lat bond as edge
    %sites...
    for splaying=1:size(pfLen,2) %of note, highest lat bond is to right of point, so pf 1 highest bond is between pfs 1 and 2 (which aligns with EB1 binding position x=2 since EB1 bindign position starts at left edge of PF1 with #1 and goes to right edge of PF13 with #14
        EB1BindingPositions(splaying+1,highestFullLatMade(splaying):pfLen(splaying)-1,1)=2;%So, the binding pocket above the hightest full lat made will always be an edge because the next dimers won't be bound to make a lattice. Also, it won't be a single edge as there are tubulin dimers in both nearby PFs
        EB1BindingPositions(splaying+1,highestFullLatMade(splaying)-1:-1:startL,1)=0;%Everything below the highest full lat will for certain be a lattice site.
        if pfLen(splaying)>highestFullLatMade(splaying)
            EB1BindingPositions(splaying+1,pfLen(splaying),1)=1;%Okay, so in the event where PFlen>the lateral bond, the very last dimer will be a single edge and therefore a value of 1.
        end
    end
    EB1BindingPositions(1,startL:highestFullLatMade(size(pfLen,2))+1,1)=0;%since there is 1.5 dimer difference, dimer 1 on pf 13 blocks dimer 1, 2, and halfway up dimer 3 (but end of dimer 3 still available for binding).
    EB1BindingPositions(1,highestFullLatMade(size(pfLen,2))+2:pfLen(1)-1,1)=2; %since there is splaying, after the fullest lat bond everything but the last point will be an edge.
    if pfLen(1)>highestFullLatMade(size(pfLen,2))+2
        EB1BindingPositions(1,pfLen(1),1)=1; %So, if the PF is longer than the highest lat bond, then the very last point should be a single edge.
    end
    %%%Added by SJG for SPlaying on 1-17-23

end

%% Determine Which Subunits protein is bound at beginning of simulation.  Protein has 50% chance of binding (b/c fixed=0.5)
if(fixed>0)%This is legacy and not used (we had fixed as 0 for all EB1 runs).
    for jx=1:filN+1
        for jy=1:startL
            if(rand<fixed)
                %tubAProps(jx,jy,3)=1;
                EB1BindingPositions(jx,jy,3)=1;
                numProt=numProt+1;
            end
        end
    end
end

%This is beginning the actual simulation steps.
for i= 2:nIterations
    %reset all event probabilities/timesteps
    pfAssoc(:)=0;%association
    %%To remove edge or lattice binding partway through the simulation
    %%Added for doing videos with edge/lat on disappearing partway through
    %%simulation. If removeEdgeOrLatBinding=0, no changes. If 1, edge binding
    %%removed, if 2, lattice binding removed. Of note, the value for i is
    %%arbitrary and should be changed as the user sees fit...
    if removeEdgeOrLatBinding>0
        if i==100000
            if removeEdgeOrLatBinding==1
                kOnProtEdge=0;
            else %removeEdgeOrLatBinding==2
                kOnProtLattice=0;
            end


        end
        if i==150000

            if removeEdgeOrLatBinding==1
                kOnProtEdge=kProtEdge*protC*(10^(9));
            else %removeEdgeOrLatBinding==2
                kOnProtLattice=kProtLattice*protC*(10^(9));
            end


        end

        if i==350000
            if removeEdgeOrLatBinding==1
                kOnProtEdge=0;
            else %removeEdgeOrLatBinding==2
                kOnProtLattice=0;
            end
        end
        if i==400000
            if removeEdgeOrLatBinding==1
                kOnProtEdge=kProtEdge*protC*(10^(9));
            else %removeEdgeOrLatBinding==2
                kOnProtLattice=kProtLattice*protC*(10^(9));
            end
        end
    end

    %%So now, just calculated all of the possible edge and lattice sites with
    %%the loops and then removed the number of proteins bound in each area.
    %%Now, need to keep track of the number of Edge Proteins and Lat Proteins.
    %smallestV and smallest keep track of smallest time step and associated
    %event.%%%SJG pick up here on 1-25-24
    smallest=[1,1,1];%[x,y,location]
    %in case of assocVect, x and y will be set to vector pos.
    smallestV=Inf;
    if protOn==1
        EdgeBindingPositions=[conditionMatrixadjusted2(EB1BindingPositions,3,0,0,fixedSeedSize,pfLen);conditionMatrixadjusted2(EB1BindingPositions,3,1,0,fixedSeedSize,pfLen);conditionMatrixadjusted2(EB1BindingPositions,3,2,0,fixedSeedSize,pfLen);...
            conditionMatrixadjusted2(EB1BindingPositions,2,0,0,fixedSeedSize,pfLen);conditionMatrixadjusted2(EB1BindingPositions,2,1,0,fixedSeedSize,pfLen);conditionMatrixadjusted2(EB1BindingPositions,2,2,0,fixedSeedSize,pfLen)];
        LatticeBindingPositions=[conditionMatrixadjusted(EB1BindingPositions(2:filN,:,:),0,0,0,fixedSeedSize,pfLen);conditionMatrixadjusted(EB1BindingPositions(2:filN,:,:),0,1,0,fixedSeedSize,pfLen);conditionMatrixadjusted(EB1BindingPositions(2:filN,:,:),0,2,0,fixedSeedSize,pfLen)];%ignored seams cause they should not have binding spots in the lattice
    else
        EdgeBindingPositions=0;
        LatticeBindingPositions=0;
    end
    for jx= 1:filN
        nN=0;%number of neighbors
        nNB=0;%number of neighbors below
        depth=max(startL+1,pfLen(jx)-20); % Doesn't let depth get below start length; 20 max length of tip to analyze

        %% Calculate Wait Time: Subunit Dissociation: 1) Determine koff 2) calculate wait time

        calculate_koff_flag = 0;  %this flags whether or not koff is actually calculated (are there any koff candidates, if not koff not assigned)

        %%%For edge and lat on rates, doing it wrong currently where calculate ALL
        %%%edge sites for each PF, but will only do this calculateion once per loop
        %%%if put up here rather than 13 times...

        %%  calculate koff %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if(jx==1)%Proto #1
            for jy=pfLen(jx):-1:highestFullLatMade(jx)+1


                % Determine koff. RULE: Subunit can only leave if NO
                % Right lateral bonds AND <2 Left lat bonds with neighbors
                if tubAProps(jx,jy,5)==0 && tubAProps(filN,jy,5)<2; %pf 13 can have 1 lateral bond because its lower bond
                    calculate_koff_flag = 1; %koff will be calculated, so will calculate an event time
                    if tubAProps(jx,jy-1,1)==1 %Tub below is GDP-tub
                        kOff=kshorten_DBelow;
                    else %Tub below is GTP-tub
                        kOff=kshorten_TBelow;
                    end

                else

                    jy = jy + 1;

                    break
                end



                if(onPenaltyBool)
                    if(pfLen(jx)<pfLen(jx+1))
                        nN=nN+1;
                    elseif(pfLen(jx)<pfLen(jx+1)+1)
                        nNB=nNB+1;
                    end
                    if(pfLen(jx)<=pfLen(filN)+1)
                        nN=nN+1;
                    elseif(pfLen(jx)<=pfLen(filN)+2)
                        nNB=nNB+1;
                    end

                end

            end

        elseif(jx==filN)%if on matrix right end (1st protofilament)
            for jy=pfLen(jx):-1:highestFullLatMade(jx)+1

                % Determine koff. RULE: Subunit can only leave if <2 Right
                % lateral bonds AND NO left lat bonds with neighbors
                if tubAProps(jx,jy,5)<2 && tubAProps(jx-1,jy,5)==0;
                    calculate_koff_flag = 1; %koff will be calculated, so will calculate an event time
                    if tubAProps(jx,jy-1,1)==1 %Tub below is GDP-tub
                        kOff=kshorten_DBelow;
                    else %Tub below is GTP-tub
                        kOff=kshorten_TBelow;
                    end
                else

                    jy = jy + 1;

                    break

                end

            end

            if(onPenaltyBool)
                if(pfLen(jx)<pfLen(jx-1))
                    nN=nN+1;
                elseif(pfLen(jx)<pfLen(jx-1)+1)
                    nNB=nNB+1;
                end
                if(pfLen(jx)<=pfLen(1)-2)
                    nN=nN+1;
                elseif(pfLen(jx)<=pfLen(1)-1)
                    nNB=nNB+1;
                end

            end


        else %if in middle
            for jy=pfLen(jx):-1:highestFullLatMade(jx)+1  %depth is the number of subunits away from the tip that we care about; defined above
                % Determine koff. RULE: Subunit can only leave if NO
                % lateral bonds with Left or Right neighbor

                if tubAProps(jx,jy,5)==0 && tubAProps(jx-1,jy,5)==0
                    calculate_koff_flag = 1; %koff will be calculated, so will calculate an event time
                    if tubAProps(jx,jy-1,1)==1 %Tub below is GDP-tub
                        kOff=kshorten_DBelow;
                    else %Tub below is GTP-tub
                        kOff=kshorten_TBelow;
                    end


                else
                    jy = jy + 1;

                    break


                end

            end


            if(onPenaltyBool)
                if(pfLen(jx)<pfLen(jx-1))
                    nN=nN+1;
                elseif(pfLen(jx)<pfLen(jx-1)+1)%why not just do greater or equal too, i believe it does the same thing here
                    nNB=nNB+1;
                end
                if(pfLen(jx)<pfLen(jx+1))
                    nN=nN+1;
                elseif(pfLen(jx)<pfLen(jx+1)+1)
                    nNB=nNB+1;
                end
            end

        end
        %assign associated timestep of dissociation for each cell

        if calculate_koff_flag == 1; %there is a subunit that is a candidate for dissociation, so koff was assigned
            disP=-log(rand)/kOff; % Wait time: subunit dissociation
            if(disP<smallestV)
                smallestV=disP; %actual wait time for off-rate

                valueToUse  = randi([jy pfLen(jx)]);%this assumes all dimers have same off rate, rather than only picking dimers that are say GTP or GDP.
                smallest=[jx,valueToUse,1]; %position of subunit with smallest off-rate wait time
                %removals=removals+1;
                %DistanceFromLatAndEnds(removals,1)=pfLen(jx);
                % DistanceFromLatAndEnds(removals,2)=highestFullLatMade(jx);
                %DistanceFromLatAndEnds(removals,3)=valueToUse;

            end
        end


        %% Calculate Wait Time: EB1 Association- set protOn=0 (off)--should skip this section

        if(fixed==0 && protOn)
            EdgeProtP=-log(rand)/(kOnProtEdge*size(EdgeBindingPositions,1));%((conditionCounter(EB1BindingPositions,2,2,0,fixedSeedSize,max(pfLen))+conditionCounter(EB1BindingPositions,2,0,0,fixedSeedSize,max(pfLen))+...

            LatticeProtP=-log(rand)/(kOnProtLattice*size(LatticeBindingPositions,1));%((conditionCounter(EB1BindingPositions(2:filN,:,:),0,1,0,fixedSeedSize,max(pfLen))+conditionCounter(EB1BindingPositions(2:filN,:,:),0,2,0,fixedSeedSize,max(pfLen))+conditionCounter(EB1BindingPositions(2:filN,:,:),0,0,0,fixedSeedSize,max(pfLen)))));%%%%This is wrong...currently overcalculating.

            %%This below could definitely be cleaner by outsourcing to a
            %%function but I think it will work so I can come back and clean
            %%up after the fact
            protP=min([EdgeProtP,LatticeProtP]);


            if protP==EdgeProtP
                proteinBindingArea=1;
            elseif protP==LatticeProtP
                proteinBindingArea=2;
            else
                proteinBindingArea=0;
            end

            if(protP<smallestV)
                smallestV=protP;
                smallest=[jx,0,4];
            end

        end
        %}
        %% Calculate Wait Time: Association of new subunits
        %Association determination.
        if(testNum==0 || i<testNum || jx~=1) %testNum=0,
            if(useBelow) %keep turned off
                if(nN==0)
                    pfAssoc(jx)=-log(rand)/(kOn-0.1*kOn*percDec1N*nNB);
                elseif(nN==1)
                    pfAssoc(jx)=-log(rand)/(kOn-kOn*percDec1N-0.1*kOn*percDec1N*nNB);
                else
                    pfAssoc(jx)=-log(rand)/(kOn-kOn*percDec2N);
                end
            else %Use This
                if(nN==0)
                    pfAssoc(jx)=-log(rand)/(kOn);
                elseif(nN==1)

                    pfAssoc(jx)=-log(rand)/(kOn*percDec1N);


                else

                    pfAssoc(jx)=-log(rand)/(kOn*percDec2N);


                end
            end
        end

        if(pfAssoc(jx)<smallestV)
            smallestV=pfAssoc(jx);
            smallest=[jx,pfLen(jx)+1,3];
        end

        %% Calculate Wait Time: Lateral Bond Formation
        %         for jy=pfLen(jx):-1:fixedSeedSize + 1

        calculate_lateral_bond_formation = 0;  %this flags whether or not lateral bond formation time is actually calculated)

        if jx==filN %13th proto


            if highestFullLatMade(jx)+1 < pfLen(jx)  %not end subunit
                if (tubAProps(jx,highestFullLatMade(jx)+1,5)==0) && (highestFullLatMade(jx)+1<=pfLen(1)) %pf 13 subunit does not have a bond and pf 1 is long enough
                    latbondP= -log(rand)/klatbondseam;
                    calculate_lateral_bond_formation = 1;
                    %if has 1 lat bond to pf#1
                elseif (tubAProps(jx,highestFullLatMade(jx)+1,5)==1) && (highestFullLatMade(jx)+1<=pfLen(1)+1)%pf 13 subunit has one lower bond and pf 1 is long enough
                    latbondP= -log(rand)/klatbondseam;
                    calculate_lateral_bond_formation = 1;
                end
            end

        else

            if highestFullLatMade(jx)+1 < pfLen(jx)
                if (highestFullLatMade(jx)+1<=pfLen(jx+1))
                    latbondP= -log(rand)/klatbond;
                    calculate_lateral_bond_formation = 1;
                end
            end

        end

        if calculate_lateral_bond_formation == 1; %will calculate a lateral bond formation event time
            if(latbondP<smallestV)
                smallestV=latbondP;
                smallest=[jx,highestFullLatMade(jx)+1,6];
            end
        end

        %% Calculate Wait Time: Lateral Bond Breakage

        subunitbreak=0;
        calculate_lateral_bond_breakage = 0;  %this flags whether or not lateral bond breakage time is actually calculated

        if highestFullLatMade(jx)>startL %Can't break lateral bonds past start length (was seed size)
            % RULE: If subunit has maximum lateral bonds (If PF 2-12: 2
            % TOTAL; If PFs 1,13: 3 TOTAL), lat. bond breakage rate =
            % kbreak/Pibreak

            calculate_lateral_bond_breakage = 1;  %will calculate a breakage time

            if (jx==filN)  %For Proto #13

                if tubAProps(jx,highestFullLatMade(jx)+1,5)==1 %if the subunit ahead has only 1of2 possible lat bonds (b/c highest FullLatMade for proto13 is highest subunit with 2 right lat bonds)

                    subunitbreak=highestFullLatMade(jx)+1; %subunit w/lat bond to be broken; will be added to Smallest variable later


                    %Determine if lat bond between subunits
                    %(pf 13, highestFullLatMade(jx)+1) and
                    %(pf 1, highestFullLatMade(jx)+1) are GXP/GXP interactions:

                    %if GDP-GTP Lat Interaction at Seam
                    if (tubAProps(jx,highestFullLatMade(jx)+1,1)==1) && (tubAProps(1,highestFullLatMade(jx)+1,1)==0)
                        latbreakP=-log(rand)/klatbreakseam_TD/LateralBreakingFoldFactor;

                        %if GTP-GDP Lat Interaction at Seam
                    elseif (tubAProps(jx,highestFullLatMade(jx)+1,1)==0) && (tubAProps(1,highestFullLatMade(jx)+1,1)==1)
                        latbreakP=-log(rand)/klatbreakseam_TD/LateralBreakingFoldFactor;

                        %if GDP-GDP Lat Interaction at Seam
                    elseif (tubAProps(jx,highestFullLatMade(jx)+1,1)==1) && (tubAProps(1,highestFullLatMade(jx)+1,1)==1)
                        latbreakP=-log(rand)/klatbreakseam_DD/LateralBreakingFoldFactor;

                        %if GTP-GTP Lat Interaction at Seam
                    else %(tubAProps(jx,highestFullLatMade(jx)+1,1)==0) && (tubAProps(1,highestFullLatMade(jx)+1,1)==0)
                        latbreakP=-log(rand)/klatbreakseam_TT/LateralBreakingFoldFactor;

                    end

                else  %highestFullLatMade(jx)+1,5) NOT equal to 1, so Look at FULL bonds (=2) at highestFullLatMade(jx); 2 RIGHT LATERAL BONDS EXIST

                    subunitbreak=highestFullLatMade(jx);

                    if (tubAProps(jx-1,highestFullLatMade(jx),5)==1) %if MAX Lateral Bonds exist-- LEFT lateral bond AND 2 RIGHT lateral bonds
                        %Determine if lat bond between subunits
                        %(pf13,highestFullLatMade(jx)) and
                        %(pf1,highestFullLatMade(jx)+1) are GXP/GXP interactions:

                        %if GDP-GTP Lat Interaction at Seam


                        if (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(1,highestFullLatMade(jx)+1,1)==0)
                            latbreakP=-log(rand)/(klatbreakseam_TD/Pibreak)/LateralBreakingFoldFactor;

                            %if GTP-GDP Lat Interaction at Seam
                        elseif (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(1,highestFullLatMade(jx)+1,1)==1)
                            latbreakP=-log(rand)/(klatbreakseam_TD/Pibreak)/LateralBreakingFoldFactor;

                            %if GDP-GDP Lat Interaction at Seam
                        elseif (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(1,highestFullLatMade(jx)+1,1)==1)
                            latbreakP=-log(rand)/(klatbreakseam_DD/Pibreak)/LateralBreakingFoldFactor;

                            %if GTP-GTP Lat Interaction at Seam
                        else %(tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(1,highestFullLatMade(jx)+1,1)==0)
                            latbreakP=-log(rand)/(klatbreakseam_TT/Pibreak)/LateralBreakingFoldFactor;
                        end

                    else %BOTH lateral bonds do NOT exist -- *NO LEFT lateral bonds* and 2 RIGHT lateral bonds

                        %if GDP-GTP Lat Interaction at Seam
                        if (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(1,highestFullLatMade(jx)+1,1)==0)
                            latbreakP=-log(rand)/klatbreakseam_TD/LateralBreakingFoldFactor;

                            %if GTP-GDP Lat Interaction at Seam
                        elseif (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(1,highestFullLatMade(jx)+1,1)==1)
                            latbreakP=-log(rand)/klatbreakseam_TD/LateralBreakingFoldFactor;

                            %if GDP-GDP Lat Interaction at Seam
                        elseif (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(1,highestFullLatMade(jx)+1,1)==1)
                            latbreakP=-log(rand)/klatbreakseam_DD/LateralBreakingFoldFactor;

                            %if GTP-GTP Lat Interaction at Seam
                        else  % (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(1,highestFullLatMade(jx)+1,1)==0)
                            latbreakP=-log(rand)/klatbreakseam_TT/LateralBreakingFoldFactor;
                        end

                    end

                end

            elseif jx==1 %For Proto #1
                subunitbreak=highestFullLatMade(jx);
                if tubAProps(filN,highestFullLatMade(jx),5)>=1  && tubAProps(filN,highestFullLatMade(jx)-1,5)==2  %Check if has MAX lat bonds: ie has 3 lat bonds (2left AND 1 right)--EXAMPLE: subunit 6 max bonds would be if bonded to (PF13,sub5), (PF13,sub6), and (PF2,sub6)
                    %if GDP-GTP Lat Interaction
                    if (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(jx+1,highestFullLatMade(jx),1)==0)
                        latbreakP=-log(rand)/(klatbreak_TD/Pibreak)/LateralBreakingFoldFactor;

                        %if GTP-GDP Lat Interaction
                    elseif (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(jx+1,highestFullLatMade(jx),1)==1)
                        latbreakP=-log(rand)/(klatbreak_TD/Pibreak)/LateralBreakingFoldFactor;

                        %if GDP-GDP Lat Interaction
                    elseif (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(jx+1,highestFullLatMade(jx),1)==1)
                        latbreakP=-log(rand)/(klatbreak_DD/Pibreak)/LateralBreakingFoldFactor;

                        %if GTP-GTP Lat Interaction
                    else% (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(jx+1,highestFullLatMade(jx),1)==0)
                        latbreakP=-log(rand)/(klatbreak_TT/Pibreak)/LateralBreakingFoldFactor;
                    end

                else %if has 1 lat bonds (1 Right lat bond Only)
                    %if GDP-GTP Lat Interaction
                    if (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(jx+1,highestFullLatMade(jx),1)==0)
                        latbreakP=-log(rand)/klatbreak_TD/LateralBreakingFoldFactor;

                        %if GTP-GDP Lat Interaction
                    elseif (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(jx+1,highestFullLatMade(jx),1)==1)
                        latbreakP=-log(rand)/klatbreak_TD/LateralBreakingFoldFactor;

                        %if GDP-GDP Lat Interaction
                    elseif (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(jx+1,highestFullLatMade(jx),1)==1)
                        latbreakP=-log(rand)/klatbreak_DD/LateralBreakingFoldFactor;

                        %if GTP-GTP Lat Interaction
                    else% (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(jx+1,highestFullLatMade(jx),1)==0)
                        latbreakP=-log(rand)/klatbreak_TT/LateralBreakingFoldFactor;
                    end

                end
            else %For Protos #2-12
                subunitbreak=highestFullLatMade(jx);
                if tubAProps(jx-1,highestFullLatMade(jx),5)==1 %if has LEFT lat bond (ie. 2 TOTAL lat bonds)
                    %if GDP-GTP Lat Interaction
                    if (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(jx+1,highestFullLatMade(jx),1)==0)
                        latbreakP=-log(rand)/(klatbreak_TD/Pibreak)/LateralBreakingFoldFactor;

                        %if GTP-GDP Lat Interaction
                    elseif (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(jx+1,highestFullLatMade(jx),1)==1)
                        latbreakP=-log(rand)/(klatbreak_TD/Pibreak)/LateralBreakingFoldFactor;

                        %if GDP-GDP Lat Interaction
                    elseif (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(jx+1,highestFullLatMade(jx),1)==1)
                        latbreakP=-log(rand)/(klatbreak_DD/Pibreak)/LateralBreakingFoldFactor;

                        %if GTP-GTP Lat Interaction
                    else% (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(jx+1,highestFullLatMade(jx),1)==0)
                        latbreakP=-log(rand)/(klatbreak_TT/Pibreak)/LateralBreakingFoldFactor;
                    end

                else %if has R lat bonds only
                    %if GDP-GTP Lat Interaction
                    if (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(jx+1,highestFullLatMade(jx),1)==0)
                        latbreakP=-log(rand)/klatbreak_TD/LateralBreakingFoldFactor;

                        %if GTP-GDP Lat Interaction
                    elseif (tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(jx+1,highestFullLatMade(jx),1)==1)
                        latbreakP=-log(rand)/klatbreak_TD/LateralBreakingFoldFactor;

                        %if GDP-GDP Lat Interaction
                    elseif (tubAProps(jx,highestFullLatMade(jx),1)==1) && (tubAProps(jx+1,highestFullLatMade(jx),1)==1)
                        latbreakP=-log(rand)/klatbreak_DD/LateralBreakingFoldFactor;

                        %if GTP-GTP Lat Interaction
                    else %(tubAProps(jx,highestFullLatMade(jx),1)==0) && (tubAProps(jx+1,highestFullLatMade(jx),1)==0)
                        latbreakP=-log(rand)/klatbreak_TT/LateralBreakingFoldFactor;
                    end

                end

            end

        end

        if calculate_lateral_bond_breakage == 1  %will calculate a breakage time

            if(latbreakP<smallestV) %Update Wait Time
                smallestV=latbreakP;

                smallest=[jx,subunitbreak,7];
            end

        end


    end




    %% Calculate wait time for EB1 dissociation


    if fixed ==0 && sum(EB1BindingPositions(:,:,3),'all')>0&&protOn==1

        GTPEdgeProtP1=-log(rand)/(kPOffEdgeGTP*size(presentProteins(presentProteins(:,8)>0&presentProteins(:,9)==0),1));

        GTPLatticeProtP1=-log(rand)/(kPOffLatticeGTP*size(presentProteins(presentProteins(:,8)==0&presentProteins(:,9)==0),1));

        GDPEdgeProtP1=-log(rand)/(kPOffEdgeGDP*size(presentProteins(presentProteins(:,8)>0&presentProteins(:,9)>0),1));

        GDPLatticeProtP1=-log(rand)/(kPOffLatticeGDP*size(presentProteins(presentProteins(:,8)==0&presentProteins(:,9)>0),1));

        protR=min([GTPEdgeProtP1,GTPLatticeProtP1,GDPEdgeProtP1,GDPLatticeProtP1]);%This is to find the fastest area for protein binding
        if protR==GTPEdgeProtP1 %This will make note of the fastest area for protein binding where 1 is in GTP edge, 2 is GTP lattice, 3 is gdp edge, 4 is gdp lattice
            proteinRemovalArea=1;
        elseif protR==GTPLatticeProtP1
            proteinRemovalArea=2;
        elseif protR==GDPEdgeProtP1
            proteinRemovalArea=3;
        elseif protR==GDPLatticeProtP1
            proteinRemovalArea=4;
        else
            proteinRemovalArea=-1;
        end


        if(protR<smallestV)
            smallestV=protR;
            smallest=[0,0,5];
        end

    end


    %% Execute most likely event
    %         smallest
    if(smallest(3)==3) % Association event-- new subunit added at tubAProps(smallest(1),smallest(2),X)
        %Update Protofilament Length
        pfLen(smallest(1))=pfLen(smallest(1))+1;
        tubAProps(smallest(1),pfLen(smallest(1)),4)=1;%This now says the dimer is present
        %Update Lateral Bonds
        %*Update Candidacy for Lateral Bond Formation:* If tubulin subunit
        %added is the immediately above the seed OR if it is immediately
        %above a laterally bonded subunit (note: tracking by lateral bond
        %to right)



        if(fixed>0 && rand<fixed) %this adds a protein with the associating tubulin subunit
            tubAProps(smallest(1),smallest(2),3)=1;
            numProt=numProt+1;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %The second portion of the or statement was added by me, not from
        %vanburen paper.
        %Oh, and what this does is tag a GTP unit if it attaches on a GDP
        %unit.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if(checkTagsBool)
            tubAProps(smallest(1),smallest(2),2)=tubAProps(smallest(1),smallest(2)-1,1)||tubAProps(smallest(1),smallest(2)-1,2);

            untagList=CheckTags(smallest(1),smallest(2),tubAProps,pfLen(:),filN,1);
            nTags=size(untagList);
            for num=1:nTags(1); %untag if near 2 neighbors
                tubAProps(untagList(num,1),untagList(num,2),2)=0;
            end
        end
        chosenEventsList(i)=1;

        %%%Below is added by SJG on 8-8-22 to deal with adjusting
        %%%edge/lattice state for tubulin dimers WRT EB1 binding.
        %Now, the lattice versus edge status along the MT needs to be adjusted. This also adjusts the hydrolysis state
        %for the binding pockets to ensure that they are update when a new
        %dimer is added in that is important for the hydrolysis state (so
        %it's at the back of a pocket (- end).
        if smallest(1)==1 %This is for the left side of PF1, along the seam
            EB1BindingPositions(smallest(1),smallest(2),2)=tubAProps(smallest(1),smallest(2),1)*2;%fixing hydrolysis states
            EB1BindingPositions(smallest(1)+1,smallest(2),2)=(tubAProps(smallest(1),smallest(2),1)+tubAProps(smallest(1)+1,smallest(2),1))*tubAProps(smallest(1)+1,smallest(2),4)+...
                (tubAProps(smallest(1),smallest(2),4)-tubAProps(smallest(1)+1,smallest(2),4))*(tubAProps(smallest(1),smallest(2),1)*2);%Fixing hydrolysis states


            for m=0:1%This is to change the edge/lat state
                for n=0:1
                    if m==0
                        if tubAProps(filN,pfLen(smallest(1))-1-n,4)==0 && tubAProps(smallest(1),pfLen(smallest(1))+1-n,4)==1
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=2;
                        elseif tubAProps(filN,pfLen(smallest(1))-1-n,4)==0 && tubAProps(smallest(1),pfLen(smallest(1))+1-n,4)==0
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=1;
                        else
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=0;
                        end
                        if tubAProps(filN,pfLen(smallest(1))-2-n,4)==1
                            EB1BindingPositions(filN+1,pfLen(smallest(1))-2-n,1)=0;
                            EB1BindingPositions(filN+1,pfLen(smallest(1))-2-n,3)=0;%This has a protein fall off of it was at the edge...
                            falloff=falloff+1;
                            j=proteinTracking(:,1)==filN+1&proteinTracking(:,2)==pfLen(smallest(1))-2-n&proteinTracking(:,5)==0;

                            proteinTracking(j,5)=timeX(i-1);
                            proteinTracking(j,6)=6;
                            proteinTracking(j,10)=6*proteinTracking(j,3);

                        end
                    else
                        if (tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4))+tubAProps(smallest(1)+m,pfLen(smallest(1))-n+1,4)+...
                                tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n+1,4)==4
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=0; %This is for in the lattice
                        elseif (tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4))+tubAProps(smallest(1)+m,pfLen(smallest(1))-n+1,4)+...
                                tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n+1,4)==3
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=3; %This is for when you have three dimers, currently classified as horizontal partners
                        elseif tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)+m,pfLen(smallest(1))-n+1,4)==2 ||...
                                tubAProps(smallest(1)+m-1,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)+m-1,pfLen(smallest(1))-n+1,4)==2
                            if m==0 && tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)==0;
                                EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=2;%This is two vertical dimers
                            elseif m==1 && tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4)==0;
                                EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=2;
                            end
                        elseif tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)==2;
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=3; %This is two horizontal dimers
                        else
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=1; %This is one tub for pocket...overcoding just in case for later
                        end
                    end
                end
            end
        elseif smallest(1)>1 && smallest(1)<filN %This will need to be verified but i think it looks at the four corners correctly
            EB1BindingPositions(smallest(1),smallest(2),2)=(tubAProps(smallest(1),smallest(2),1)+tubAProps(smallest(1)-1,smallest(2),1))*...
                tubAProps(smallest(1)-1,smallest(2),4)+(tubAProps(smallest(1),smallest(2),4)-tubAProps(smallest(1)-1,smallest(2),4))*...
                (tubAProps(smallest(1),smallest(2),1)*2);%to set up the hydrolysis state for binding pocket

            EB1BindingPositions(smallest(1)+1,smallest(2),2)=(tubAProps(smallest(1),smallest(2),1)+tubAProps(smallest(1)+1,smallest(2),1))*...
                tubAProps(smallest(1)+1,smallest(2),4)+(tubAProps(smallest(1),smallest(2),4)-tubAProps(smallest(1)+1,smallest(2),4))*...
                (tubAProps(smallest(1),smallest(2),1)*2);%to set up the hydrolysis state for binding pocket
            for m=0:1
                for n=0:1
                    if (tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4))+tubAProps(smallest(1)+m,pfLen(smallest(1))-n+1,4)+...
                            tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n+1,4)==4
                        EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=0; %This is for in the lattice
                    elseif (tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4))+tubAProps(smallest(1)+m,pfLen(smallest(1))-n+1,4)+...
                            tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n+1,4)==3
                        EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=3; %This is for when you have three dimers, currently classified as horizontal partners
                    elseif tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)+m,pfLen(smallest(1))-n+1,4)==2||...
                            tubAProps(smallest(1)+m-1,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)+m-1,pfLen(smallest(1))-n+1,4)==2
                        if m==0 && tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)==0;
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=2;%This is two vertical dimers
                        elseif m==1 && tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4)==0;
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=2;
                        end%This is two vertical dimers
                    elseif tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)==2;
                        EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=3; %This is two horizontal dimers
                    else
                        EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=1; %This is one tub for pocket...overcoding just in case for later
                    end
                end
            end
        else %This is for the PF 13, along the seam and along the 12 PF edge
            EB1BindingPositions(smallest(1),smallest(2),2)=(tubAProps(smallest(1),smallest(2),1)+tubAProps(smallest(1)-1,smallest(2),1))*...
                tubAProps(smallest(1)-1,smallest(2),4)+(tubAProps(smallest(1),smallest(2),4)-tubAProps(smallest(1)-1,smallest(2),4))*...
                (tubAProps(smallest(1),smallest(2),1)*2);%This is for hydrolysis state
            EB1BindingPositions(smallest(1)+1,smallest(2),2)=tubAProps(smallest(1),smallest(2),1)*2;%This is for hydrolysis state.
            if splayingOccurs==1
                EB1BindingPositions(1,smallest(2)+1,3)=0;%ADDED BY SJG on 1-17-23 to remove any EB1 bound when the site is no longer an edge site. The plus 1 could be wrong...

            end


            for m=0:1
                for n=0:1
                    if m==0
                        if (tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4))+tubAProps(smallest(1)+m,pfLen(smallest(1))-n+1,4)+...
                                tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n+1,4)==4
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=0; %This is for in the lattice
                        elseif (tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4))+tubAProps(smallest(1)+m,pfLen(smallest(1))-n+1,4)+...
                                tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n+1,4)==3
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=3; %This is for when you have three dimers, currently classified as horizontal partners
                        elseif tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)+m,pfLen(smallest(1))-n+1,4)==2||...
                                tubAProps(smallest(1)+m-1,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)+m-1,pfLen(smallest(1))-n+1,4)==2
                            if m==0 && tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)==0;
                                EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=2;%This is two vertical dimers
                            elseif m==1 && tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4)==0;
                                EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=2;
                            end%This is two vertical dimers
                        elseif tubAProps(smallest(1)+m,pfLen(smallest(1))-n,4)+tubAProps(smallest(1)-1+m,pfLen(smallest(1))-n,4)==2;
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=3; %This is two horizontal dimers
                        else
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=1; %This is one tub for pocket...overcoding just in case for later
                        end
                    else
                        if tubAProps(1,pfLen(smallest(1))+2-n,4)==0 && tubAProps(filN,pfLen(smallest(1))+1-n,4)==1;
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=2;
                        elseif tubAProps(1,pfLen(smallest(1))+2-n,4)==0 && tubAProps(filN,pfLen(smallest(1))+1-n,4)==0;
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=1;
                        else
                            EB1BindingPositions(smallest(1)+m,pfLen(smallest(1))-n,1)=0;
                        end
                        if tubAProps(1,pfLen((smallest(1)))+1-n,4)==1
                            EB1BindingPositions(1,pfLen((smallest(1)))+1-n,1)=0;
                            EB1BindingPositions(1,pfLen((smallest(1)))+1-n,3)=0; %This is for protein to fall off when pf1 pocket is messed up by addition of dimer on pf13
                            falloff=falloff+1;
                            j=proteinTracking(:,1)==1 & proteinTracking(:,2)==pfLen(smallest(1))+1-n & proteinTracking(:,5)==0;

                            proteinTracking(j,5)==0;
                            proteinTracking(j,5)=timeX(i-1);
                            proteinTracking(j,6)=6;
                            proteinTracking(j,10)=6*proteinTracking(j,3);

                        end
                    end
                end
            end
            %%%find the indices that need to be changed in
            %%%proteinTracking...
        end
        for m=0:1
            for n=0:1
                jj=proteinTracking(:,1)==(smallest(1)+m) & proteinTracking(:,2)==(smallest(2)-n) & proteinTracking(:,5)==0;
                ii=find(proteinTracking(:,1)==smallest(1)+m & proteinTracking(:,2)==smallest(2)-n & proteinTracking(:,5)==0);%gets index in proteinTracking to change as needed for all spots affected by tubulin addition

                if(ii>0)
                    proteinTracking(ii,9)=EB1BindingPositions(smallest(1)+m,smallest(2)-n,2);%should adjust all 4 spots for GTP state (don't need all four, only need two but still
                    proteinTracking(ii,8)=EB1BindingPositions(smallest(1)+m,smallest(2)-n,1);%should adjust all 4 spots for edge/latt state
                    proteinTracking(jj,12)=EB1BindingPositions(smallest(1)+m,smallest(2)-n,2);
                    proteinTracking(jj,11)=EB1BindingPositions(smallest(1)+m,smallest(2)-n,1);

                end
            end
        end
    elseif(smallest(3)==1)  %Dissociation Event

        if smallest(2) < startL; % dont depolymerize below start length for each protofilament

            smallest(2) = startL;

        end

        %% protein dissociation
        %%ADDED by SJG on 8-8-22 for EB1.
        tubAProps(smallest(1),smallest(2):pfLen(smallest(1)),4)=0;%This makes the dimers be marked as not present for determining binding pockets

        for k=smallest(2)-1:pfLen(smallest(1))%might be able to do this with the other loops to determine lattice/edge but not sure so did seperately and can try to combine later
            for m=0:1
                for n=0:1
                    if smallest(1)==1 %first PF
                        if m==0
                            if tubAProps(smallest(1)+m,k-n,4)+tubAProps(smallest(1)+m,k+1-n,4)==0

                                EB1BindingPositions(smallest(1),k-n,3)=0;
                                falloff=falloff+1;
                                j=proteinTracking(:,5)==0 & proteinTracking(:,1)==smallest(1)+m&proteinTracking(:,2)==k-n;

                                proteinTracking(j,5)=timeX(i-1);
                                proteinTracking(j,6)=5;
                                proteinTracking(j,10)=5+6*(proteinTracking(j,3)-1);%makes it either 5 or 11 for edge bindingn or lat binding
                            end%No else here because want the state to stay the same whether protein bound or not
                        else
                            if tubAProps(smallest(1)+m,k-n,4)+tubAProps(smallest(1)+m,k+1-n,4)+...
                                    tubAProps(smallest(1)+m-1,k-n,4)+tubAProps(smallest(1)+m-1,k+1-n,4)<2%%%SAM EDIT

                                EB1BindingPositions(smallest(1)+m,k-n,3)=0;
                                falloff=falloff+1;
                                %protLost=protLost+1;
                                j=proteinTracking(:,5)==0&proteinTracking(:,1)==smallest(1)+m&proteinTracking(:,2)==k-n;

                                proteinTracking(j,5)=timeX(i-1);
                                proteinTracking(j,6)=5;
                                proteinTracking(j,10)=5+6*(proteinTracking(j,3)-1);

                            end
                        end
                    elseif smallest(1)>1 && smallest(1)<filN %middle PFs
                        if tubAProps(smallest(1)+m,k-n,4)+tubAProps(smallest(1)+m,k+1-n,4)+...
                                tubAProps(smallest(1)+m-1,k-n,4)+tubAProps(smallest(1)+m-1,k+1-n,4)<2%%%SAM EDIT

                            EB1BindingPositions(smallest(1)+m,k-n,3)=0;
                            falloff=falloff+1;
                            j=proteinTracking(:,5)==0&proteinTracking(:,1)==smallest(1)+m&proteinTracking(:,2)==k-n;

                            proteinTracking(j,5)=timeX(i-1);
                            proteinTracking(j,6)=5;
                            proteinTracking(j,10)=5+6*(proteinTracking(j,3)-1);

                        end%No else here because don't want to change protein state if still a place to bind
                    else %This will be for the last PF.
                        if m==0
                            if tubAProps(smallest(1)+m,k-n,4)+tubAProps(smallest(1)+m,k+1-n,4)+...
                                    tubAProps(smallest(1)+m-1,k-n,4)+tubAProps(smallest(1)+m-1,k+1-n,4)<2%SAM EDIT

                                EB1BindingPositions(smallest(1)+m,k-n,3)=0;
                                %protLost=protLost+1;
                                falloff=falloff+1;
                                j=proteinTracking(:,5)==0&proteinTracking(:,1)==smallest(1)+m&proteinTracking(:,2)==k-n;

                                proteinTracking(j,5)=timeX(i-1);
                                proteinTracking(j,6)=5;
                                proteinTracking(j,10)=5+6*(proteinTracking(j,3)-1);

                            end
                        else
                            if tubAProps(smallest(1)+m-1,k-n,4)+tubAProps(smallest(1)+m-1,k+1-n,4)==0%+...


                                EB1BindingPositions(filN+1,k-n,3)=0;%This should work even though I set it up differently for where in EB1 binding positions relative to other loops.

                                falloff=falloff+1;
                                j=proteinTracking(:,5)==0&proteinTracking(:,1)==smallest(1)+m&proteinTracking(:,2)==k-n;

                                proteinTracking(j,5)=timeX(i-1);
                                proteinTracking(j,6)=5;
                                proteinTracking(j,10)=5+6*(proteinTracking(j,3)-1);

                            end
                        end
                    end
                end
            end

        end


        %% Tubulin Dissociation
        %%Again, added by SJG on 8-8-22.


        tubAProps(smallest(1),smallest(2):pfLen(smallest(1)),:)=0; %reset all Tub Properties to zero
        pfLen(smallest(1))=smallest(2)-1; %recalculate protofilament length
        chosenEventsList(i)=-1;

        if(smallest(2)<=highestFullHydro)  %if shortened below the maximum full hydro
            highestFullHydro=smallest(2)-1; %resets max hydrolength
        end
        %%%edge versus latteral when dimer falls off. Also adjusts hydrollysis
        %%%state for where dimers fell off.
        for x=smallest(2):pfLen(smallest(1))

            if smallest(1)==1 %This is for the left side of PF1, along the seam
                EB1BindingPositions(smallest(1),x,2)=tubAProps(smallest(1),x,1)*2;%This is for changing hydrolysis states
                EB1BindingPositions(smallest(1)+1,x,2)=(tubAProps(smallest(1),x,1)+tubAProps(smallest(1)+1,x,1))*tubAProps(smallest(1)+1,x,4)+...
                    (tubAProps(smallest(1),x,4)-tubAProps(smallest(1)+1,x,4))*(tubAProps(smallest(1),x,1)*2);

                for m=0:1%%This is for changing lattice/edge states
                    for n=0:1
                        if m==0
                            if tubAProps(filN,x-1-n,4)==0 && tubAProps(smallest(1),x+1-n,4)==1
                                EB1BindingPositions(smallest(1)+m,x-n,1)=2;
                            elseif tubAProps(filN,x-1-n,4)==0 && tubAProps(smallest(1),x+1-n,4)==0 && tubAProps(smallest(1),x-n,4)==1
                                EB1BindingPositions(smallest(1)+m,x-n,1)=1;

                            else
                                EB1BindingPositions(smallest(1)+m,x-n,1)=0;
                            end
                            if tubAProps(filN,x-n-2,4)==1 && tubAProps(smallest(1),x-n,4)==0 && tubAProps(filN,x-n-1,4)==1
                                EB1BindingPositions(filN+1,x-n-2,1)=2;
                            elseif tubAProps(filN,x-n-2,4)==1 && tubAProps(smallest(1),x-n,4)==0 && tubAProps(filN,x-n-1,4)==0
                                EB1BindingPositions(filN+1,x-n-2,1)=1;
                            end
                            %ADJUST LAST PF STATE AS WELL...i think above
                            %fixes this issue
                            %this if statement setup not validated yet...
                        else
                            if (tubAProps(smallest(1)+m,x-n,4))+tubAProps(smallest(1)+m,x-n+1,4)+...
                                    tubAProps(smallest(1)-1+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n+1,4)==4 || ...
                                    (tubAProps(smallest(1)+m,x-n,4))+tubAProps(smallest(1)+m,x-n+1,4)+...
                                    tubAProps(smallest(1)-1+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n+1,4)==0
                                EB1BindingPositions(smallest(1)+m,x-n,1)=0; %This is for in the lattice
                            elseif (tubAProps(smallest(1)+m,x-n,4))+tubAProps(smallest(1)+m,x-n+1,4)+...
                                    tubAProps(smallest(1)-1+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n+1,4)==3
                                EB1BindingPositions(smallest(1)+m,x-n,1)=3; %This is for when you have three dimers, currently classified as horizontal partners
                            elseif tubAProps(smallest(1)+m,x-n,4)+tubAProps(smallest(1)+m,x-n+1,4)==2 &&...
                                    tubAProps(smallest(1)-1+m,x-n,4)==0;
                                EB1BindingPositions(smallest(1)+m,x-n,1)=2;%This is two vertical dimers
                            elseif tubAProps(smallest(1)+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n,4)==2;
                                EB1BindingPositions(smallest(1)+m,x-n,1)=3; %This is two horizontal dimers
                            else
                                EB1BindingPositions(smallest(1)+m,x-n,1)=1; %This is one tub for pocket...overcoding just in case for later
                            end
                        end
                    end
                end
            elseif smallest(1)>1 && smallest(1)<filN %This will need to be verified but i think it looks at the four corners correctly
                EB1BindingPositions(smallest(1),x,2)=(tubAProps(smallest(1),x,1)+tubAProps(smallest(1)-1,x,1))*...
                    tubAProps(smallest(1)-1,x,4)+(tubAProps(smallest(1),x,4)-tubAProps(smallest(1)-1,x,4))*...
                    (tubAProps(smallest(1),x,1)*2);%This is for hydro state

                EB1BindingPositions(smallest(1)+1,x,2)=(tubAProps(smallest(1),x,1)+tubAProps(smallest(1)+1,x,1))*...
                    tubAProps(smallest(1)+1,x,4)+(tubAProps(smallest(1),x,4)-tubAProps(smallest(1)+1,x,4))*...
                    (tubAProps(smallest(1),x,1)*2);%This is for hydro state
                for m=0:1%%This is for edge/lattice state
                    for n=0:1
                        if tubAProps(smallest(1)+m,x-n,4)+tubAProps(smallest(1)+m,x-n+1,4)+...
                                tubAProps(smallest(1)+m-1,x-n,4)+tubAProps(smallest(1)+m-1,x-n+1,4)==0
                            EB1BindingPositions(smallest(1),x,:)=0;%%%MADE EDITS HERE...
                        else
                            if (tubAProps(smallest(1)+m,x-n,4))+tubAProps(smallest(1)+m,x-n+1,4)+...
                                    tubAProps(smallest(1)-1+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n+1,4)==4|| ...
                                    (tubAProps(smallest(1)+m,x-n,4))+tubAProps(smallest(1)+m,x-n+1,4)+...
                                    tubAProps(smallest(1)-1+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n+1,4)==0
                                EB1BindingPositions(smallest(1)+m,x-n,1)=0; %This is for in the lattice
                            elseif (tubAProps(smallest(1)+m,x-n,4))+tubAProps(smallest(1)+m,x-n+1,4)+...
                                    tubAProps(smallest(1)-1+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n+1,4)==3
                                EB1BindingPositions(smallest(1)+m,x-n,1)=3; %This is for when you have three dimers, currently classified as horizontal partners
                            elseif tubAProps(smallest(1)+m,x-n,4)+tubAProps(smallest(1)+m,x-n+1,4)==2 &&...
                                    tubAProps(smallest(1)-1+m,x-n,4)==0;
                                EB1BindingPositions(smallest(1)+m,x-n,1)=2;%This is two vertical dimers
                            elseif tubAProps(smallest(1)+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n,4)==2;
                                EB1BindingPositions(smallest(1)+m,x-n,1)=3; %This is two horizontal dimers
                            else
                                EB1BindingPositions(smallest(1)+m,x-n,1)=1; %This is one tub for pocket...overcoding just in case for later
                            end
                        end
                    end
                end
            else %This is for the PF 13, along the seam and along the 12 PF edge
                EB1BindingPositions(smallest(1),x,2)=(tubAProps(smallest(1),x,1)+tubAProps(smallest(1)-1,x,1))*...
                    tubAProps(smallest(1)-1,x,4)+(tubAProps(smallest(1),x,4)-tubAProps(smallest(1)-1,x,4))*...
                    (tubAProps(smallest(1),x,1)*2);%This is for hydro state
                EB1BindingPositions(smallest(1)+1,x,2)=tubAProps(smallest(1),x,1)*2;%This is for hydro state
                for m=0:1%This is for lattice/edge determination
                    for n=0:1
                        if m==0
                            if (tubAProps(smallest(1)+m,x-n,4))+tubAProps(smallest(1)+m,x-n+1,4)+...
                                    tubAProps(smallest(1)-1+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n+1,4)==4|| ...
                                    (tubAProps(smallest(1)+m,x-n,4))+tubAProps(smallest(1)+m,x-n+1,4)+...
                                    tubAProps(smallest(1)-1+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n+1,4)==0
                                EB1BindingPositions(smallest(1)+m,x-n,1)=0; %This is for in the lattice
                            elseif (tubAProps(smallest(1)+m,x-n,4))+tubAProps(smallest(1)+m,x-n+1,4)+...
                                    tubAProps(smallest(1)-1+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n+1,4)==3
                                EB1BindingPositions(smallest(1)+m,x-n,1)=3; %This is for when you have three dimers, currently classified as horizontal partners
                            elseif tubAProps(smallest(1)+m,x-n,4)+tubAProps(smallest(1)+m,x-n+1,4)==2 &&...
                                    tubAProps(smallest(1)-1+m,x-n,4)==0;
                                EB1BindingPositions(smallest(1)+m,x-n,1)=2;%This is two vertical dimers
                            elseif tubAProps(smallest(1)+m,x-n,4)+tubAProps(smallest(1)-1+m,x-n,4)==2;
                                EB1BindingPositions(smallest(1)+m,x-n,1)=3; %This is two horizontal dimers
                            else
                                EB1BindingPositions(smallest(1)+m,x-n,1)=1; %This is one tub for pocket...overcoding just in case for later
                            end
                        else
                            if tubAProps(1,x+2-n,4)==0 && tubAProps(filN,x+1-n,4)==1;
                                EB1BindingPositions(smallest(1)+m,x-n,1)=2;
                            elseif tubAProps(1,x+2-n,4)==0 && tubAProps(filN,x+1-n,4)==0 && tubAProps(filN,x-n,4)==1;
                                EB1BindingPositions(smallest(1)+m,x-n,1)=1;
                            else
                                EB1BindingPositions(smallest(1)+m,x-n,1)=0;
                            end
                            if tubAProps(1,x+1-n,4)==1 && tubAProps(filN,x-n,4)==0 && tubAProps(1,x+2-n,4)==1
                                EB1BindingPositions(1,x+1-n,1)=2;
                            elseif tubAProps(1,x+1-n,4)==1 && tubAProps(filN,x-n,4)==0 && tubAProps(1,x+2-n,4)==0
                                EB1BindingPositions(1,x+1-n,1)=1;
                            end %this if statement setup not validated yet...
                            %ADJUST PF1 WHEN DImer ON PF13 FALLS OFF
                        end
                    end
                end
            end
            %this is to adjust the protein GTP and edge state after tubulin
            %falls off. It will do it for each dimer that falls off due to
            %using k so will have overlap when multiple dimers fall off but
            %seemed easiest to code for. Can be expanded out thought to be
            %more efficient if needed.
            for m=0:1
                for n=0:1
                    jj=proteinTracking(:,1)==(smallest(1)+m) & proteinTracking(:,2)==(x-n) & proteinTracking(:,5)==0;
                    ii=find(proteinTracking(:,1)==smallest(1)+m & proteinTracking(:,2)==x-n & proteinTracking(:,5)==0);%gets index in proteinTracking to change as needed for all spots affected by tubulin addition

                    if(ii>0)
                        proteinTracking(ii,9)=EB1BindingPositions(smallest(1)+m,x-n,2);%should adjust all 4 spots for GTP state (don't need all four, only need two but still
                        proteinTracking(ii,8)=EB1BindingPositions(smallest(1)+m,x-n,1);%should adjust all 4 spots for edge/latt state
                        proteinTracking(jj,12)=EB1BindingPositions(smallest(1)+m,x-n,2);
                        proteinTracking(jj,11)=EB1BindingPositions(smallest(1)+m,x-n,1);

                    end
                end
            end
        end



    elseif(smallest(3)==6) %lateral bond formation
        numlatbonds = numlatbonds +1;


        if (smallest(1)==filN) %PF#13-each subunit n may have 0,1,or2 lat. interactions with subunits n+1 and n+2 in PF#1
            if tubAProps(smallest(1),smallest(2),5)==0
                tubAProps(smallest(1),smallest(2),5)=1; %Lateral Bond (Right) now exists
                highestFullLatMade(smallest(1)) = smallest(2)-1; %(shouldnt change)
                if splayingOccurs==1
                    %%ADded by SJG to have protein fall off when lateral bond
                    %%appears on 1-18-23
                    if EB1BindingPositions(1,smallest(2)+1,3)==1
                        EB1BindingPositions(1,smallest(2)+1,3)=0;%
                        usableIndex=proteinTracking(:,1)==1&proteinTracking(:,2)==highestFullLatMade(smallest(1))+2&proteinTracking(:,5)==0;
                        proteinTracking(usableIndex,5)=timeX(i-1);
                        %%Added by SJG on 1-18-23
                    end
                    %Second part added as kept having proteins still bound two
                    %beind the tip, not sure why and thought this would deal
                    %with that....
                    if EB1BindingPositions(1,smallest(2),3)==1
                        EB1BindingPositions(1,smallest(2),3)=0;
                        usableIndex=proteinTracking(:,1)==1&proteinTracking(:,2)==highestFullLatMade(smallest(1))+1&proteinTracking(:,5)==0;
                        proteinTracking(usableIndex,5)=timeX(i-1);
                    end
                end

            elseif tubAProps(smallest(1),smallest(2),5)==1
                tubAProps(smallest(1),smallest(2),5)=2; %second Lateral Bond (Right) now exists
                %                 tubAProps(smallest(1),smallest(2),4)=0; %NO longer a candidate for lat bond formation
                %Update Candidacy (TubAProps, layer4): after lat bond forms, if the next tubulin subunit
                %exists in proto, then it is now a Candidate.
                %                 if pfLen(smallest(1))>smallest(2)
                %                     tubAProps(smallest(1),smallest(2)+1,4)=1;
                %                 end
                %Track highest Right lateral bond made in each protofilament
                highestFullLatMade(smallest(1)) = smallest(2);
                if splayingOccurs==1
                    if EB1BindingPositions(size(pfLen,2)+1,smallest(2),3)==1
                        EB1BindingPositions(size(pfLen,2)+1,smallest(2),3)=0;
                        usableIndex=proteinTracking(:,1)==size(pfLen,2)+1&proteinTracking(:,2)==highestFullLatMade(smallest(1))&proteinTracking(:,5)==0;
                        proteinTracking(usableIndex,5)=timeX(i-1);
                    end
                    %Second part added as kept having proteins still bound two
                    %beind the tip, not sure why and thought this would deal
                    %with that....
                    if EB1BindingPositions(size(pfLen,2)+1,smallest(2)-1,3)==1
                        EB1BindingPositions(size(pfLen,2)+1,smallest(2)-1,3)=0;
                        usableIndex=proteinTracking(:,1)==size(pfLen,2)+1&proteinTracking(:,2)==highestFullLatMade(smallest(1))-1&proteinTracking(:,5)==0;
                        proteinTracking(usableIndex,5)=timeX(i-1);
                    end
                    %%Added by SJG on 1-18-23
                end


            end


        else %PF#1-12
            tubAProps(smallest(1),smallest(2),5)=1; %Lateral Bond (Right) now exists
            %             tubAProps(smallest(1),smallest(2),4)=0; %NO longer a candidate for lat bond formation
            %Update Candidacy (TubAProps, layer4): after lat bond forms, if the next tubulin subunit
            %exists in proto, then it is now a Candidate.
            %             if pfLen(smallest(1))>smallest(2)
            %                 tubAProps(smallest(1),smallest(2)+1,4)=1;
            %             end
            %Track highest Right lateral bond made in each protofilament
            highestFullLatMade(smallest(1)) = smallest(2);
        end

    elseif(smallest(3)==7) %lateral bond breakage
        numlatbonds = numlatbonds-1;
        if (smallest(1)==filN) %Proto #13

            tubAProps(smallest(1),smallest(2),5)=tubAProps(smallest(1),smallest(2),5)-1; %Lateral Bond (Right) now DOES NOT exist

            if tubAProps(smallest(1),smallest(2),5)==1 %if only 1of2 lat bonds is now present, update highestFullLatMade.
                highestFullLatMade(smallest(1)) = smallest(2)-1;
            end



        else %Proto #1-12
            tubAProps(smallest(1),smallest(2),5)=0; %Lateral Bond (Right) now DOES NOT exist

            highestFullLatMade(smallest(1)) = smallest(2)-1;
        end




    elseif(smallest(3)==4)%protein adds on ADDED by SJG on 8-8-22, copied from Version 10 of Rebecca's function.
        numProt=numProt+1;
        if proteinBindingArea==1%%%%%ERROR IN THIS AREA FOR EB1 BInding #


            jz=randi(size(EdgeBindingPositions,1));
            EB1BindingPositions(EdgeBindingPositions(jz,1),EdgeBindingPositions(jz,2),3)=1;
            a=1;
        elseif proteinBindingArea==2
            %Adjusted so that it only looks up to the length of each PF so
            %that it does not get to areas where max PF length is longer
            %than multiple pfs and therefore has a binding spot since 0 is
            %default state and is the state used to define binding spot...
            %bindingPositions=LatticeBindingPositions;%possiblebindings
            %LatticeBindingPositions=[conditionMatrixadjusted(EB1BindingPositions(2:filN,:,:),0,0,0,fixedSeedSize,pfLen);conditionMatrixadjusted(EB1BindingPositions(2:filN,:,:),0,1,0,fixedSeedSize,pfLen);conditionMatrixadjusted(EB1BindingPositions(2:filN,:,:),0,2,0,fixedSeedSize,pfLen];%ADDED by SJG on 8-10-22 so this is only calculated if binding will occur rather than every iteration.
            jz=randi(size(LatticeBindingPositions,1));%%%Got an error here cause no binding posiitions, make sure consistent here and at start..
            EB1BindingPositions(LatticeBindingPositions(jz,1),LatticeBindingPositions(jz,2),3)=1;%check back in here but should be fixed i think...
        else %proteinBindingArea==0
        end

        %don't think i need to sort by rows for this method.
        %tubAProps(smallest(1),jy,3)=1;
        chosenEventsList(i)=2;

        if proteinBindingArea==2
            proteinTracking(bindcount,1)=LatticeBindingPositions(jz,1);%%removed +1 here...+1;%This is supposed to give the X position, off by one because the lattice protien finding matrix only looks at edges 2-13 and is offset by one so this fixes the offset
            proteinTracking(bindcount,2)=LatticeBindingPositions(jz,2);
            proteinTracking(bindcount,3)=proteinBindingArea;
            proteinTracking(bindcount,4)=timeX(i-1);
            proteinTracking(bindcount,8)=EB1BindingPositions(LatticeBindingPositions(jz,1),LatticeBindingPositions(jz,2),1);
            proteinTracking(bindcount,9)=EB1BindingPositions(LatticeBindingPositions(jz,1),LatticeBindingPositions(jz,2),2);
            proteinTracking(bindcount,11)=EB1BindingPositions(LatticeBindingPositions(jz,1),LatticeBindingPositions(jz,2),1);
            proteinTracking(bindcount,12)=EB1BindingPositions(LatticeBindingPositions(jz,1),LatticeBindingPositions(jz,2),2);
            bindcount=bindcount+1;%for debugging
        else


            proteinTracking(bindcount,1)=EdgeBindingPositions(jz,1);
            proteinTracking(bindcount,2)=EdgeBindingPositions(jz,2);
            proteinTracking(bindcount,3)=proteinBindingArea;
            proteinTracking(bindcount,4)=timeX(i-1);
            proteinTracking(bindcount,8)=EB1BindingPositions(EdgeBindingPositions(jz,1),EdgeBindingPositions(jz,2),1);
            proteinTracking(bindcount,9)=EB1BindingPositions(EdgeBindingPositions(jz,1),EdgeBindingPositions(jz,2),2);
            proteinTracking(bindcount,11)=EB1BindingPositions(EdgeBindingPositions(jz,1),EdgeBindingPositions(jz,2),1);
            proteinTracking(bindcount,12)=EB1BindingPositions(EdgeBindingPositions(jz,1),EdgeBindingPositions(jz,2),2);
            bindcount=bindcount+1;%for debugging
        end
        protein=protein+1;
    else%protein comes off (smallest(3)==5)ADDED by SJG on 8-8-22, copied from Version 10 of Rebecca's function.
        if proteinRemovalArea==1
            possibleRemovals=presentProteins(presentProteins(:,8)>0&presentProteins(:,9)==0,:);%[conditionMatrix(EB1BindingPositions,3,0,1,fixedSeedSize,max(pfLen));conditionMatrix(EB1BindingPositions,2,0,1,fixedSeedSize,max(pfLen))];

            dim=size(possibleRemovals,1);
            jz=randi([1,dim]);
            while jz>size(possibleRemovals,1)
                jz=randi([1,dim]);

            end
            EB1BindingPositions(possibleRemovals(jz,1),possibleRemovals(jz,2),3)=0;
        elseif proteinRemovalArea==2
            possibleRemovals=presentProteins(presentProteins(:,8)==0&presentProteins(:,9)==0,:);%conditionMatrix(EB1BindingPositions,0,0,1,fixedSeedSize,max(pfLen));

            dim=size(possibleRemovals,1);
            jz=randi([1,dim]);
            while jz>size(possibleRemovals,1)
                jz=randi([1,dim]);

            end
            EB1BindingPositions(possibleRemovals(jz,1),possibleRemovals(jz,2),3)=0;
        elseif proteinRemovalArea==3
            possibleRemovals=presentProteins(presentProteins(:,8)>0&presentProteins(:,9)>0,:);%[conditionMatrix(EB1BindingPositions,3,1,1,fixedSeedSize,max(pfLen));conditionMatrix(EB1BindingPositions,3,2,1,fixedSeedSize,max(pfLen));...

            dim=size(possibleRemovals,1);
            jz=randi([1,dim]);
            while jz>size(possibleRemovals,1)
                jz=randi([1,dim]);

            end
            EB1BindingPositions(possibleRemovals(jz,1),possibleRemovals(jz,2),3)=0;
        elseif proteinRemovalArea==4
            possibleRemovals=presentProteins(presentProteins(:,8)==0&presentProteins(:,9)>0,:);%[conditionMatrix(EB1BindingPositions,0,1,1,fixedSeedSize,max(pfLen));conditionMatrix(EB1BindingPositions,0,2,1,fixedSeedSize,max(pfLen))];

            dim=size(possibleRemovals,1);
            jz=randi([1,dim]);
            while jz>size(possibleRemovals,1)
                jz=randi([1,dim]);

            end
            EB1BindingPositions(possibleRemovals(jz,1),possibleRemovals(jz,2),3)=0;
        else %proteinRemovalArea==0
        end
        kk=proteinTracking(:,1)==possibleRemovals(jz,1)&proteinTracking(:,2)==possibleRemovals(jz,2)&proteinTracking(:,5)==0;

        proteinTracking(kk,10)=proteinRemovalArea+6*(proteinTracking(kk,3)-1);%makes it so that you know the removal area

        chosenEventsList(i)=-2;

        j=proteinTracking(:,1)==possibleRemovals(jz,1)&proteinTracking(:,2)==possibleRemovals(jz,2)&proteinTracking(:,5)==0;

        proteinTracking(j,5)=timeX(i-1);
        proteinTracking(j,6)=proteinRemovalArea;

    end

    timeX(i)=timeX(i-1)+smallestV;%sum up total time elapsed
    pfLenI(i,:)=pfLen(:);

    highestFullH(i)=highestFullHydro;
    numProtI(i)=numProt;

    if(max(pfLen(:))>matrixSize(2)-5)
        tubAProps(:,end+1:end+50,:)=0;
        matrixSize(2)=matrixSize(2)+50;
    end


    TrackSmallest(i,:)=smallest;

    %% GTP Hydrolysis -- NOW separate from other wait time calculations.

    if(useHydro) %% GTP Hydrolysis

        for jx= 1:filN
            for jy2=pfLen(jx)-1:-1:highestFullHydro+1 % step backwards from pfLength-1 (non-terminal) to seed

                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %Should check to see if faster this way, or with if statements
                %first.
                %determine hydrolisis occurance
                if(~tubAProps(jx,jy2,1))&& (jy2 > startL) %if not hydrolyzed (i.e. GtP-tub; =0)
                    hydP=-log(rand)/kHyd;
                    if(hydP<smallestV) %if GTP hydrolysis wait time shorter than time that passed in i-iteration, then execute GTP hydrolysis
                        tubAProps(jx,jy2,1)=1; %hydrolyze
                        tubAProps(jx,jy2,2)=0; %untag if hydrolzed

                        %update highestFullHydro
                        for check=highestFullHydro+1:jy2
                            if(isequal(tubAProps(:,check,1),onesCol))
                                highestFullHydro=check;
                            else
                                break
                            end
                        end

                        %%Added by Sam 8-8-22 copied from Rebecca's version
                        %%10
                        if jx==1 %for dimers on first PF
                            EB1BindingPositions(jx,jy2,2)=tubAProps(jx,jy2,1)*2;
                            EB1BindingPositions(jx+1,jy2,2)=(tubAProps(jx,jy2,1)+tubAProps(jx+1,jy2,1))*tubAProps(jx+1,jy2,4)+...
                                (tubAProps(jx,jy2,4)-tubAProps(jx+1,jy2,4))*(tubAProps(jx,jy2,1)*2);
                        elseif jx>1 && jx<filN %for dimers in rest of lattice
                            EB1BindingPositions(jx,jy2,2)=(tubAProps(jx,jy2,1)+tubAProps(jx-1,jy2,1))*...
                                tubAProps(jx-1,jy2,4)+(tubAProps(jx,jy2,4)-tubAProps(jx-1,jy2,4))*...
                                (tubAProps(jx,jy2,1)*2);

                            EB1BindingPositions(jx+1,jy2,2)=(tubAProps(jx,jy2,1)+tubAProps(jx+1,jy2,1))*...
                                tubAProps(jx+1,jy2,4)+(tubAProps(jx,jy2,4)-tubAProps(jx+1,jy2,4))*...
                                (tubAProps(jx,jy2,1)*2);
                        else %for dimers on last PF
                            EB1BindingPositions(jx,jy2,2)=(tubAProps(jx,jy2,1)+tubAProps(jx-1,jy2,1))*...
                                tubAProps(jx-1,jy2,4)+(tubAProps(jx,jy2,4)-tubAProps(jx-1,jy2,4))*...
                                (tubAProps(jx,jy2,1)*2);
                            EB1BindingPositions(jx+1,jy2,2)=tubAProps(jx,jy2,1)*2;
                        end
                        %%ADDEd by Sam8-8-22 copied from Rebecca's version
                        %%10
                        for m=0:1
                            jj=proteinTracking(:,1)==(jx+m) & proteinTracking(:,2)==(jy2) & proteinTracking(:,5)==0;
                            ii=find(proteinTracking(:,1)==jx+m & proteinTracking(:,2)==jy2 & proteinTracking(:,5)==0);%gets index in proteinTracking to change as needed for all spots affected by tubulin addition

                            if(ii>0)
                                proteinTracking(ii,9)=EB1BindingPositions(jx+m,jy2,2);%should adjust all 4 spots for GTP state (don't need all four, only need two but still
                                proteinTracking(jj,12)=EB1BindingPositions(jx+m,jy2,2);


                            end
                        end
                    end
                end
            end
        end
    end

    presentProteins=[];
    presentProteins=proteinTracking(proteinTracking(:,1)>0,:);
    presentProteins=presentProteins(presentProteins(:,5)==0,:);
% % if GTPEvL==1 %%SJG 1-25-24 return once I find right EvL script to
% compare to. 
% %     if mod(i,cortime)==0 %while written clumsily, this here will after every cortime timesteps, will count the number of EB1 bound to GTP edge sites and to GTP lattice sites and record them so they can be used later. This will just give the total #. 
% %         GTPLat(i,cortime)=0;
% %         GTPEdge(i,cortime)=0;
% %         for protofilament=1:filN+1
% %             for mtLength=max(pfLen)
% %                 if EB1BindingPositions(protofilament,mtLength,3)==1 %EB1 bound
% %                     if EB1BindingPositions(protofilament,mtLength,2)==0 %GTP state
% %                         if EB1BindingPositions(protofilament,mtLength,1)==0
% %                             GTPLat(i/cortime)=GTPLat(i/cortime)+1;
% %                         elseif EB1BindingPositions(protofilament,mtLength,1)>1
% %                             GTPEdge(i/cortime)=GTPEdge(i/cortime)+1;
% %                 
% %                         end
% %                     end
% %                 end
% %             end
% %         end
% % 
% %     end
% % end



    %%%ADDED by SJG on 8-8-22; copied from Rebecca's version 10.
    if mod(i,cortime)== 0 && movie_on==1;
        fittedMTLength(i/cortime,1)=ceil(mean(pfLen(:,:)));
        i %Again, to show the output if you have no movie...SJG added 12-28-23
        if movie_on==1
            tubAPLOT=[];
            tubAPLOT = tubAProps;

            for j=1:filN
                tubAPLOT(j,pfLen(j)+1:end,1)=3;
            end
            %format numbers in 14-row matrix for plotting
            tubAProps2=zeros(14,length(tubAPLOT(j,:,1)));
            tubAProps2(1:filN,:)=tubAPLOT(:,:,1);
            tubAProps2(filN+1,:)=4*ones(1,length(tubAPLOT(1,:,1)));
            pcolor(tubAProps2(:,:,1))

            xlabel('Tubulin Subunits')
            ylabel('Protofilament')
            title('GTP Hydrolysis')
            ylim([0 filN+1])
            xlim([0 max(pfLen)+1])

            hold on
            plot(highestFullLatMade,[1.5:filN+0.5],'*r','MarkerSize',10)
            plot(pfLenI(end,:),[1.5:filN+0.5],'ow','MarkerSize',10)

            %%added by Sam
            index=presentProteins(:,3)==1;
            totalEdgeProteins=presentProteins(index,:);
            index=presentProteins(:,3)==2;
            totalLatticeProteins=presentProteins(index,:);
            if size(totalEdgeProteins)>0
                plot(totalEdgeProteins(:,2)+1,totalEdgeProteins(:,1),'*g','MarkerSize',10);%Sam adjusted what was plotted to a variable i created for TPX2 locations (currently) which can be changed to EB1 once the settings for EB1 are created
            end
            if size(totalLatticeProteins)>0
                plot(totalLatticeProteins(:,2)+1,totalLatticeProteins(:,1),'*m','MarkerSize',10);
            end
            %%added by Sam

            hold off

            frame = getframe(gcf);   %captures CURRENT image
            writeVideo(aviobj,frame);
        end
        rowSums=sum(EB1BindingPositions(:,:,3));
        for k=1:(ceil(size(EB1BindingPositions,2)/binSize))-1   %This will only do bins that are up to the bin at the greatest length where proteins can be.

            binValues(i/1000,k)=sum(rowSums(1,binSize*(k-1)+1:binSize*k));
        end
        binValues(i/1000,k+1)=sum(rowSums(1,binSize*(k)+1:end));

        outputData(i/cortime,1)=max(pfLen);%length
        outputData(i/cortime,2)=max(pfLen)-min(pfLen);%pfdiff
        outputData(i/cortime,3)=sum(rowSums);%#proteins
        outputData(i/cortime,4)=sum(rowSums(1,min(pfLen):max(pfLen)));%proteins in taper
        if min(pfLen)>31
            outputData(i/cortime,5)=sum(rowSums(1,min(pfLen)-31:min(pfLen)-1));%#proteins 30 subunits from back of taper without overlap
        end
        outputData(i/cortime,6)=((sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==0 &EB1BindingPositions(:,1:max(pfLen),3)==1&EB1BindingPositions(:,1:max(pfLen),1)==2))+sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==0 &EB1BindingPositions(:,1:max(pfLen),3)==1&EB1BindingPositions(:,1:max(pfLen),1)==3))));
        outputData(i/cortime,7)=(sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==0 &EB1BindingPositions(:,1:max(pfLen),3)==1&EB1BindingPositions(:,1:max(pfLen),1)==0)));%/sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==0&EB1BindingPositions(:,1:max(pfLen),1)==0)))*100;%get percent of GTP lattice bound out of total GTP lattice positions %was way off
        outputData(i/cortime,8)=outputData(i/cortime,6)/(outputData(i/cortime,6)+outputData(i/cortime,7));%%gtp edge proteins

        outputData(i/cortime,9)= ((sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==1 &EB1BindingPositions(:,1:max(pfLen),3)==1&EB1BindingPositions(:,1:max(pfLen),1)==2))+sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==1 &EB1BindingPositions(:,1:max(pfLen),3)==1&EB1BindingPositions(:,1:max(pfLen),1)==3))+...
            sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==2 &EB1BindingPositions(:,1:max(pfLen),3)==1&EB1BindingPositions(:,1:max(pfLen),1)==2)+sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==2 &EB1BindingPositions(:,1:max(pfLen),3)==1&EB1BindingPositions(:,1:max(pfLen),1)==3)))));
        outputData(i/cortime,10)=((sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==1 &EB1BindingPositions(:,1:max(pfLen),3)==1&EB1BindingPositions(:,1:max(pfLen),1)==0))+sum(sum(EB1BindingPositions(:,1:max(pfLen),2)==2 &EB1BindingPositions(:,1:max(pfLen),3)==1&EB1BindingPositions(:,1:max(pfLen),1)==0))));
        outputData(i/cortime,11)=outputData(i/cortime,9)/(outputData(i/cortime,9)+outputData(i/cortime,10));%%GDPedge
        %Want to have where the proteins bind noted
        outputData(i/cortime,12)=timeX(i);%time
% %         if GTPEvL==1 %%1-25-24-SJG return to once find right EvL script
% %         outputData(i/cortime,13)=GTPLat(i/cortime);
% %         outputData(i/cortime,14)=GTPEdge(i/cortime);
% % 
% %         end

        if (i>1000)
            testingValues(i/cortime,1)=(mean(pfLen)-mean(pfLenI(i-1000)))/(timeX(i)-timeX(i-1000));
            testingValues(i/cortime,2)=timeX(i)-timeX(i-1000);
        end
        proteinPlaces{i/cortime,1}=[conditionMatrixadjusted2(EB1BindingPositions,1,1,1,6,pfLen);conditionMatrixadjusted2(EB1BindingPositions,0,1,1,6,pfLen);...
            conditionMatrixadjusted2(EB1BindingPositions,2,1,1,6,pfLen);conditionMatrixadjusted2(EB1BindingPositions,3,1,1,6,pfLen);...
            conditionMatrixadjusted2(EB1BindingPositions,1,0,1,6,pfLen);conditionMatrixadjusted2(EB1BindingPositions,2,0,1,6,pfLen);...
            conditionMatrixadjusted2(EB1BindingPositions,3,0,1,6,pfLen);conditionMatrixadjusted2(EB1BindingPositions,0,0,1,6,pfLen);...
            conditionMatrixadjusted2(EB1BindingPositions,1,2,1,6,pfLen);conditionMatrixadjusted2(EB1BindingPositions,0,2,1,6,pfLen);
            conditionMatrixadjusted2(EB1BindingPositions,2,2,1,6,pfLen);conditionMatrixadjusted2(EB1BindingPositions,3,2,1,6,pfLen)];
        proteinPlacesWithEvL{i/cortime,1}=[presentProteins(:,1:3)];
        protosLengths(i/cortime,:)=pfLen;
    end



    %%%ADDED by SJG on 8-8-22; this will control the lateral bond breaking
    %%%with extreme tapering. Essentially, if the shortest two PFs are
    %%%greater than a given value, then the lateral bonds will break very
    %%%quickly leading to MT depolymerization. The current value of pibreak
    %%%of 0.1 was arbitrarily picked; originally tried a value of 1 but
    %%%this didn't depolymerize quick enough (because lateral bonds didn't
    %%%break fast enough).
    smallestpfLens=mink(pfLen,2);
    if max(pfLen)-smallestpfLens(1)>LongestTaperLength && max(pfLen)-smallestpfLens(2)>LongestTaperLength
        LateralBreakingFoldFactor=1000;%%will essentially divide all lateral bond break rates by this value when conditions are met.  was 1000 to quickly make all lat bonds break.
        kshorten_TBelow = 2; %1/sec %originally was 0.02
        kshorten_DBelow = 2000; %1/sec %originally was 20
    else
        LateralBreakingFoldFactor=1;%this is just resetting the pibreak to the pibreak from the startr of the simulation.
        kshorten_TBelow = 0.2; %1/sec %originally was 0.02, for paper, was 0.2SJG
        kshorten_DBelow = 200; %1/sec %originally was 20, for paper, was 200 .SJG
    end
    if splayingOccurs==1
        for splaying=1:size(pfLen,2)-1 %of note, highest lat bond is to right of point, so pf 1 highest bond is between pfs 1 and 2 (which aligns with EB1 binding position x=2 since EB1 bindign position starts at left edge of PF1 with #1 and goes to right edge of PF13 with #14
            if pfLen(splaying)<pfLen(splaying+1)
                if pfLen(splaying+1)>highestFullLatMade(splaying)
                    EB1BindingPositions(splaying+1,pfLen(splaying+1),1)=1;%Okay, so in the event where PFlen>the lateral bond, the very last dimer will be a single edge and therefore a value of 1.
                end
                EB1BindingPositions(splaying+1,highestFullLatMade(splaying):pfLen(splaying+1)-1,1)=2;%So, the binding pocket above the hightest full lat made will always be an edge because the next dimers won't be bound to make a lattice. Also, it won't be a single edge as there are tubulin dimers in both nearby PFs
                %so if the right PF is longer, this will do everything based off of
                %that
                EB1BindingPositions(splaying+1,pfLen(splaying+1)+1:pfLen(splaying+1)+10,1)=0;%Okay, so this is to clean up past the PF and set those to zero which was not currently happening...see if this solves the issues...

            else
                if pfLen(splaying)>highestFullLatMade(splaying)
                    EB1BindingPositions(splaying+1,pfLen(splaying),1)=1;%Okay, so in the event where PFlen>the lateral bond, the very last dimer will be a single edge and therefore a value of 1.
                end
                EB1BindingPositions(splaying+1,highestFullLatMade(splaying):pfLen(splaying)-1,1)=2;%So, the binding pocket above the hightest full lat made will always be an edge because the next dimers won't be bound to make a lattice. Also, it won't be a single edge as there are tubulin dimers in both nearby PFs
                %if the left PF is longer, then will use it to determine edges...
                EB1BindingPositions(splaying+1,pfLen(splaying)+1:pfLen(splaying)+10,1)=0;%Okay, so this is to clean up past the PF and set those to zero which was not currently happening...see if this solves the issues...
            end
            EB1BindingPositions(splaying+1,highestFullLatMade(splaying)-1:-1:startL,1)=0;%Everything below the highest full lat will for certain be a lattice site.

        end
        %Need PF 13..
        EB1BindingPositions(size(pfLen,2)+1,startL:highestFullLatMade(size(pfLen,2))-1,1)=0;%essentially, all up to highest full lat are zero....
        EB1BindingPositions(size(pfLen,2)+1,pfLen(size(pfLen,2)),1)=1;%tip will always be 1....
        EB1BindingPositions(size(pfLen,2)+1,highestFullLatMade(size(pfLen,2)):pfLen(size(pfLen,2))-1,1)=2;%making after highest full lat and before pfLen as an edge...
        %Below is PF 1 on the left...
        EB1BindingPositions(1,startL:highestFullLatMade(size(pfLen,2))+1,1)=0;%since there is 1.5 dimer difference, dimer 1 on pf 13 blocks dimer 1, 2, and halfway up dimer 3 (but end of dimer 3 still available for binding).
        EB1BindingPositions(1,highestFullLatMade(size(pfLen,2))+2:pfLen(1)-1,1)=2; %since there is splaying, after the fullest lat bond everything but the last point will be an edge.
        if pfLen(1)>highestFullLatMade(size(pfLen,2))+2
            EB1BindingPositions(1,pfLen(1),1)=1; %So, if the PF is longer than the highest lat bond, then the very last point should be a single edge.
        end
        %%%Added by SJG for SPlaying on 1-17-23 %This appears to be working as of
        %%%end of day on 1-17-23. Just need to have EB fall off when bound at
        %%%lattice site on seam.

        %%Okay, for debugging,will add a stop if sum of EB bound below lat is
        %%greater than 0.

    end

end
chosenEventList=TrackSmallest;
shortestFilLength(:)=min(pfLenI,[],2); %shortest protofilament length at each iteration
longestFilLength(:)=max(pfLenI,[],2);
extension_length = longestFilLength - shortestFilLength;
avgFilLength(:)=mean(pfLenI,2);
stdev(:)=std(pfLenI,0,2);
if(plotBool)
    figure(1)
    plot(timeX(:)./60,(avgFilLength(:)-startL)*8,'Color','g'); %subunits 8nm long: time (min) vs filament length (nm)
    xlabel('Time (min)')
    ylabel('MT Length (nm)')
    figure(2)
    plot(timeX(:)./60,extension_length*8,'Color','g');
    xlabel('Time (min)')
    ylabel('Plus-End Extension Length (nm)')
    figure(3)
    histogram(extension_length*8)
    xlabel('Extension Length (nm)')
    ylabel('Count')

    figure(4)
    plot(timeX(:)./60,stdev*8,'Color','m');
    hold on
    plot(timeX(:)./60,(longestFilLength(:)-startL)*8,'Color','b'); %subunits 8nm long: time (min) vs filament length (nm)
    hold off
    xlabel('Time (min)')
    ylabel('Length (nm)')
    legend('Tip standard deviation', 'MT length')

    figure(5)
    plot(timeX(:)./60,extension_length*8,'Color','g');
    hold on
    plot(timeX(:)./60,(longestFilLength(:)-startL)*8,'Color','b'); %subunits 8nm long: time (min) vs filament length (nm)
    hold off
    xlabel('Time (min)')
    ylabel('Plus-End Extension Length (nm)')
    legend('extension length', 'MT length')

    figure(6)
    for j=1:filN
        tubAProps(j,pfLen(j)+1:end,1)=3;
    end
    %format numbers in 14-row matrix for plotting
    tubAProps2=zeros(14,length(tubAProps(j,:,1)));
    tubAProps2(1:filN,:)=tubAProps(:,:,1);
    tubAProps2(filN+1,:)=4*ones(1,length(tubAProps(1,:,1)));
    pcolor(tubAProps2(:,:,1))

    xlabel('Tubulin Subunits')
    ylabel('Protofilament')
    title('GTP Hydrolysis')
    ylim([0 filN+1])
    xlim([0 max(pfLen)])

    hold on
    plot(highestFullLatMade,[1.5:filN+0.5],'*r','MarkerSize',10)
    plot(pfLenI(end,:),[1.5:filN+0.5],'ow','MarkerSize',10)
    hold off

    index=presentProteins(:,3)==1;
    totalEdgeProteins=presentProteins(index,:);
    index=presentProteins(:,3)==2;
    totalLatticeProteins=presentProteins(index,:);
    if size(totalEdgeProteins)>0
        plot(totalEdgeProteins(:,2)+1,totalEdgeProteins(:,1),'*g','MarkerSize',10);%Sam adjusted what was plotted to a variable i created for TPX2 locations (currently) which can be changed to EB1 once the settings for EB1 are created
    end
    if size(totalLatticeProteins)>0
        plot(totalLatticeProteins(:,2)+1,totalLatticeProteins(:,1),'*m','MarkerSize',10);
    end
end
%ADDED by SJG on 8-8-22 copied from Rebecca's Verson 10 code.
if movie_on == 1;

    close(aviobj);
end
proteinTracking(:,7)=0;
for q=1:size(proteinTracking,1)
    proteinTracking(q,7)=proteinTracking(q,5)-proteinTracking(q,4);
end
ii=proteinTracking(:,5)==0&proteinTracking(:,3)==1;%this is for still bound that bound to edge
proteinTracking(ii,10)=13;

jj=proteinTracking(:,5)==0&proteinTracking(:,3)==2;%this is for bound to lat that are still bound
proteinTracking(jj,10)=14;%This is for still bound that bound to lattice
toc
end