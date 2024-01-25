%The goal is to make a function that will subset a matrix based on
%conditions that I want. For this, you would pass a matrix that is
%positions and each position has a given condition and it will return the
%positions that have the desired conditions. The assumption here is that
%the matrix is 3D matrix where dimension 1 is X, dimension 2 is y, and
%dimension 3 is conditions which demonstrate if lattice/edge bound, dimer
%hydrolysis state, and eb1 bound for each of the three positions in the
%third dimension. 
function[matrix]=conditionMatrixadjusted2(StartMatrix,latticeVsEdge,DimerHydroState,EB1Bound,seedSize,pfLenVector);
dimensions=size(StartMatrix);
matrix=[];
counter=1;

for p=1:dimensions(1)
    if p==1 %When the PF is on either end, you only need to check that PF to determine where binding spots are
        if p==1
            yValue=pfLenVector(p);%This is the length of PF1 for edge 1
        else
            yValue=pfLenVector(p-1);%This gives the length of PF13 for edge 14.
        end
        for j=seedSize:yValue
            if StartMatrix(p,j,1)==latticeVsEdge && StartMatrix(p,j,2)==DimerHydroState && StartMatrix(p,j,3)==EB1Bound
                matrix(counter,1)=p;
                matrix(counter,2)=j;
                counter=counter+1;
            else
                counter=counter+0;
            end
        end
    elseif p==2 || p==3 || p==4 || p==5 || p==6 || p==7 || p==8 || p==9 || p==10 || p==11 || p==12 || p==13
        a=pfLenVector(p);
        b=pfLenVector(p-1);
        yValue=max([a,b]);%for middle PFs, need to check the dimer to the left and right to determine how far y to check for edges
        for j=seedSize:yValue
            if StartMatrix(p,j,1)==latticeVsEdge && StartMatrix(p,j,2)==DimerHydroState && StartMatrix(p,j,3)==EB1Bound
                matrix(counter,1)=p;
                matrix(counter,2)=j;
                counter=counter+1;
            else
                counter=counter+0;
            end
        end
    else
        if p==dimensions(1)%When the PF is on either end, you only need to check that PF to determine where binding spots are
            if p==1
                yValue=pfLenVector(p);%This is the length of PF1 for edge 1
            else
                yValue=pfLenVector(p-1);%This gives the length of PF13 for edge 14.
            end
            for j=seedSize:yValue
                if StartMatrix(p,j,1)==latticeVsEdge && StartMatrix(p,j,2)==DimerHydroState && StartMatrix(p,j,3)==EB1Bound
                    matrix(counter,1)=p;
                    matrix(counter,2)=j;
                    counter=counter+1;
                else
                    counter=counter+0;
                end
            end
        end
    end
end    
end
