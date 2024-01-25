%The goal is to make a function that will subset a matrix based on
%conditions that I want. For this, you would pass a matrix that is
%positions and each position has a given condition and it will return the
%positions that have the desired conditions. The assumption here is that
%the matrix is 3D matrix where dimension 1 is X, dimension 2 is y, and
%dimension 3 is conditions which demonstrate if lattice/edge bound, dimer
%hydrolysis state, and eb1 bound for each of the three positions in the
%third dimension. 
function[matrix]=conditionMatrixadjusted(StartMatrix,latticeVsEdge,DimerHydroState,EB1Bound,seedSize,pfLenVector);
dimensions=size(StartMatrix);
matrix=[];
counter=1;
for i=1:dimensions(1)
    for j=seedSize:pfLenVector(i+1)
        if StartMatrix(i,j,1)==latticeVsEdge && StartMatrix(i,j,2)==DimerHydroState && StartMatrix(i,j,3)==EB1Bound
            matrix(counter,1)=i+1;
            matrix(counter,2)=j;
            counter=counter+1;
        else
            counter=counter+0;
        end
    end
end
end