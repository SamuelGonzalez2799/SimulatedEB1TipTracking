function result=CheckTags(newX,newY,tagA,len,n,recurse)
xys=[];
success=0;
if(newX==1)
    if(...
            (len(n)+2>len(newX) && isequal(tagA(n,newY-2:newY-1,1),[0,0]) )...
            ||...
            ((len(newX+1)>len(newX)) && isequal(tagA(newX+1,newY:newY+1,1),[0,0]))...
            ||...
            ((len(n)+2==len(newX)&&len(newX+1)==len(newX)) && (isequal(tagA(n,newY-2,1),0)&&(isequal(tagA(newX+1,newY,1),0))))...
            )
        xys=[newX,newY];
        tagA(newX,newY,2)=0;
        i=1;
        success=1;
        while(tagA(newX,newY-i,2)==1)
            xys(end+1,:)=[newX,newY-i];
            i=i+1;
        end
    end
    %%%%
    %these need to check adjacent filaments for tagged tubulin that can be
    %untagged. May or may not need to recurse... not sure.
    if(recurse || success)
        if(tagA(newX+1,newY,2))
            xys=[xys;CheckTags(newX+1,newY,tagA,len,n,0)];
        elseif(tagA(newX+1,newY-1,2))
            xys=[xys;CheckTags(newX+1,newY-1,tagA,len,n,0)];
        end
        if(tagA(n,newY-1,2))
            xys=[xys;CheckTags(n,newY-1,tagA,len,n,0)];
        elseif(tagA(n,newY-2,2))
            xys=[xys;CheckTags(n,newY-2,tagA,len,n,0)];
        end
    end
    
    
    
elseif(newX==n)
    if(...
            ((len(newX-1)>len(newX)) && isequal(tagA(newX-1,newY:newY+1,1),[0,0]))...
            ||...
            ((len(1)-2>len(newX)) && isequal(tagA(1,newY+1:newY+2,1),[0,0]))...
            ||...
            ((len(newX-1)==len(newX)&&len(1)-2==len(newX)) && (isequal(tagA(newX-1,newY,1),0)&&(isequal(tagA(1,newY+1,1),0))))...
            )
        xys=[newX,newY];
        tagA(newX,newY,2)=0;
        i=1;
        success=1;
        while(tagA(newX,newY-i,2)==1)
            xys(end+1,:)=[newX,newY-i];
            i=i+1;
        end
    end
    
    if(recurse || success)
        if(tagA(1,newY+2,2))
            xys=[xys;CheckTags(1,newY+2,tagA,len,n,0)];
        elseif(tagA(1,newY+1,2))
            xys=[xys;CheckTags(1,newY+1,tagA,len,n,0)];
        end
        if(tagA(newX-1,newY,2))
            xys=[xys;CheckTags(newX-1,newY,tagA,len,n,0)];
        elseif(tagA(newX-1,newY-1,2))
            xys=[xys;CheckTags(newX-1,newY-1,tagA,len,n,0)];
        end
    end
    
    
    
    
else
    if(...
            ((len(newX-1)>len(newX)) && isequal(tagA(newX-1,newY:newY+1,1),[0,0]))...
            ||...
            ((len(newX+1)>len(newX)) && isequal(tagA(newX+1,newY:newY+1,1),[0,0]))...
            ||...
            ((len(newX-1)==len(newX)&&len(newX+1)==len(newX)) && isequal(tagA([newX-1,newX+1],newY,1),[0;0]))...
            )
        xys=[newX,newY];
        tagA(newX,newY,2)=0;
        i=1;
        success=1;
        while(tagA(newX,newY-i,2)==1)
            xys(end+1,:)=[newX,newY-i];
            i=i+1;
        end
    end
    
    if(recurse || success)
        if(tagA(newX+1,newY,2))
            xys=[xys;CheckTags(newX+1,newY,tagA,len,n,0)];
        elseif(tagA(newX+1,newY-1,2))
            xys=[xys;CheckTags(newX+1,newY-1,tagA,len,n,0)];
        end
        if(tagA(newX-1,newY,2))
            xys=[xys;CheckTags(newX-1,newY,tagA,len,n,0)];
        elseif(tagA(newX-1,newY-1,2))
            xys=[xys;CheckTags(newX-1,newY-1,tagA,len,n,0)];
        end
    end
    
    
end

result=xys;

end