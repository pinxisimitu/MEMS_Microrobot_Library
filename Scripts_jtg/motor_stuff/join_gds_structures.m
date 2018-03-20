% combine the elements of gsStruct1 and gsStruct2
% elements from gsStruct2 are taken and added to gsStruct1
function gs=join_gds_structures(gsStruct1,gsStruct2)
    for i=1:numel(gsStruct2)
        gsStruct1(end+1)=gsStruct2(i);
    end
    gs=gsStruct1;
end