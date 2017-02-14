function ratio = neighbourhoodAreaRatio(AF, AFNew, TRIV, nVerts)
    I = [TRIV(:,1);TRIV(:,2);TRIV(:,3)];
    J = ones(length(I),1);
    S = [AF;AF;AF];
    AV = sparse(I,J,S,nVerts,1);
    
    S = [AFNew;AFNew;AFNew];
    AVNew = sparse(I,J,S,nVerts,1);
    
    ratio = full(AV)./full(AVNew);
end