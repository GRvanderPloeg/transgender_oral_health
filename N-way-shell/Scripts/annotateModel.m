function annotatedModel=annotateModel(X, model, metadata)
%
% annotateModel: orders and adds metadata to an existing PARAFAC model.
% 
% See also:
% 'parafac', 'sortParafacComponents', 'calcVarExplained'
%
% annotatedModel=annotateModel(X, model, metadata)
%
% Orders and adds metadata to an existing PARAFAC model.
% Uses the raw input data to determine the variance explained of each
% PARAFAC component and then orders then accordingly.
% Attaches metadata to the loading matrix.
%
% ---------------------- INPUT ---------------------
%
% X                 This is the input array and needs to have the same
%                   number of modes as the PARAFAC model.
%
% model             The PARAFAC model as output of the parafac() function.
%                   Optionally you can make a similar object by creating
%                   a cell object and filling each element with the 
%                   loadings of one mode at a time.
%
% metadata          Metadata in a cell object, with the metadata 
%                   corresponding to one mode in each cell element.
%
% ---------------------- OUTPUT ---------------------
%
% annotatedModel    The PARAFAC model in a cell object, with each cell
%                   being a mode. The metadata is attached to the right
%                   side of every loading matrix. This causes the loadings
%                   of an N component model to still exist in the first N
%                   columns of each loading matrix.

% Created by: G.R. van der Ploeg (g.r.ploeg@uva.nl)

orderedModel = sortParafacComponents(X, model);
numModes = size(orderedModel, 2);
annotatedModel = {};

for i=1:numModes
    annotatedModel{i} = [orderedModel{i} metadata{i}];
end