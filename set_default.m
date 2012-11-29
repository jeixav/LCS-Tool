function input = set_default(input,default)

% default = struct('graphs',true,'progress',true,'stats',true);

lFieldnames = fieldnames(default);
nFields = numel(lFieldnames);

for iField = 1:nFields
    if ~isfield(input,lFieldnames{iField})
        input.(lFieldnames{iField}) = ...
            default.(lFieldnames{iField});
    end
end
