function input = set_default(input,default)

lFieldnames = fieldnames(default);
nFields = numel(lFieldnames);

for iField = 1:nFields
    if ~isfield(input,lFieldnames{iField})
        input.(lFieldnames{iField}) = ...
            default.(lFieldnames{iField});
    end
end
