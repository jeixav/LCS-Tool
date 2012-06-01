function position = impose_periodic_bc(position,domain,periodic_bc_spec)

for i = 1:size(periodic_bc_spec,2)
    if periodic_bc_spec(i)
        length = domain(i,2) - domain(i,1);
        position(:,i,:) = mod(position(:,i,:) + .5*length,length) - .5*length;
    end
end

end

