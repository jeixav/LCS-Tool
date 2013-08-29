function position = impose_periodic_bc(position,domain,periodicBcSpec)

for i = 1:size(periodicBcSpec,2)
    if periodicBcSpec(i)
        length = domain(i,2) - domain(i,1);
        position(i,:) = mod(position(i,:) - domain(i,1),length) + domain(i,1);
    end
end

end

