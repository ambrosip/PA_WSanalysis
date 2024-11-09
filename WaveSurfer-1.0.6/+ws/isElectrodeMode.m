function result = isElectrodeMode(thing)
    result = isequal(thing,'vc') || ...
             isequal(thing,'cc') || ...
             isequal(thing,'i_equals_zero') ;
end