function result = isTriggerType(thing)
    result = isequal(thing,'builtin') || ...
             isequal(thing,'counter') || ...
             isequal(thing,'external') ;
end