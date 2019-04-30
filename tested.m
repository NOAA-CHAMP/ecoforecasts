tested()

    global station;
    global x;

    evalin('base', 'global x; x=station.factories.spawning_seatemp.fuzzies;');

    openvar x

    evalin('base', 'station.factories.spawning_seatemp.fuzzies = x;');

return;
