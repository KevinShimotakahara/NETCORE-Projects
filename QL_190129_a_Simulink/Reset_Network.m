
function Reset_Network()

    %% Variables
    global eNBs
    global SBSs
    global UDs
    
    clear eNBs
    clear SBSs
    clear UDs
    % ******************************************************************* %
    
    %% Reset network parameters
    global Seed;
    
    Seed = 2 * Seed;
    Main_Init();
    % ******************************************************************* %

end



