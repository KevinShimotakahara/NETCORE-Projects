function appLayer(u)
    switch u(1)
        case 1
            runDevice1App(u(2),u(3))
        case 2
            runDevice2App(u(2),u(3))
        case 3
            runDevice3App(u(2),u(3))
        case 4
            runDevice4App(u(2),u(3))
        case 5
            runDevice5App(u(2))
        case 7
            runDevice7App(u(2))
        otherwise
    end
end