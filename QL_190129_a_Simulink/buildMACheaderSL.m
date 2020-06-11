%% MAC
function MACheaderSL = buildMACheaderSL(SourceID,DestID,L,LCID,CE,varargin) %may need TBS for padding subheader
    if CE
        %This is a sidelink BSR message
        MACheaderSL = [0 0 0 LCID];
        %SL BSR CE
        DestIndex = num2binHeaderComponent(varargin{1},4);
        LGCID = num2binHeaderComponent(varargin{2},2);
        buffSize = num2binHeaderComponent(varargin{3},6);
        SLbsrCE = [DestIndex LGCID buffSize 0 0 0 0];
        MACheaderSL = [MACheaderSL SLbsrCE];
    else
        %This is sidelink data transmission 

        %Build header components
        SourceIDheaderComponent = num2binHeaderComponent(SourceID,24);
        DestIDheaderComponent = num2binHeaderComponent(DestID,16);
        %useless bits reserved for forward compatibility:
        useless = zeros(1,8);
        if L < 2^7
            %size of L field should be 7
            LheaderComponent = num2binHeaderComponent(L,7);
            MACheaderSL = [useless SourceIDheaderComponent DestIDheaderComponent 0 0 0 LCID 0 LheaderComponent];
        elseif L < 2^15
            %size of L field should be 15
            LheaderComponent = num2binHeaderComponent(L,15);
            MACheaderSL = [useless SourceIDheaderComponent DestIDheaderComponent 0 0 0 LCID 1 LheaderComponent];
        elseif L < 2^16
            %size of L field should be 16
            LheaderComponent = num2binHeaderComponent(L,16);
            MACheaderSL = [useless SourceIDheaderComponent DestIDheaderComponent 0 1 0 LCID LheaderComponent];
        else
            fprintf("ERROR: L exceeds what can be expressed in header.\n")
        end

    end
end

function binHeaderComponent = num2binHeaderComponent(number,lengthHeaderComponent)
    %This funcition automates the task of converting a decimal value into a
    %binary value represented by an array of zeros and ones, then sticking
    %it into a header of size lengthHeaderComponent, adding any padding to
    %the MSBs of the header if needed.
    numberchar = dec2bin(number);

    paddingLength = lengthHeaderComponent - length(numberchar);
    binHeaderComponent = zeros(1,lengthHeaderComponent);

    for k = paddingLength+1:lengthHeaderComponent
        binHeaderComponent = str2double(numberchar(k-paddingLength));
    end
end