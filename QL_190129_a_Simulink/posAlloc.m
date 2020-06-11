
function [Pos] = posAlloc(enbRadius, NodePos)

    % Generate (x, y) point
    xPos    = abs(normrnd((enbRadius/4), (enbRadius/4)));
    yPos    = abs(normrnd((enbRadius/4), (enbRadius/4)));
    xSign   = 2 * randi([0, 1], 1, 1) - 1;
    ySign   = 2 * randi([0, 1], 1, 1) - 1;
    xPos    = xSign * xPos;
    yPos    = ySign * yPos;

    xPos = xPos + NodePos(1);
    yPos = yPos + NodePos(2);

    Pos = [xPos, yPos];
    
end



