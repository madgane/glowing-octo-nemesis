
function userLoc = getPointInRhombus(hexSide,eastRotAngle)

RAD_90 = (pi / 2);
RAD_30 = (pi / 6);
RAD_120 = (2 * pi / 3);

nRhombus = 3;
rhombusSide = hexSide * cos(RAD_30);
userLoc = complex(rand,rand) * rhombusSide;

K = tan(RAD_30);K1 = sec(RAD_30);
skewMat = [1 0 ; -K K1];
userLoc = skewMat * [real(userLoc) ; imag(userLoc)];
userLoc = userLoc(1,1) + sqrt(-1) * userLoc(2,1);

xRhombus = randi(nRhombus,1,1) - 1;
userLoc = userLoc * exp(sqrt(-1) * ((xRhombus * RAD_120) - RAD_90 - eastRotAngle));

% userLoc = userLoc * exp(sqrt(-1) * xRhombus * RAD_120);
% 
% if inclSector
%     userLoc = (userLoc  + exp(-sqrt(-1) * RAD_90) * hexSide) * exp(sqrt(-1) * RAD_120 * (cSector - 1)) * exp(-sqrt(-1) * eastRotAngle);
% else
%     userLoc = userLoc * exp(-sqrt(-1) * eastRotAngle);
% end

end