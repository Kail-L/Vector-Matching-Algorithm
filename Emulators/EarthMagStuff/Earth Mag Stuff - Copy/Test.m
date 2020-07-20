
load('SatPosAtt.mat')

t = RefVecs(:,1);
r = RefVecs(:,2:4);

BiCaverly = zeros(3,length(t));
BiCalc = zeros(3,length(t));

for i = 1:length(t)
    BiCaverly(:,i) = EarthMagField(r(i,:)',t(i));
    BiCalc(:,i) = EarthMag3rdOrder(r(i,:),t(i));
    
end

%% Plot
plot3(BiCaverly(1,:),BiCaverly(2,:),BiCaverly(3,:));
hold on;
grid on;
plot3(BiCalc(1,:),BiCalc(2,:),BiCalc(3,:));
legend('Caverly solution', 'Our solution');