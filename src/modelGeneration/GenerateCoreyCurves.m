function [ SW_corey_m,krw_model,kro_model ] = ...
    GenerateCoreyCurves( NbSimu, Swir, Swor,krw_end,kro_end,no,nw)
%GENERATECOREYCURVES Summary of this function goes here
%   Detailed explanation goes here

RelPermEntries = 25;
SW_corey_m = zeros(NbSimu,RelPermEntries);
krw_model = zeros(NbSimu,RelPermEntries);
kro_model = zeros(NbSimu,RelPermEntries);

for i=1:NbSimu
    SW_corey_m(i,:) = linspace(Swir(i), ...
        1 - Swor(i),RelPermEntries);
end
 
for i=1:NbSimu
     krw_model(i,:) = krw_end(i) .* ...
         ((SW_corey_m(i,:)-Swir(i))./...
         (1-Swir(i)-...
         Swor(i))).^nw(i);
     kro_model(i,:) = kro_end(i) .* ...
         ((1 - SW_corey_m(i,:) - Swor(i))./...
         (1 - Swir(i) - Swor(i))).^...
         no(i);
     
     plot(SW_corey_m(i,:), krw_model(i,:),'b-', SW_corey_m(i,:), ...
         kro_model(i,:), 'r-');
     xlabel('S_w');
     ylabel('Relative Permeability');
     xlim([0 1]); ylim([0 1]); grid on; hold on;
end

end
