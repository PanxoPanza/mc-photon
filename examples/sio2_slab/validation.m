clear all
addpath matlab_exta_library\
%--------------------------------------------------------------------------
% compute exact results
%--------------------------------------------------------------------------
tfilm = 100;                % film thickness in um
ll = linspace(0.3,1.0,100);   % wavelength range in um
N = sqrt(2.25);             % refractive index

[Rf, Tf] = Fresnel(0,N,1.0,'p');

R = Rf + Tf.^2.*Rf.*exp(-8*pi*imag(N)*tfilm./ll)./...
        (1 - Rf.^2.*exp(-8*pi*imag(N)*tfilm./ll));

T = Tf.^2.*exp(-4*pi*imag(N)*tfilm./ll)./...
    (1 - Rf.^2.*exp(-8*pi*imag(N)*tfilm./ll));

figure, hold on
% plot exact values
plot(ll,R,'-r','DisplayName', 'R_{tot} (exact)')    % total reflectance
plot(ll,T,'-b','DisplayName', 'T_{tot} (exact)')    % total transmittance
hold off

%--------------------------------------------------------------------------
%           Plot mc-photon results
%--------------------------------------------------------------------------
% extract data from mc-photon
dataRtot = ReadFile('Rtot_glass_slab.pow');     % Total reflectance
dataTtot = ReadFile('Ttot_glass_slab.pow');     % Total transmittance
dataTspec = ReadFile('Tspec_glass_slab.pow');   % Specular transmittance

% plot results from mc-photon
hold on
plot(dataRtot(:,3),dataRtot(:,4),'o r','DisplayName', ...
    'R_{tot} (mc-photon)','MarkerSize',3.0,'MarkerFaceColor','r')
plot(dataTtot(:,3),dataTtot(:,4),'o b','DisplayName', ...
    'T_{tot} (mc-photon)','MarkerSize',3.0,'MarkerFaceColor','b')
plot(dataTspec(:,3),dataTspec(:,4),'o b','DisplayName', ...
    'T_{spec} (mc-photon)','MarkerSize',3.0,'MarkerFaceColor','w')
hold off

%--------------------------------------------------------------------------
% Format plot labels and axis
%--------------------------------------------------------------------------
% format legend
lh = legend('show');
lh.Box = 'off';
lh.FontSize = 14;

% format figure axis
ax = gca;
ax.FontSize = 14;
ax.FontName = 'Calibri Light';
ax.XLabel.String = 'Wavelength (\mum)';
ax.YLabel.String = 'Reflectance / Transmittance';
ax.YLim = [0 1];
ax.Title.String = 'Glass slab (0.1 mm thick)';

rmpath matlab_exta_library\ % remove path