%--------------------------------------------------------------------------
% Testing a multilayered film at 0.3 - 1.0 um wavelength range. The
% multilayered film consist of:
%   
%   50 nm   silver
%   500 nm  dielectric (N = 1.5)
%   200 nm  silver
%
% The results consider the reflectance and transmittance computed from
% transfer matrix method.
%
% F. Ramirez 2022
%--------------------------------------------------------------------------

data = ReadFile('exact_solution.txt');
lambda = data(:,1);
R = data(:,2);
T = data(:,3);

figure,
hold on
plot(lambda,R*100,'r','DisplaYName','R')
plot(lambda,T*100,'b','DisplaYName','T')
hold off

ax = gca;
ax.FontName = 'Calibri';
ax.FontSize = 14;
ax.XLim = [0.3,1.0];
ax.YLim = [0 100];
ax.XLabel.String = 'Wavelength (\mum)';
ax.YLabel.String = 'Reflectance/Transmittance (%)';
legend('show');
ax.Legend.Box = 'off';

data = ReadFile('Rt_multilayer_test.pow');
lam0 = data(:,3); Rsim = data(:,4)*100;
data = ReadFile('Tt_multilayer_test.pow');
Tsim = data(:,4)*100;

hold on,
plot(lam0,Rsim,'o r','DisplayName','R (mc-photon)',...
        'MarkerSize',3.0,'MarkerFaceColor','r')
plot(lam0,Tsim,'o b','DisplayName','T (mc-photon)',...
        'MarkerSize',3.0,'MarkerFaceColor','b')
hold off