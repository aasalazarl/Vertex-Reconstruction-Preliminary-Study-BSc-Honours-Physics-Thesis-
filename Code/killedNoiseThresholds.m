function killedNoiseThresholds

threshold = [750 200 55 55]
numberOfSiPMs = [1 2 3 4]

unc_threshold = [50 50 5 5] 

realThreshold = threshold./50 % The signal was being amplified by x5 x10.
unc_realThreshold = unc_threshold./50

errorbar(numberOfSiPMs, realThreshold, unc_realThreshold, '.')
%title('Noise killing by number of SiPMs in coincidence')
xlabel('Number of SiPMs')
ylabel('Threshold [mV]')
hold on

XX = [0.93 1 1.25 1.5 1.75 2 2.25 2.5 2.75 3 3.25 3.5 3.75 4 4.07]
YY = spline(numberOfSiPMs,realThreshold, XX)
plot(XX, YY)
legend('Data points', 'Interpolation')

print(gcf, 'killingThreshold', '-dpng', '-r300')

end