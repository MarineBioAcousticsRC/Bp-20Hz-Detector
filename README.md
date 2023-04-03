# Bp-20Hz-Detector
Acoustic energy detector for fin whale 20 Hz calls

Calculates the difference in acoustic energy between the signal (22 Hz) and background noise levels (average energy between 10 and 34 Hz) from a long-term spectral average (LTSA) in [Triton](https://github.com/MarineBioAcousticsRC/Triton.git).

Outputs an XML file with the average ratio for each day of the timeseries

## Notes:
1. If you have 200 kHz data, use XMLFin3PowerDetectDay.m and see the [Wiki](https://github.com/alrice12/Bp-20Hz-Detector/wiki) for a tutorial on running the code. This code can handle small gaps in data.
2. If you have 320 kHz data, use XMLFinPowerDetectDaySOCAL_320k.m. This code has **not** been modified to handle data gaps, which can result in days becoming offset. Procedure is the same as above but results should be scrutinized if code has not been modified first.
3. Tested with MATLAB 2013b and 2016b

## References:

1. Širović A, Hildebrand JA, Wiggins SM, McDonald MA, Moore SE, Thiele D (2004) Seasonality of blue and fin whale calls and the influence of sea ice in the Western Antarctic Peninsula. Deep-Sea Res Part II-Top Stud Oceanogr 51:2327–2344. https://doi.org/10.1016/j.dsr2.2004.08.005

2. Širović A, Rice A, Chou E, Hildebrand JA, Wiggins SM, Roch MA (2015) Seven years of blue and fin whale call abundance in the Southern California Bight. Endanger Species Res 28:61–76. https://doi.org/10.3354/esr00676
