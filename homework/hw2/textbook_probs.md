### Answers to textbook problems

#### 4.2
*Todo*: attach figure

#### 4.3
*Todo*: attach figure

#### 4.4
*Todo*: attach figure

#### 5.1
*Todo*: attach figure

#### 5.5
a) Using a double-sided 2DFT sequence on a signal with a slight frequency deviation from the expected frequency, f\_0 + \Delta\_f, the FFT of the impulse (i.e the measured signal by the scanner) would have a linear frequency gradient along k\_x. Thus, when the transform is taken back and the image reconstructed, the delta function would appear shifted in x direction by distance as seen below:

b) If a 1-sided sequence is used, then the projections are taken radially from the origin in k-space. Therefore, the gradient which appears equally in all slices of the double sided sequence now appears inconsistently along each slice. The signal for each slice, s(t), will effectively be multiplied by an exponential term, i.e. s(t)' = s(t)\*exp(-2\*\pi\*x\*\Delta\_f), which will make points with higher x values more affected. Since the projection sequence is symmetric, the recovered image will be effectively just blurred.

#### 5.6
Since the points are separated by distance L, the frequency of the reconstructed sinusoid is L/2. Half of the period is then 1/L.

a) In order for the reconstruction to be 0, the waveform must be sampled at the nodes of the sinusoid. Therefore, if \Delta\_ky = nL and subsequently FOV\_y = 1/nL for any integer n greater than or equal to 1, the points sampled will fall at one of the nodes of the sinusoid.

b) In order for the image to be a single impulse, the same intensity must be observed at every point. Therefore, if \Delta\_ky = 2nL and subsequently FOV\_y = 1/2nL for any integer n greater than or equal to 1, the points sampled will always fall at the same point of the sinusoid.

c) In either case, the shifting wouldn't make a difference in terms of the pattern which is observed. However, in the case presented in b), the magnitude of the impulse will vary based on the magnitude of the shift.


#### 5.7

