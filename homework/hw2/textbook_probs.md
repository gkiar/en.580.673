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
a) Using a double-sided 2DFT sequence on a signal with a slight frequency deviation from the expected frequency, f\_0 + \Delta\_f, the FFT of the impulse (i.e the measured signal by the scanner) would have a linear frequency gradient along k\_x. Thus, when the transform is taken back and the image reconstructed, the delta function would appear shifted in x direction by distance of 1/\Delta\_f.

b) If a 1-sided sequence is used, then the projections are taken radially from the origin in k-space. Therefore, the gradient which appears equally in all slices of the double sided sequence now appears inconsistently along each slice. The signal for each slice, s(t), will effectively be multiplied by an exponential term, exp(-i\*kx\*\Delta\_f), which will make points with higher k\_x values more affected. Since the projection sequence is symmetric, the recovered image will be effectively just blurred.

#### 5.6
a)

b)

c)
#### 5.7

