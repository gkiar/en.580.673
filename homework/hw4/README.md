##Homework 4

### Part 1: Image Reconstruction

1 **Reconstructing the Data**
  - [**``Reconstruct Function``**](./code/Reconstruct.m)
  - [**``Reconstructed Images``**](./results/part1/scans/)
  
2 **Measure SNR**
  - [**``Location of ROIs``**](./results/part1/rois.png)
  - SNR Table:

| Scan | Signal Mean | Noise STD |   SNR  |
|------|-------------|-----------|--------|
| 6    | 5.0967e+05  | 1.5762e+03|323.3520|
| 7    |5.3141e+05   | 1.5694e+03|338.6021|
| 8    |5.3959e+05   | 1.9027e+03|283.5944|
| 12   |5.4627e+04   | 1.9693e+03| 27.7393|
| 13   |5.4677e+04   | 1.9168e+03| 28.5246|
| 14   |5.5428e+04   | 1.5864e+03| 34.9400|
| 15   |5.7323e+04   | 2.4103e+03| 23.7827|
| 17   |2.2688e+05   | 3.7176e+03| 61.0292|
| 18   |6.0563e+04   | 1.8177e+03| 33.3188|
| 19   |3.0125e+04   | 1.0083e+03| 29.8768|
| 20   |5.9963e+04   | 3.8026e+03| 15.7689|
| 21   |5.8248e+04   | 3.8738e+03| 15.0365|
| 25   |   | | |

3 **Water Content**
  - The body of the lobster, the center of the coconut and pretty much the whole fish (sans eye) contain a lot of water.
  - The tissues which had the highest difference in proton density seemed to be the fleshy part of the lobster and it's shell. The shell contained almost no signal, but the flesh contained lots.
  - In my opinion, the PD sequence was best for lobster, the T2W signal was best for the fish, and the T1W was best for the coconut. The T1W and T2W signals were quite similar here, but I noticed subtl differences in favour of the selections made above.

4 **Big ass image**
  - The results for this are above (i.e. image with other images, and SNR in previous table).
  
### Part 2: The Noise in MRI

1 **Visualzing Noise Data**
  - [**``Phasor Images for Each Coil``**](./results/part2/noise_phasors.png)
  - The noise appears to me as approximately Gaussian. The mean is near zero (though slightly higher) which makes sense to me as we expected 0 centered noise but there are non idealities in the gradients which I imagine increase the noise baseline. The standard deviation value seems relatively uninformative to me here, as the images are on somewhat of an arbitrary intensity scale determined by the gain of the scanner. One important difference I did notice is that noise is non-uniform across the coils. I expect this is due to the location of the coils within the scanner, and heterogeneity of the B0 field causes the coils to react differently in the central, ideal, region of the scanner and the periphery, where that ideality breaks down.
  - We also notice that the signal is discrete, due to the hardware limitations of the ADC.

2 **Decomposing Noise Distributions** 
  - [**``Noise Historams for Each Coil``**](./results/part2/noise_hists.png)
  - Again, the noise appears to be Gaussian to me. They do not all have the same amplitude, but as they all contain the same number of samples (i.e. the vectors are the same length), that means that those with lower amplitude have a higher standard deviation of their noise. I again think that this difference is due to the location of the coils within the scanner, and that they experience different magnetic fields based on their location which leads to different noise characteristics.

3 **Image Space Noise**
  - [**``Noise Histograms for Each Coil``**](./results/part2/noise_fft_hists.png)
  - The distribution does still look Guassian, though of a different amplitude. The scaling factor inherent to the Fourier Transform is 1/sqrt(length(noise_data)), so when applying this we obtain a more similarly scaled distribution.

4 **Complex Histograms after FFT**
  - [**``Magnitude Noise images``**](./results/part2/magnitude_hists.png)
  - Magnitude of noise looks Raleigh/Rician, not Gaussian, when we take the magnitude.

5 **Imaginary and Real Noise Relationship**
  - We can predict parameters of the Rician distribution in image space based on parameters of the Gaussian distribution in k-space. The relationship is given by equation 3 [here](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC2254141/).

6 **Distribution of noise**
  - Using the parameters and equations we didn't seem to get values that very closely matched those observed - it seems that the equations over estimated the noise.

|Observed Mean | Observed STD | Calc'd Mean | Calc'd STD |
|--------------|--------------|-------------|------------|
|    0.0458    | 0.0254       |     6.3491  | 11.0145    |
|    0.0564    | 0.0314       |     7.8277  | 16.7420    |
|    0.0498    | 0.0276       |     6.9068  | 13.0345    |
|    0.0636    | 0.0354       |     8.8270  | 21.2895    |
|    0.0666    | 0.0371       |     9.2428  | 23.3429    |
|    0.0501    | 0.0278       |     6.9485  | 13.1924    |


### Part 3: The Signal in MRI

1 **Averaging for SNR**
  - This somewhat closely follows the pattern we expected. Moving from 1 average to 2, we expect a sqrt(2) times increase, which would bring SNR up to 45, slighly above the measured SNR. Then from 2 to 3 we expect a sqrt(3)/sqrt(2) increase, which would bring us to 51, our exact measurement. Lastly, from 3 to 4 would be a swrt(4)/sqrt(3) gain, bringing us to 60, a tiny bit higher than that observed.
  - SNR Table:

|Avgs|    Signal  |   Noise    |   SNR   |
|----|------------|------------|---------|
| 1  | 5.3800e+04 | 1.6711e+03 | 32.1934 |
| 2  | 1.1033e+05 | 2.6214e+03 | 42.0891 |
| 3  | 1.6251e+05 | 3.1725e+03 | 51.2242 |
| 4  | 2.1834e+05 | 3.7476e+03 | 58.2629 |


2 **Double oversampling**
  - First, we define SNR \proportional dx dy dz \sqrt( Nx Ny Nz Nav dt ). Allow us to now walk through the steps given. When you double the FOV in the x direction, you are changing the SNR by 2, as dx = FOVx/Nx. Then, when you double Nx, you multiply SNR by \sqrt(2)/2, leaving you now \sqrt(2) up from where you were initially. Finally, doubling the bandwidth is the same has halfing dt, so you are now multiplying by 1/\sqrt(2), leaving you with no change in your SNR. The advantage of oversampling then becomes obvious that while you do not lose out on SNR, you are able to double the number of samples in the readout direction, and remove aliasing along the x direction.
  - As we're now covering twice as much of k-space in the x direction, the area under our readout gradient has to double. As you're doubling Nx, and dt is being halved, your readout gradient will take the same amount of time. This means that in order to get our larger field of view, then, the amplitude of the readout gradient must double.

3 **Aliasing in PE and FE**
  - This is because we double oversample in the readout direction, so when we reduce the FOV we preserve the image region. We cannot double oversample in the phase encode direction, so we notice ghosting/aliasing/wrap around.

4 **Resolution and SNR**
  - The SNR relationship does not seem to hold in these images. I think it is because there is a "lower bound" on SNR of the images, and as the SNR is sufficiently low/close to this posited "lower bound" here, single changes will have a minimal effect. I imagine if the images were of higher quality, the SNR relationship would hold when a similar change occured.
  - In the case of image 25, the relationship does more closely hold as the image is of way higher quality. This backs up our theory about the previous part, as well.

| Scan | Signal Mean | Noise STD |   SNR  |
|------|-------------|-----------|--------|
| 20   |5.9963e+04   | 3.8026e+03| 15.7689|
| 21   |5.8248e+04   | 3.8738e+03| 15.0365|
| 25   |             |           |        |


5 **Bandwidth and SNR**
  - Image bandwidths were calculated with the equation: BW1/ps1 = BW2/ps2. The image bandwidth is the pixel bandwidth times the number of pixels.
  - We notice that the relationship doesn't hold like we would expect. I maintain my theory of crappy images breaking down the relationship.

| Scan | Pixel Shift | Pixel Bandwidth |  Image Bandwidth  |     dt   |   SNR   |
|------|-------------|-----------------|-------------------|----------|---------|
| 12   |    1        |  434 Hz         |   156240   Hz     | 6.40e-6  | 27.7393 |
| 13   |    2        |  868 Hz         |   312480   Hz     | 3.20e-6  | 28.5246 |
| 14   |    1.5      |  651 Hz         |   234360   Hz     | 4.27e-6  | 34.9400 |
| 15   |    0.5      |  217 Hz         |    78120   Hz     | 1.28e-5  | 23.7827 |
