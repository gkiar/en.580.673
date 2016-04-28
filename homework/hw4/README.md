##Homework 4

### Part 1: Image Reconstruction

1. **Reconstructing the Data**
  - [**``Reconstruct Function``**]()
  - [**``Reconstructed Images``**]()
  
2. **Measure SNR**
  - [**``Location of ROIs``**]()
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

3. **Water Content**
  - The body of the lobster, the center of the coconut and pretty much the whole fish (sans eye) contain a lot of water.
  - The tissues which had the highest difference in proton density seemed to be the fleshy part of the lobster and it's shell. The shell contained almost no signal, but the flesh contained lots.
  - In my opinion, the PD sequence was best for lobster, the T2W signal was best for the fish, and the T1W was best for the coconut. The T1W and T2W signals were quite similar here, but I noticed subtl differences in favour of the selections made above.

4. **Big ass image**
  - The results for this are above (i.e. image with other images, and SNR in previous table).
  
### Part 2: The Noise in MRI

1. **Visualzing Noise Data**
  - [**``Phasor Images for Each Coil``**]()
  - The noise appears to me as approximately Gaussian. The mean isn't quite centered at 0, but slightly higher, which makes sense to me as you can't quite have "negative intensity" (can you?!). The standard deviation value seems relatively uninformative to me here, as the images are on somewhat of an arbitrary intensity scale determined by the gain of the scanner. One important difference I did notice is that noise is non-uniform across the coils. I expect this is due to the location of the coils within the scanner, and heterogeneity of the B0 field causes the coils to react differently in the central, ideal, region of the scanner and the periphery, where that ideality breaks down.

2. **Decomposing Noise Distributions**
  - 
