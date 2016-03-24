**Question 1:**
> As there is no gradient applied, all pixels have the same phase as one another. Thus, we cannot calculate individual pixel intensities and only get a measurement for the sum of all pixels. The equation becomes:

Insert picture

**Question 2:**
> The encoding matrix takes the shape Nfe x Nx. When we add more frequency encodes, we are adding more rows to the table but they are linearly dependent upon the previous, thus the rank does not increase. If there was relaxation considered in this model, it would only change this if each pixel had different tissue properties.

**Question 3:**
> For answering these questions we'll first define what the rank of these matrices ends up being: rank=min(Npe, Ny)*min(Nfe, Nx)

> a) As we increase Ndt and Npe, we will increase rank until we have reached the number of pixels in the image, at which point the rank is Np. Then, further increase will not affect the rank but will reduce the condition number.

> b) As we decrease Ndt and Npe, we will maintain full rank (Np) until we reach the number of pixels in the X or Y directions for frequency or phase encoding, respectively. Beyond that, we begin to lose rank and then aliasing will occur as multiple pixels will have the same phase accrual.

> c) At Ndt = Npe = 10 (i.e. Np) the rank will be full and the image will be reconstructed. However, the condition number will be very high with only this many samples, meaning that the computation of the inverse is very unstable. Once the number of phase and frequency encodes reaches 19 and upward, the condition number drops considerably, yielding a stable reconstruction.

**Question 4:**
> As can be seen in the equation defined at the start of Question 3, the rank would equal N_pe up until the value of Ny, then the rank would remain Ny. The reconstructed image would have 1 X-pixel

**Question 5:**
> You need twice as many samples to have a low condition number. Rank is equivalent to that in the solely real case, but condition number is super large making the matrix difficult to invert.

**Question 6:**
> Adding noise makes things worse.

| Sigma | Error |
|-------|-------|
|0.1000 | 1.0160|
|0.2000 | 3.2022|
|0.3000 | 9.0983|
|0.4000 |11.5341|
|0.5000 |16.4340|

**Question 7:**
> rank= min(Npe, Ny)*min(Ndt, Nx) == Nx*Ny for recovery. Therefore not interchangeable.

**Question 8:**
> The dimensions are Npe*Ndt x Nx*Ny, which means that it gets pretty fucking big with big images to use this method.

**Question 9:**
> Our encoding matrix grows quickly, but our condition number drops (especially useful if there is noise). Only a problem if my computer sucks and can't handle the rapid growth?

**Question 10:**
> If we worked with continuous values instead of pixels we would instead need to assign continuous functions to calculate the phase accrued at each point (i.e. an infinite number for Ndt and Npe). This would be WAY too computationally complex to actually solve in practice.

**Question 11:**
> Make encoding matrix in latex and attach it

**Question 12:**
> tic toc around line 499 and then tic toc around fourier matrix and compare
