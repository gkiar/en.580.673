**Question 1:**
> No, we cannot calculate the thingies. Each pixel has the same phase so the magnetization becomes: <br/> M(t) =  exp(-i\*phi(t))\*sum(p=1,Npixels,Ap)

**Question 2:**
> The matrix is still rank deficient... We add rows but they are not linearly independent from the previous so the rank does not increase. Relaxation wouldn't change anything unless each pixel had different tissue properties

**Question 3:**
> a) Beyond N_dt=N_px & Npe=N_py, nothing, rank is N_p

> b) Until N_dt=N_px & Npe=N_py, nothing, then aliasing because we start losing rank

> c) High cond but full rank at 10 and up, cond number drops significantly at 19. below 10, image is not recovered

**Question 4:**
> Rank would equal N_pe up until the value of Ny, then the rank would remain Ny. We would effectively recover an image with Nx = 1 (all x pixels blurred together for each row). rank= min(Npe, Ny)*min(Nfe, Nx)

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
