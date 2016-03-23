**Question 1:**
> No, we cannot calculate the thingies. Each pixel has the same phase so the magnetization becomes: <br/> M(t) =  exp(-i\*phi(t))\*sum(p=1,Npixels,Ap)

**Question 2:**
> The matrix is still rank deficient... We add rows but they are not linearly independent from the previous so the rank does not increase. Relaxation wouldn't change anything unless each pixel had different tissue properties

**Question 3:**
> a) Beyond N_dt=N_px & Npe=N_py, nothing, rank is N_p

> b) Until N_dt=N_px & Npe=N_py, nothing, then aliasing because we start losing rank

> c) High cond but full rank at 10 and up, cond number drops significantly at 19. below 10, image is not recovered
