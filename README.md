# SphereSurfaceHistogram

<a href="https://www.repostatus.org/#inactive"><img src="https://www.repostatus.org/badges/latest/inactive.svg" alt="Project Status: Inactive – The project has reached a stable, usable state but is no longer being actively developed; support/maintenance will be provided as time allows." /></a> [![](https://github.com/ffreyer/SphereSurfaceHistogram.jl/workflows/Documentation/badge.svg)](https://ffreyer.github.io/SphereSurfaceHistogram.jl/dev)

This package does histogram binning on a unit sphere, i.e. binning of unit vectors. It also provides some plotting utilities.

Binning a single vector takes on the order of `30ns`.

---

Here's an example plot. The left side shows a histogram of 1e6 random unit vectors, sampled correctly. The right side shows `normalize(2rand(3) .- 1.0)`, also one million samples. 

![](https://github.com/ffreyer/SphereSurfaceHistogram.jl/blob/master/docs/src/assets/main_example.png)
