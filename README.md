# SphereSurfaceHistogram

[![github-ci](https://github.com/ffreyer/SphereSurfaceHistogram.jl/workflows/Tests/badge.svg)](https://github.com/ffreyer/SphereSurfaceHistogram.jl/actions?query=workflow%3A%22Tests%22)
[![codecov](https://codecov.io/gh/ffreyer/SphereSurfaceHistogram.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/ffreyer/SphereSurfaceHistogram.jl)
[![](https://img.shields.io/badge/docs-dev-blue.svg)](https://ffreyer.github.io/SphereSurfaceHistogram.jl/dev)


This package does histogram binning on a unit sphere, i.e. binning of unit vectors. It also provides some plotting utilities.

Binning a single vector takes on the order of `25ns` for `10_000` bins and `35ns` for `1_000_000` bins on my machine.

---

Here's an example plot. The left side shows a histogram of 1e6 random unit vectors, sampled correctly. The right side shows `normalize(2rand(3) .- 1.0)`, also one million samples. 

![example](https://user-images.githubusercontent.com/10947937/110827691-9f0bfb80-8296-11eb-9caa-b18742f41c8f.png)
