# SSHBinner

As mentioned before a `SSHBinner` can be created using `SSHBinner(N_bins)`, where `N_bins` is the requested number of bins. Typically the binner will have slightly more or less bins as it optimizes the shape and size of each bin.

```@docs
SSHBinner(::Real)
```

Beyond this a binner can also be created from the bins of another binner. This is useful for saving and loading it. Instead of saving the whole structure, you may save just `binner.bins` and recover the full object using `SSHBinner(bins)`.

```@docs
SSHBinner(::Vector)
```

You can also pick between two algorithms for generating bins. This can be done using the `method` keyword argument in `SSHBinner` to pass either

```@docs
partition_sphere1
partition_sphere2
```

Note that both methods are optimized by another function, so that the area of each bin is equal.
