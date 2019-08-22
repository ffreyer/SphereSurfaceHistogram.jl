"""
    to_hue(B::SSHBinner[, s=1.0, v=1.0, a=1.0])

Maps each bin filling to a hue relative to the maximum bin filling.
"""
function to_hue(B::SSHBinner, s=1.0, v=1.0, a=1.0)
    # hue is cyclic, so the maximum should be < 1.0 (360)
    bins = 0.8 * B.bins ./ maximum(B.bins)
    map(bins) do k
        HSVA(360k, s, v, a)
    end
end


"""
    to_hue2(B::SSHBinner[, s=1.0, v=1.0, a=1.0])

Maps each bin filling to a hue relative to the maximum and minimum bin filling.
"""
function to_hue2(B::SSHBinner, s=1.0, v=1.0, a=1.0)
    # hue is cyclic, so the maximum should be < 1.0 (360)
    bins = 0.8 * (B.bins .- minimum(B.bins)) ./ (maximum(B.bins) - minimum(B.bins))
    map(bins) do k
        HSVA(360k, s, v, a)
    end
end

function to_alpha_hue(B::SSHBinner, s=1.0, v=1.0)
    # hue is cyclic, so the maximum should be < 1.0 (360)
    bins = 0.8 * (B.bins .- minimum(B.bins)) ./ (maximum(B.bins) - minimum(B.bins))
    map(bins) do k
        HSVA(360k, s, v, min(1.0, 2k))
    end
end


"""
    to_alpha(B::SSHBinner[, r=0.2, g=0.4, b=0.8])

Maps each bin filling to the alpha level of a color (given through r, g, b).
"""
function to_alpha(B::SSHBinner, r=0.2, g=0.4, b = 0.8)
    bins = B.bins ./ maximum(B.bins)
    map(bins) do k
        RGBA(r, g, b, k)
    end
end
