################################################################################
### SSHBinner utilities (things from Base and more)
################################################################################


function Base.show(io::IO, B::SSHBinner)
    print(io, "SSHBinner with $(length(B.bins)) bins")
end

Base.summary(io::IO, B::SSHBinner) = print(io, "$(typeof(B)){$(length(B))}")

# "Return the number of elements in the collection."
# would be sum(B.bins)
# Maybe length(B.bins) makes more sense?
Base.length(B::SSHBinner) = length(B.bins)

# TODO: What should this be exactly?
# - length per theta? -> phi_divisions
# - summarized length per theta? (1x 1, 4x 4, 20x 32, ...)
Base.size(B::SSHBinner) = B.phi_divisions

# Clears bins
function Base.empty!(B::SSHBinner)
    B.bins .= 0
    B
end

Base.maximum(B::SSHBinner) = maximum(B.bins)
Base.minimum(B::SSHBinner) = minimum(B.bins)
end
