# Helper functions
readline_skip_text!(io, T) = parse(T, first(split(strip(readline(io)))))

function read_vec3!(io, T)
    xs = split(strip(readline(io)))
    @assert length(xs) == 3 "Expected 3 components for a 3-vector."
    SVector{3,T}(parse.(T, xs))
end


# lattice vecs are columns of cell matrix
to_frac_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L = mod.(cell \ position, L(1.0))

# lattice vecs are columns of cell matrix
to_cart_coords(cell::AbstractMatrix{L}, position::AbstractVector{L}) where L = cell * position
