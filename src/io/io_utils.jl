# Low-level helpers for reading structured text files (force constants, POSCAR, etc.)

readline_skip_text!(io, T) = parse(T, first(split(strip(readline(io)))))

@inline function read_svec3!(io, ::Type{T}; conv = T(1.0)) where T
    xs = split(strip(readline(io)))
    return SVector{3,T}(conv .* parse.(T, xs))
end

@inline function read_vec3!(io, out::AbstractVector{T}; conv = T(1.0)) where T
    out .= conv .* parse.(T, split(strip(readline(io))))
    return out
end

@inline function read_mat3_rows!(io, ::Type{T}; conv = T(1.0)) where T
    r1 = read_svec3!(io, T; conv = conv)
    r2 = read_svec3!(io, T; conv = conv)
    r3 = read_svec3!(io, T; conv = conv)
    M = hcat(r1, r2, r3)'   # r1 is first row, etc.
    return SMatrix{3, 3, T, 9}(M)
end

@inline function read_tensor3!(io, ::Type{T}; conv = T(1.0)) where T
    M = @MArray zeros(T, 3, 3, 3)
    for ii in 1:3, jj in 1:3
        @views read_vec3!(io, M[ii, jj, :]; conv = conv)
    end
    return SArray{Tuple{3,3,3}, T, 3, 27}(M)
end

@inline function read_tensor4!(io, ::Type{T}; conv = T(1.0)) where T
    M = @MArray zeros(T, 3, 3, 3, 3)
    for ii in 1:3, jj in 1:3, kk in 1:3
        @views read_vec3!(io, M[ii, jj, kk, :]; conv = conv)
    end
    return SArray{Tuple{3,3,3,3}, T, 4, 81}(M)
end
