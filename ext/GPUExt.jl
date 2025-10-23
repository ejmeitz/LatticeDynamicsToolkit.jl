module GPUExt

using TDEPToolkit
using TensorOperations
using KernelAbstractions
using StaticArrays

#!TODO UPDATE TO FLOAT32

const GPU_FT = Float32



struct IFC2_GPU{T}
    a1::Int32
    a2::Int32
    ifcs::SMatrix{3,3,T,9}
end

struct IFC3_GPU{T}
    a1::Int32
    a2::Int32
    a3::Int32
    ifcs::SArray{Tuple{3,3,3},T,3,27}
end

struct IFC4_GPU{T}
    a1::Int32
    a2::Int32
    a3::Int32
    a4::Int32
    ifcs::SArray{Tuple{3,3,3,3},T,4,81}
end

function pack_ifc2_flat(fc2, na::Integer)

    n2 = 0
    @inbounds for a1 in 1:na
        n2 += length(get_interactions(fc2, a1))
    end
    buf = Vector{IFC2_GPU{GPU_FT}}(undef, n2)

    i = 1
    @inbounds for a1 in 1:na
        for pair in get_interactions(fc2, a1)
            a2 = Int32(pair.idxs[2])
            K  = SMatrix{3,3,GPU_FT,9}(pair.ifcs)
            buf[i] = IFC2_GPU{GPU_FT}(Int32(a1), a2, K)
            i += 1
        end
    end
    return buf
end


function pack_ifc3_flat(fc3, na::Integer)
    n3 = 0
    @inbounds for a1 in 1:na
        n3 += length(get_interactions(fc3, a1))
    end
    buf = Vector{IFC3_GPU{GPU_FT}}(undef, n3)

    i = 1
    @inbounds for a1 in 1:na
        for trip in get_interactions(fc3, a1)
            a2 = Int32(trip.idxs[2]); a3 = Int32(trip.idxs[3])
            M  = SArray{Tuple{3,3,3},GPU_FT,3,27}(trip.ifcs)
            buf[i] = IFC3_GPU{GPU_FT}(Int32(a1), a2, a3, M)
            i += 1
        end
    end
    return buf
end


function pack_ifc4_flat(fc4, na::Integer)
    n4 = 0
    @inbounds for a1 in 1:na
        n4 += length(get_interactions(fc4, a1))
    end
    buf = Vector{IFC4_GPU{GPU_FT}}(undef, n4)

    i = 1
    @inbounds for a1 in 1:na
        for quat in get_interactions(fc4, a1)
            a2 = Int32(quat.idxs[2]); a3 = Int32(quat.idxs[3]); a4 = Int32(quat.idxs[4])
            Q  = SArray{Tuple{3,3,3,3},GPU_FT,4,81}(quat.ifcs)
            buf[i] = IFC4_GPU{GPU_FT}(Int32(a1), a2, a3, a4, Q)
            i += 1
        end
    end
    return buf
end


sort_by_a1!(v::Vector{<:Union{IFC2_GPU,IFC3_GPU,IFC4_GPU}}) = sort!(v; by = x->x.a1)


@kernel function e2_interactions_kernel_batched(U, pairs2, e2_atom)
    p = @index(Global, 1)  # interaction id
    s = @index(Global, 2)  # sample id
    if p > length(pairs2) || s > size(U,2); return; end
    pr = pairs2[p]
    v  = pr.K * U[pr.a2, s]
    e  = dot(U[pr.a1, s], v)
    Atomix.@atomic e2_atom[pr.a1, s] += e
end

@kernel function e3_interactions_kernel_batched(U, trips3, e3_atom)
    T = eltype(first(U))
    p = @index(Global, 1); s = @index(Global, 2)
    if p > length(trips3) || s > size(U,2); return; end
    tr = trips3[p]
    u2 = U[tr.a2, s]; u3 = U[tr.a3, s]
    c1=c2=c3=zero(T)
    @inbounds for j in 1:3
        u2j = u2[j]
        for k in 1:3
            t = u2j * u3[k]
            c1 += tr.M[1,j,k] * t
            c2 += tr.M[2,j,k] * t
            c3 += tr.M[3,j,k] * t
        end
    end
    e = U[tr.a1, s][1]*c1 + U[tr.a1, s][2]*c2 + U[tr.a1, s][3]*c3
    Atomix.@atomic e3_atom[tr.a1, s] += e
end

@kernel function e4_interactions_kernel_batched(U, quads4, e4_atom)
    T = eltype(first(U))
    p = @index(Global, 1); s = @index(Global, 2)
    if p > length(quads4) || s > size(U,2); return; end
    q = quads4[p]
    u2 = U[q.a2, s]; u3 = U[q.a3, s]; u4 = U[q.a4, s]
    c1=c2=c3=zero(T)
    @inbounds for j in 1:3
        u2j = u2[j]
        for k in 1:3
            t = u2j * u3[k]
            for l in 1:3
                prod = t * u4[l]
                c1 += q.Q[1,j,k,l] * prod
                c2 += q.Q[2,j,k,l] * prod
                c3 += q.Q[3,j,k,l] * prod
            end
        end
    end
    u1 = U[q.a1, s]
    e  = u1[1]*c1 + u1[2]*c2 + u1[3]*c3
    Atomix.@atomic e4_atom[q.a1, s] += e
end

#! UPDATE FOR GENERIC BACKEND
function energies_interaction_mode_batched(
    U_d::CuArray{SVector{3,Float64},2},
    ifc2_flat, ifc3_flat, ifc4_flat; do3::Bool, do4::Bool
)
    na, Ns = size(U_d)
    e2_atom = CUDA.zeros(Float64, na, Ns)
    e3_atom = do3 ? CUDA.zeros(Float64, na, Ns) : CuArray{Float64}(undef, 0, 0)
    e4_atom = do4 ? CUDA.zeros(Float64, na, Ns) : CuArray{Float64}(undef, 0, 0)

    if length(ifc2_flat) > 0
        e2_interactions_kernel_batched(CUDADevice())(
            U_d, ifc2_flat, e2_atom; ndrange=(length(ifc2_flat), Ns)
        )
    end
    if do3 && length(ifc3_flat) > 0
        e3_interactions_kernel_batched(CUDADevice())(
            U_d, ifc3_flat, e3_atom; ndrange=(length(ifc3_flat), Ns)
        )
    end
    if do4 && length(ifc4_flat) > 0
        e4_interactions_kernel_batched(CUDADevice())(
            U_d, ifc4_flat, e4_atom; ndrange=(length(ifc4_flat), Ns)
        )
    end

    e2 = 0.5    .* vec(sum(e2_atom; dims=1))
    e3 = do3 ? (1/6)  .* vec(sum(e3_atom; dims=1)) : zeros(Float64, Ns)
    e4 = do4 ? (1/24) .* vec(sum(e4_atom; dims=1)) : zeros(Float64, Ns)
    return (e2, e3, e4)
end


end