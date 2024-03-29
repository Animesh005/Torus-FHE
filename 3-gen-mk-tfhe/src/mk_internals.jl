# LWE
using CUDA

"""
A structure representing an encrypted bit in multi-key TFHE.
"""
# struct MKLweSample{T}

#     params :: LweParams
#     a :: T # Array{Torus32, 2} # masks from all parties: (n, parties)
#     b :: Union{Torus32} # the joined phase
#     current_variance :: Float64 # average noise of the sample

#     MKLweSample(params::LweParams, a, b::Torus32, cv::Float64) =
#         new{typeof(a)}(params, a, b, cv)

#     function MKLweSample(params::LweParams, parties::Int)
#         new{Array{Torus32, 2}}(params, Array{Torus32}(undef, params.size, parties), Torus32(0), 0.)
#     end

# end

struct MKLweSample

    params :: LweParams
    a :: Array{Torus32, 2} # masks from all parties: (n, parties)
    b :: Union{Torus32} # the joined phase
    current_variance :: Float64 # average noise of the sample

    MKLweSample(params::LweParams, a::Array{Torus32, 2}, b::Torus32, cv::Float64) =
        new(params, a, b, cv)

    function MKLweSample(params::LweParams, parties::Int)
        new(params, Array{Torus32}(undef, params.size, parties), Torus32(0), 0.)
    end

end


Base.:-(x::MKLweSample, y::MKLweSample) =
    MKLweSample(x.params, x.a .- y.a, x.b - y.b, x.current_variance + y.current_variance)

Base.:-(x::MKLweSample) =
    MKLweSample(x.params, .- x.a, - x.b, x.current_variance)

Base.:+(x::MKLweSample, y::MKLweSample) =
    MKLweSample(x.params, x.a .+ y.a, x.b + y.b, x.current_variance + y.current_variance)


Base.:*(x::Torus32, y::MKLweSample) =
    MKLweSample(y.params, x .* y.a, x * y.b, x*x * y.current_variance)


### GPU Version ###

struct MKLweSampleGPU

    params :: LweParams
    a :: CuArray{Torus32, 2} # masks from all parties: (n, parties)
    b :: Union{Torus32} # the joined phase
    current_variance :: Float64 # average noise of the sample

    MKLweSampleGPU(params::LweParams, a::CuArray{Torus32, 2}, b::Torus32, cv::Float64) =
        new(params, a, b, cv)

    function MKLweSampleGPU(params::LweParams, parties::Int)
        new(params, CuArray{Torus32}(undef, params.size, parties), Torus32(0), 0.)
    end
end

Base.:-(x::MKLweSampleGPU, y::MKLweSampleGPU) =
    MKLweSampleGPU(x.params, x.a .- y.a, x.b - y.b, x.current_variance + y.current_variance)

Base.:-(x::MKLweSampleGPU) =
    MKLweSampleGPU(x.params, .- x.a, - x.b, x.current_variance)

Base.:+(x::MKLweSampleGPU, y::MKLweSampleGPU) =
    MKLweSampleGPU(x.params, x.a .+ y.a, x.b + y.b, x.current_variance + y.current_variance)


Base.:*(x::Torus32, y::MKLweSampleGPU) =
    MKLweSampleGPU(y.params, x .* y.a, x * y.b, x*x * y.current_variance)
    

function mk_lwe_phase(sample::MKLweSample, lwe_keys::Array{LweKey, 1})
    parties = length(lwe_keys)
    phases = [
        lwe_phase(LweSample(sample.params, sample.a[:,i], Torus32(0), 0.), lwe_keys[i])
        for i in 1:parties]
    sample.b + reduce(+, phases)
end


function mk_lwe_noiseless_trivial(mu::Torus32, params::LweParams, parties::Int)
    MKLweSample(params, zeros(typeof(mu), params.size, parties), mu, 0.)
end

# function mk_lwe_noiseless_trivial_gpu(mu::Torus32, params, parties::Int)
#     MKLweSampleGPU(params, CUDA.zeros(typeof(mu), params.size, parties), mu, 0.)
# end


# RLWE


struct MKRLweSample
    params :: RLweParams
    a :: Array{TorusPolynomial, 1} # mask (mask_size = 1, so length = parties)
    b :: TorusPolynomial
    current_variance :: Float64

    function MKRLweSample(
            params::RLweParams, a::Array{<:TorusPolynomial, 1}, b::TorusPolynomial,
            current_variance::Float64)
        new(params, a, b, current_variance)
    end
end


Base.:+(x::MKRLweSample, y::MKRLweSample) =
    MKRLweSample(x.params, x.a .+ y.a, x.b + y.b, x.current_variance + y.current_variance)


Base.:-(x::MKRLweSample, y::MKRLweSample) =
    MKRLweSample(x.params, x.a .- y.a, x.b - y.b, x.current_variance + y.current_variance)


# result = (0, ..., 0, mu)
function mk_rlwe_noiseless_trivial(mu::TorusPolynomial, rlwe_params::RLweParams, parties::Int)
    p_degree = rlwe_params.polynomial_degree

    MKRLweSample(
        rlwe_params,
        [torus_polynomial(zeros(typeof(mu.coeffs[1]), p_degree)) for i in 1:parties],
        mu,
        0.)
end


function mk_mul_by_monomial(sample::MKRLweSample, ai::Int32)
    MKRLweSample(
        sample.params,
        mul_by_monomial.(sample.a, ai),
        mul_by_monomial(sample.b, ai),
        sample.current_variance)
end


function mk_rlwe_extract_sample(x::MKRLweSample)
    # Iterating over parties here, not mask elements! (mask_size = 1)
    # TODO: (issue #2) correct for mask_size > 1
    a = hcat([reverse_polynomial(p).coeffs for p in x.a]...)
    b = x.b.coeffs[1]
    lwe_params = LweParams(x.params.polynomial_degree * x.params.mask_size)
    MKLweSample(lwe_params, a, b, 0.) # TODO: (issue #7) calculate the current variance
end


# Internal keys


struct SharedKey
    tgsw_params :: TGswParams
    rlwe_params :: RLweParams
    a :: Array{TorusPolynomial, 1}

    function SharedKey(rng, tgsw_params::TGswParams, rlwe_params::RLweParams)
        decomp_length = tgsw_params.decomp_length
        p_degree = rlwe_params.polynomial_degree
        rand_uni = rlwe_params.is32 ? rand_uniform_torus32 : rand_uniform_torus64
        a = [torus_polynomial(rand_uni(rng, p_degree)) for i in 1:decomp_length]
        new(tgsw_params, rlwe_params, a)
    end
end

""" Structure that will store the Common Random Polynomials that are being used by the 3rd Generation MKTFHE protocol."""
struct CRP_3gen
    tgsw_params :: TGswParams
    rlwe_params :: RLweParams
    a :: Array{TorusPolynomial, 1}
    transformed_a :: Array{<:Union{TransformedTorusPolynomial, TransformedTorusPolynomial64},1}

    function CRP_3gen(rng, tgsw_params::TGswParams, rlwe_params::RLweParams, a_same::Bool = false)
        decomp_length = tgsw_params.decomp_length
        p_degree = rlwe_params.polynomial_degree
        rand_uni = rlwe_params.is32 ? rand_uniform_torus32 : rand_uniform_torus64
        if a_same==false
            a = [torus_polynomial(rand_uni(rng, p_degree)) for i in 1:decomp_length]
        else
            crp=torus_polynomial(rand_uni(rng, p_degree))
            a = [crp for i in 1:decomp_length]
        end
        transformed_a=forward_transform.(a,rlwe_params.is32)
        new(tgsw_params, rlwe_params, a,transformed_a)
    end
end




struct TransformedSharedKey
    tgsw_params :: TGswParams
    rlwe_params :: RLweParams

    a :: Array{<:Union{TransformedTorusPolynomial, TransformedTorusPolynomial64}, 1}

    TransformedSharedKey(tgsw_params, rlwe_params, a) =
        new(tgsw_params, rlwe_params, a)
end


function forward_transform(sample::SharedKey)
    TransformedSharedKey(
        sample.tgsw_params,
        sample.rlwe_params,
        forward_transform.(sample.a, sample.rlwe_params.is32))
end


struct PublicKey
    tgsw_params :: TGswParams
    rlwe_params :: RLweParams
    b :: Array{TorusPolynomial, 1}

    function PublicKey(
            rng, rlwe_key::RLweKey, alpha::Float64, shared::SharedKey, tgsw_params::TGswParams, wo_FFT :: Int64 = 0)

        rlwe_params = rlwe_key.params
        decomp_length = tgsw_params.decomp_length
        p_degree = rlwe_params.polynomial_degree
        is32 = rlwe_params.is32

        # The right part of the public key is `b_i = e_i + a*s_i`,
        # where `a` is shared between the parties
        # TODO: (issue #2) [1] was omitted in the original!
        #       It works while k=1, but will fail otherwise
        # TODO: (issue #2) this is basically tgsw_encrypt_zero() for mask_size=1

        if wo_FFT==0
            b = rlwe_params.is32 ?
                [(
                    transformed_mul(rlwe_key.key[1], shared.a[i], is32)
                    + torus_polynomial(rand_gaussian_torus32(rng, zero(Int32), alpha, p_degree)))
                for i in 1:decomp_length] :
                [(
                    transformed_mul(rlwe_key.key[1], shared.a[i], is32)
                    + torus_polynomial(rand_gaussian_torus64(rng, zero(Int64), alpha, p_degree)))
                for i in 1:decomp_length]
        else
            b = rlwe_params.is32 ?
                [(
                    rlwe_key.key[1] * shared.a[i]
                    + torus_polynomial(rand_gaussian_torus32(rng, zero(Int32), alpha, p_degree)))
                for i in 1:decomp_length] :
                [(
                    rlwe_key.key[1] * shared.a[i]
                    + torus_polynomial(rand_gaussian_torus64(rng, zero(Int64), alpha, p_degree)))
                for i in 1:decomp_length]
        end
        new(tgsw_params, rlwe_params, b)
    end
    


    """ In 3Gen, We add a new constructor of individual public key with the CRP not the Shared Key, but it's basically the same thing."""
    function PublicKey(
        rng, rlwe_key::RLweKey, alpha::Float64, crp::CRP_3gen, tgsw_params::TGswParams,wo_FFT ::Int64 = 0) 

        rlwe_params = rlwe_key.params
        decomp_length = tgsw_params.decomp_length
        p_degree = rlwe_params.polynomial_degree
        is32 = rlwe_params.is32

        if wo_FFT==0
        
            b = rlwe_params.is32 ?
                [(
                    transformed_mul(rlwe_key.key[1], crp.a[i], is32)
                    + torus_polynomial(rand_gaussian_torus32(rng, zero(Int32), alpha, p_degree)))
                for i in 1:decomp_length] :
                [(
                    transformed_mul(rlwe_key.key[1], crp.a[i], is32)
                    + torus_polynomial(rand_gaussian_torus64(rng, zero(Int64), alpha, p_degree)))
                for i in 1:decomp_length]
        else
            b = rlwe_params.is32 ?
            [(
                rlwe_key.key[1]*crp.a[i]
                + torus_polynomial(rand_gaussian_torus32(rng, zero(Int32), alpha, p_degree)))
            for i in 1:decomp_length] :
            [(
                rlwe_key.key[1]*crp.a[i]
                + torus_polynomial(rand_gaussian_torus64(rng, zero(Int64), alpha, p_degree)))
            for i in 1:decomp_length]
        end

        new(tgsw_params, rlwe_params, b)
    end
        
end


struct TransformedPublicKey
    tgsw_params :: TGswParams
    rlwe_params :: RLweParams

    b :: Array{<:Union{TransformedTorusPolynomial, TransformedTorusPolynomial64}, 1}

    TransformedPublicKey(tgsw_params, rlwe_params, b) =
        new(tgsw_params, rlwe_params, b)
end


function forward_transform(sample::PublicKey)
    TransformedPublicKey(
        sample.tgsw_params,
        sample.rlwe_params,
        forward_transform.(sample.b, sample.rlwe_params.is32))
end


""" Structure that contains the common public key between parties in 3rd Generation MKTFHE. 
The Common Public Key is obtained by adding individual keys that have used the same CRP.
"""
struct CommonPubKey_3gen
    params :: SchemeParameters_3gen
    parties::Int
    b :: Array{TorusPolynomial,1} # decomposed key, because we always used the decomposed common public key
    transformed_b :: Array{<:Union{TransformedTorusPolynomial, TransformedTorusPolynomial64},1}

    function CommonPubKey_3gen(pubkeys::Array{PublicKey,1},params::SchemeParameters_3gen,parties::Int)
        b=Array{TorusPolynomial,1}(undef,params.gsw_decomp_length)
        
        b = pubkeys[1].b

        for i in 2:parties
            b = b .+ pubkeys[i].b
        end
        
        transformed_b=forward_transform.(b,rlwe_parameters(params).is32)

        new(params,parties,b,transformed_b)
    end

end

# TGSW


# Uni-encrypted TGSW sample
struct MKTGswUESample

    tgsw_params :: TGswParams
    rlwe_params :: RLweParams

    c0 :: Array{TorusPolynomial, 1}
    c1 :: Array{TorusPolynomial, 1}
    d0 :: Array{TorusPolynomial, 1}
    d1 :: Array{TorusPolynomial, 1}
    f0 :: Array{TorusPolynomial, 1}
    f1 :: Array{TorusPolynomial, 1}

    current_variance :: Float64 # avg variance of the sample

    function MKTGswUESample(
            tgsw_params::TGswParams, rlwe_params::RLweParams,
            c0::Array{<: TorusPolynomial, 1}, c1::Array{<: TorusPolynomial, 1},
            d0::Array{<: TorusPolynomial, 1}, d1::Array{<: TorusPolynomial, 1},
            f0::Array{<: TorusPolynomial, 1}, f1::Array{<: TorusPolynomial, 1},
            current_variance::Float64)

        decomp_length = tgsw_params.decomp_length
        p_degree = rlwe_params.polynomial_degree

        @assert size(c0) == (decomp_length,)
        @assert size(c1) == (decomp_length,)
        @assert size(d0) == (decomp_length,)
        @assert size(d1) == (decomp_length,)
        @assert size(f0) == (decomp_length,)
        @assert size(f1) == (decomp_length,)

        new(tgsw_params, rlwe_params, c0, c1, d0, d1, f0, f1, current_variance)
    end
end


# Encrypt an integer value
# In the paper: RGSW.UniEnc
# Similar to tgsw_encrypt()/rlwe_encrypt(), except the public key is supplied externally.
function mk_tgsw_encrypt(
    rng, message::Union{Int32, Int64, IntPolynomial}, alpha::Float64,
    rlwe_key::RLweKey, shared_key::SharedKey, public_key::PublicKey)

tgsw_params = shared_key.tgsw_params
rlwe_params = rlwe_key.params
p_degree = rlwe_params.polynomial_degree
decomp_length = tgsw_params.decomp_length
is32 = rlwe_params.is32

if is32
    rand_bool = rand_uniform_bool
    rand_uni_torus = rand_uniform_torus32
    rand_gau_torus = rand_gaussian_torus32
    zero = Int32(0)
else
    rand_bool = rand_uniform_bool64
    rand_uni_torus = rand_uniform_torus64
    rand_gau_torus = rand_gaussian_torus64
    zero = Int64(0)
end


# The shared randomness
r = int_polynomial(rand_bool(rng, p_degree))

# C = (c0,c1) \in T^2dg, with c0 = s_party*c1 + e_c + m*g
c1 = [torus_polynomial(rand_uni_torus(rng, p_degree)) for i in 1:decomp_length]
# TODO: (issue #2) it was just key, not key[1] in the original. Hardcoded mask_size=1? Check!
c0 = [
        (torus_polynomial(rand_gau_torus(rng, zero, alpha, p_degree))
        + transformed_mul(rlwe_key.key[1], c1[i], is32)
        + message * tgsw_params.gadget_values[i])
        for i in 1:decomp_length]

# D = (d0, d1) = r*[Pkey_party | Pkey_parties] + [E0|E1] + [0|m*g] \in T^2dg
d1 = [
        (torus_polynomial(rand_gau_torus(rng, zero, alpha, p_degree))
        + transformed_mul(r, shared_key.a[i], is32)
        + message * tgsw_params.gadget_values[i])
        for i in 1:decomp_length]
d0 = [
        (torus_polynomial(rand_gau_torus(rng, zero, alpha, p_degree))
        + transformed_mul(r, public_key.b[i], is32))
        for i in 1:decomp_length]

# F = (f0,f1) \in T^2dg, with f0 = s_party*f1 + e_f + r*g
f1 = [torus_polynomial(rand_uni_torus(rng, p_degree)) for i in 1:decomp_length]
# TODO: (issue #2) it was just key, not key[1] in the original. Hardcoded mask_size=1? Check!
f0 = [
        (torus_polynomial(rand_gau_torus(rng, zero, alpha, p_degree))
        + transformed_mul(rlwe_key.key[1], f1[i], is32)
        + r * tgsw_params.gadget_values[i])
        for i in 1:decomp_length]

MKTGswUESample(tgsw_params, rlwe_params, c0, c1, d0, d1, f0, f1, alpha^2)
end



struct MKTransformedTGswUESample

    tgsw_params :: TGswParams
    rlwe_params :: RLweParams

    d :: Array{<:Union{TransformedTorusPolynomial, TransformedTorusPolynomial64}, 1}
    f0 :: Array{<:Union{TransformedTorusPolynomial, TransformedTorusPolynomial64}, 1}
    f1 :: Array{<:Union{TransformedTorusPolynomial, TransformedTorusPolynomial64}, 1}

    current_variance :: Float64 # avg variance of the sample

    MKTransformedTGswUESample(tgsw_params, rlwe_params, d, f0, f1, current_variance) =
        new(tgsw_params, rlwe_params, d, f0, f1, current_variance)
end


function forward_transform(sample::MKTGswUESample)
    MKTransformedTGswUESample(
        sample.tgsw_params,
        sample.rlwe_params,
        forward_transform.(sample.d1, sample.rlwe_params.is32),
        forward_transform.(sample.f0, sample.rlwe_params.is32),
        forward_transform.(sample.f1, sample.rlwe_params.is32),
        sample.current_variance)
end


function UniProduct_old(
    acc::MKRLweSample, unienc::MKTransformedTGswUESample, pk::Array{TransformedPublicKey},
    sk::TransformedSharedKey, party::Int, parties::Int)

    tgsw_params = unienc.tgsw_params
    rlwe_params = acc.params
    declen = tgsw_params.decomp_length
    is32 = rlwe_params.is32

    dec_a = hcat(decompose.(acc.a, Ref(tgsw_params)))
    dec_b = decompose(acc.b, tgsw_params)

    tr_dec_a = is32 ? Array{TransformedTorusPolynomial}(undef, parties, declen) :
                      Array{TransformedTorusPolynomial64}(undef, parties, declen)
    tr_dec_b = is32 ? Array{TransformedTorusPolynomial}(undef, declen) :
                      Array{TransformedTorusPolynomial64}(undef, declen)

    for i = 1 : parties
        tr_dec_a[i, :] = forward_transform.(dec_a[i], is32)
    end
    tr_dec_b = forward_transform.(dec_b, is32)

    u = Array{TorusPolynomial}(undef, parties)
    for i = 1 : parties
        u[i] = inverse_transform(sum([tr_dec_a[i, l] * unienc.d[l] for l in 1: declen]), is32)
    end
    u0 = inverse_transform(sum([tr_dec_b[l] * unienc.d[l] for l in 1:declen]), is32)

    v = Array{TorusPolynomial}(undef, parties)
    for i = 1 : parties
        v[i] = inverse_transform(sum([tr_dec_a[i, l] * pk[i].b[l] for l in 1: declen]), is32)
    end
    v0 = -inverse_transform(sum([tr_dec_b[l] * sk.a[l] for l in 1: declen]), is32)

    dec_v = hcat(decompose.(v, Ref(tgsw_params)))
    dec_v0 = decompose(v0, tgsw_params)

    tr_dec_v = is32 ? Array{TransformedTorusPolynomial}(undef, parties, declen) :
                      Array{TransformedTorusPolynomial64}(undef, parties, declen)
    for i = 1 : parties
        tr_dec_v[i, :] = forward_transform.(dec_v[i], is32)
    end
    tr_dec_v0 = forward_transform.(dec_v0, is32)

    w0 = [
            inverse_transform(sum([tr_dec_v[j, l] * unienc.f0[l] for l in 1: declen]), is32)
        for j in 1:parties]
    w00 = inverse_transform(sum([tr_dec_v0[l] * unienc.f0[l] for l in 1: declen]), is32)

    w1 = [
            inverse_transform(sum([tr_dec_v[j, l] * unienc.f1[l] for l in 1: declen]), is32)
        for j in 1:parties]
    w10 = inverse_transform(sum([tr_dec_v0[l] * unienc.f1[l] for l in 1: declen]), is32)

    anew = u
    anew[party] += sum(w1) + w10
    bnew = u0 + sum(w0) + w00

    MKRLweSample(rlwe_params, anew, bnew, .0)
end






# TGSW sample expanded for all the parties.
# The full matrix has size `(parties+1)x(parties+1)` and its elements are
# `decomp_length`-vectors of `TorusPolynomial`.
# Since the C matrix returned by RGSW.Expand in the paper is sparse and has repeating elements,
# we are keeping only nonzero elements of it, namely:
# To save memory, we are only keeping the distinct elements of it, namely:
# c0 == C_{1,1},
# c1 == C_{1,party+1}
# x == C_{i,1}, i = 2...parties+1
# y == C_{i,party+1}, i == 2...parties+1
struct MKTGswExpSample

    tgsw_params :: TGswParams
    rlwe_params :: RLweParams

    x :: Array{TorusPolynomial, 2}
    y :: Array{TorusPolynomial, 2}
    c0 :: Array{TorusPolynomial, 1}
    c1 :: Array{TorusPolynomial, 1}

    current_variance :: Float64 # avg variance of the sample

    function MKTGswExpSample(
            tgsw_params::TGswParams, rlwe_params::RLweParams,
            x::Array{<: TorusPolynomial, 2},
            y::Array{<: TorusPolynomial, 2},
            c0::Array{<: TorusPolynomial, 1},
            c1::Array{<: TorusPolynomial, 1},
            current_variance::Float64)

        decomp_length = tgsw_params.decomp_length
        p_degree = rlwe_params.polynomial_degree
        parties = size(x, 2)

        @assert size(x) == (decomp_length, parties)
        @assert size(y) == (decomp_length, parties)
        @assert size(c0) == (decomp_length,)
        @assert size(c1) == (decomp_length,)

        new(tgsw_params, rlwe_params, x, y, c0, c1, current_variance)
    end
end


struct MKTransformedTGswExpSample

    tgsw_params :: TGswParams
    rlwe_params :: RLweParams

    x :: Array{TransformedTorusPolynomial, 2}
    y :: Array{TransformedTorusPolynomial, 2}
    c0 :: Array{TransformedTorusPolynomial, 1}
    c1 :: Array{TransformedTorusPolynomial, 1}

    current_variance :: Float64 # avg variance of the sample

    MKTransformedTGswExpSample(tgsw_params, rlwe_params, x, y, c0, c1, current_variance) =
        new(tgsw_params, rlwe_params, x, y, c0, c1, current_variance)
end


function forward_transform(sample::MKTGswExpSample)
    MKTransformedTGswExpSample(
        sample.tgsw_params,
        sample.rlwe_params,
        forward_transform.(sample.x),
        forward_transform.(sample.y),
        forward_transform.(sample.c0),
        forward_transform.(sample.c1),
        sample.current_variance)
end


# In the paper: RGSW.Expand
function mk_tgsw_expand(sample::MKTGswUESample, party::Int, public_keys::Array{PublicKey, 1})

    rlwe_params = sample.rlwe_params
    tgsw_params = sample.tgsw_params
    p_degree = sample.rlwe_params.polynomial_degree
    decomp_length = tgsw_params.decomp_length
    parties = length(public_keys)
    is32 = rlwe_params.is32

    # Calculate the distinct elements of the expanded matrix C.
    # See the description of `MKTGswUESample` for details.

    zero_p = is32 ? torus_polynomial(zeros(Torus32, p_degree)) : torus_polynomial(zeros(Torus64, p_degree))

    # g^{-1}(b_i[j] - b_party[j]) = [u_0, ...,u_dg-1]
    decomp = [
        (i == party ?
            nothing : # unused
            decompose(public_keys[i].b[j] - public_keys[party].b[j], tgsw_params)
            )
        for j in 1:decomp_length, i in 1:parties
        ]

    x = [
        (sample.d0[j] + (i == party ?
            zero_p :
            # add xi[j] = <g^{-1}(b_i[j] - b_party[j]), f0>
            sum(transformed_mul(decomp[j,i][l], sample.f0[l], is32) for l in 1:decomp_length)))
        for j in 1:decomp_length, i in 1:parties
        ]

    y = [
        (i == party ?
            sample.d1[j] :
            # add yi[j] = <g^{-1}(b_i[j] - b_party[j]), f1>
            sum(transformed_mul(decomp[j,i][l], sample.f1[l], is32) for l in 1:decomp_length)
            )
        for j in 1:decomp_length, i in 1:parties
        ]

    # TODO: (issue #7) calculate the current variance correctly
    MKTGswExpSample(tgsw_params, rlwe_params, x, y, sample.c0, sample.c1, 0.0)
end


function mk_tgsw_extern_mul(
        sample::MKRLweSample, exp_sample::MKTransformedTGswExpSample, party::Int, parties::Int)

    rlwe_params = sample.params
    tgsw_params = exp_sample.tgsw_params
    p_degree = rlwe_params.polynomial_degree
    decomp_length = tgsw_params.decomp_length
    is32 = rlwe_params.is32

    dec_a = hcat(decompose.(sample.a, Ref(tgsw_params))...)
    dec_b = decompose(sample.b, tgsw_params)

    # It would be faster to do all the summations in transformed space and then transform back.
    # But, because of the transform used (FFT), the precision is limited,
    # and we are close to that limit.
    # Elements in `dec_a` and `dec_b` use `bs_log2_base` bits.
    # Elements in `exp_sample` arrays have full range of Int32 (32 bits).
    # So if we want to do all the calculations in this function in transformed space, that takes
    # `bs_log2_base + 32 + log2(p_degree) + ceil(log2(decomp_length * (parties+1)))` bits
    # == 8 + 32 + 10 + ceil(log2(4 * 3)) = 54 bits - not good.

    tr_dec_a = forward_transform.(dec_a, is32)
    tr_dec_b = forward_transform.(dec_b, is32)

    # c'_i = g^{-1}(a_i)*d1, i<parties, i!=party
    # c'_party = \sum g^{-1}(a_i)*yi + g^{-1}(b)*c1
    a = [
        (i == party ?
            (sum([inverse_transform(tr_dec_a[l,j] * exp_sample.y[l,j], is32)
                for l in 1:decomp_length, j in 1:parties])
            + sum([inverse_transform(tr_dec_b[l] * exp_sample.c1[l], is32)
                for l in 1:decomp_length])) :
            sum([inverse_transform(tr_dec_a[l,i] * exp_sample.y[l,party], is32)
                for l in 1:decomp_length]))
        for i in 1:parties]

    # c'_parties = \sum g^{-1}(a_i)*xi + g^{-1}(b)*c0
    b = (sum([inverse_transform(tr_dec_a[l,i] * exp_sample.x[l,i], is32)
            for l in 1:decomp_length, i in 1:parties])
        + sum([inverse_transform(tr_dec_b[l] * exp_sample.c0[l], is32)
            for l in 1:decomp_length]))

    # TODO: (issue #7) calculate current_variance
    MKRLweSample(rlwe_params, a, b, 0.)
end


# Keyswitch


function mk_keyswitch(ks::Array{KeyswitchKey, 1}, sample::MKLweSample)
    parties = length(ks)
    result_parts = [
        keyswitch(ks[p], LweSample(sample.params, sample.a[:,p], Int32(0), 0.))
        for p in 1:parties]

    out_lwe_params = ks[1].out_lwe_params

    result = mk_lwe_noiseless_trivial(sample.b, out_lwe_params, parties)
    result + MKLweSample(
        out_lwe_params,
        hcat([part.a for part in result_parts]...),
        reduce(+, [part.b for part in result_parts]),
        0.0)
end


""" Modified version of key switch to take into account the switch from a ciphertext of size 2 to a ciphertext of size #parties. """
function mk_keyswitch_3gen(ks::Array{KeyswitchKey, 1}, sample::LweSample)
    parties = length(ks)
    result_parts = [
        keyswitch(ks[p], LweSample(sample.params, sample.a, Int32(0), 0.))
        for p in 1:parties]

    out_lwe_params = ks[1].out_lwe_params

    result = mk_lwe_noiseless_trivial(sample.b, out_lwe_params, parties)
    result + MKLweSample(
        out_lwe_params,
        hcat([part.a for part in result_parts]...),
        reduce(+, [part.b for part in result_parts]),
        0.0)
end


# Bootstrap


# A part of the bootstrap key generated independently by each party
# (since it involves their secret keys).
struct BootstrapKeyPart

    tgsw_params :: TGswParams
    rlwe_params :: RLweParams
    key_uni_enc :: Array{MKTGswUESample, 1}
    public_key :: PublicKey

    function BootstrapKeyPart(
            rng, lwe_key::LweKey, rlwe_key::RLweKey, alpha::Float64,
            shared_key::SharedKey, public_key::PublicKey,wo_FFT::Int64=0)

        tgsw_params = shared_key.tgsw_params
        rlwe_params = shared_key.rlwe_params
        in_out_params = lwe_key.params
        if wo_FFT==0
            uni_enc = [
                mk_tgsw_encrypt(rng, lwe_key.key[j], alpha, rlwe_key, shared_key, public_key)
                for j in 1:in_out_params.size]
        else
            uni_enc = [
                mk_tgsw_encrypt_wo_FFT(rng, lwe_key.key[j], alpha, rlwe_key, shared_key, public_key)
                for j in 1:in_out_params.size]
        end
        new(tgsw_params, rlwe_params, uni_enc, public_key)
    end
end


struct MKBootstrapKey
    tgsw_params :: TGswParams
    rlwe_params :: RLweParams
    pk :: Array{TransformedPublicKey, 1}
    crs :: TransformedSharedKey
    bk :: Array{MKTransformedTGswUESample, 2}

    function MKBootstrapKey(parts::Array{BootstrapKeyPart, 1}, shared_key::SharedKey)

        parties = length(parts)

        transformed_public_keys = forward_transform.([part.public_key for part in parts])

        transformed_samples = [
            forward_transform(parts[i].key_uni_enc[j])
            for j in 1:length(parts[1].key_uni_enc), i in 1:parties]

        transformed_crs = forward_transform(shared_key)

        new(parts[1].tgsw_params, parts[1].rlwe_params,
            transformed_public_keys, transformed_crs, transformed_samples)
    end
end


function mk_mux_rotate(
        accum::MKRLweSample, sample::MKTransformedTGswUESample,
        public_keys :: Array{TransformedPublicKey}, sk::TransformedSharedKey,
        barai::Int32, party::Int, parties::Int)
    # ACC = BKi*[(X^barai-1)*ACC]+ACC
    temp_result = mk_mul_by_monomial(accum, barai) - accum
    accum + UniProduct_old(temp_result, sample, public_keys, sk, party, parties)
end



function mk_blind_rotate(accum::MKRLweSample, bk::MKBootstrapKey, bara::Array{Int32, 2})
    n, parties = size(bara)
    for i in 1:parties
        for j in 1:n
            baraij = bara[j,i]
            if baraij == 0
                continue
            end
            accum = mk_mux_rotate(accum, bk.bk[j,i], bk.pk, bk.crs, baraij, i, parties)
        end
    end
    accum
end


function mk_blind_rotate_and_extract(
        v::TorusPolynomial, bk::MKBootstrapKey, barb::Int32, bara::Array{Int32, 2})
    parties = size(bara, 2)
    testvectbis = mul_by_monomial(v, -barb)
    acc = mk_rlwe_noiseless_trivial(testvectbis, bk.rlwe_params, parties)
    acc = mk_blind_rotate(acc, bk, bara)
    bk.rlwe_params.is32 ? mk_rlwe_extract_sample(acc) : mk_rlwe_extract_sample_64(acc)
end


function mk_bootstrap_wo_keyswitch(bk::MKBootstrapKey, mu::Union{Torus32, Torus64}, x::MKLweSample)

    p_degree = bk.rlwe_params.polynomial_degree

    barb = decode_message(x.b, p_degree * 2)
    bara = decode_message.(x.a, p_degree * 2)

    #the initial testvec = [mu,mu,mu,...,mu]
    testvect = torus_polynomial(repeat([mu], p_degree))

    mk_blind_rotate_and_extract(testvect, bk, barb, bara)
end


function mk_bootstrap(bk::MKBootstrapKey, ks::Array{KeyswitchKey, 1}, mu::Union{Torus32, Torus64}, x::MKLweSample)
    u = mk_bootstrap_wo_keyswitch(bk, mu, x)
    mk_keyswitch(ks, u)
end




""" TGSW without FFT : """




function mk_tgsw_encrypt_wo_FFT(
    rng, message::Union{Int32, Int64, IntPolynomial}, alpha::Float64,
    rlwe_key::RLweKey, shared_key::SharedKey, public_key::PublicKey)

tgsw_params = shared_key.tgsw_params
rlwe_params = rlwe_key.params
p_degree = rlwe_params.polynomial_degree
decomp_length = tgsw_params.decomp_length
is32 = rlwe_params.is32

if is32
    rand_bool = rand_uniform_bool
    rand_uni_torus = rand_uniform_torus32
    rand_gau_torus = rand_gaussian_torus32
    zero = Int32(0)
else
    rand_bool = rand_uniform_bool64
    rand_uni_torus = rand_uniform_torus64
    rand_gau_torus = rand_gaussian_torus64
    zero = Int64(0)
end


# The shared randomness
r = int_polynomial(rand_bool(rng, p_degree))

# C = (c0,c1) \in T^2dg, with c0 = s_party*c1 + e_c + m*g
c1 = [torus_polynomial(rand_uni_torus(rng, p_degree)) for i in 1:decomp_length]
# TODO: (issue #2) it was just key, not key[1] in the original. Hardcoded mask_size=1? Check!
c0 = [
        (torus_polynomial(rand_gau_torus(rng, zero, alpha, p_degree))
        + rlwe_key.key[1] * c1[i]
        + message * tgsw_params.gadget_values[i])
        for i in 1:decomp_length]

# D = (d0, d1) = r*[Pkey_party | Pkey_parties] + [E0|E1] + [0|m*g] \in T^2dg
d1 = [
        (torus_polynomial(rand_gau_torus(rng, zero, alpha, p_degree))
        + r * shared_key.a[i]
        + message * tgsw_params.gadget_values[i])
        for i in 1:decomp_length]
d0 = [
        (torus_polynomial(rand_gau_torus(rng, zero, alpha, p_degree))
        + r * public_key.b[i])
        for i in 1:decomp_length]

# F = (f0,f1) \in T^2dg, with f0 = s_party*f1 + e_f + r*g
f1 = [torus_polynomial(rand_uni_torus(rng, p_degree)) for i in 1:decomp_length]
# TODO: (issue #2) it was just key, not key[1] in the original. Hardcoded mask_size=1? Check!
f0 = [
        (torus_polynomial(rand_gau_torus(rng, zero, alpha, p_degree))
        + rlwe_key.key[1] * f1[i]
        + r * tgsw_params.gadget_values[i])
        for i in 1:decomp_length]

MKTGswUESample(tgsw_params, rlwe_params, c0, c1, d0, d1, f0, f1, alpha^2)
end