"""
    non_zero_elements(t)

Return an iterator over non-zero elements of `t` as tuples whose first element
is a tuple of indices and whose second element is the value of `t` at those indices.
"""
function non_zero_elements(t)
    ((ind, t[ind...]) for ind in Iterators.product(axes(t)...) if !iszero(t[ind...]))
end

"""
    non_zero_elements_python(t)

Like `non_zero_elements`, but use Python's order for product of iterators. Also use zero-based indexing in
reported indices.
"""
function non_zero_elements_python(t)
    ((reverse(ind .- 1), t[reverse(ind)...]) for ind in Iterators.product(axes(t)...) if !iszero(t[reverse(ind)...]))
end

"""
    phys_to_chem(two_body)

Convert the rank-four tensor `two_body` representing two-body integrals from physicist
index order to chemist index order: i,j,k,l -> i,l,j,k

See `chem_to_phys`, `check_two_body_symmetries`.
"""
function phys_to_chem(two_body)
    @tullio two_body_out[i,l,j,k] := two_body[i,j,k,l]
    return two_body_out
end

"""
    chem_to_phys(two_body)

Convert the rank-four tensor `two_body` representing two-body integrals from chemist
index order to physicist index order: i,j,k,l -> i,k,l,j

See `phys_to_chem`, `check_two_body_symmetries`.

# Note
Either `chem_to_phys` or `phys_to_chem` generates the cyclic group of order three.
Denote `chem_to_phys` by `g` and `phys_to_chem` by `h`. The elements `g`, `h`, `I` form
a group with `gh = hg = I`, `g^2=h`, and `h^2=g`.
"""
function chem_to_phys(two_body)
    @tullio two_body_out[i,k,l,j] := two_body[i,j,k,l]
    ## Following is tranform used in qiskit_nature. It is the same as
    ## above, because of symmetry.
    #  @tullio two_body_out[l,j,i,k] := two_body[i,j,k,l]
    return two_body_out
end

"""
    to_chem(two_body_tensor)

Permute `two_body_tensor` to chemists' index order, if it is not already
in this order and return the new tensor.

See `phys_to_chem`, `chem_to_phys`.
"""
function to_chem(two_body_tensor)
    if check_two_body_symmetries(two_body_tensor)
        return two_body_tensor
    end
    transformed_tensor = phys_to_chem(two_body_tensor)
    if check_two_body_symmetries(transformed_tensor)
        return transformed_tensor
    end
    transformed_tensor = phys_to_chem(transformed_tensor)
    if check_two_body_symmetries(transformed_tensor)
        return transformed_tensor
    end
    throw(ArgumentError("Unable to permute `two_body_tensor` to chemists' index order"))
end

"""
    to_phys(two_body_tensor)

Permute `two_body_tensor` to physicists' index order, if it is not already
in this order and return the new tensor.

# Throws
- `ArgumentError`: If the algorithm fails to permute `two_body_tensor` to physicists' order.
"""
function to_phys(two_body_tensor)
    index_order = find_index_order(two_body_tensor)
    if index_order == :physicist
        return two_body_tensor
    end
    if index_order == :chemist
        return chem_to_phys(two_body_tensor)
    end
    if index_order == :intermediate
        return phys_to_chem(two_body_tensor)
    end
    throw(ArgumentError("Unable to permute `two_body_tensor` to physicists' index order"))
end

## This routine (and many others here) are quite inefficient. But, they are probably
## more clear as written, and furthermore they are unlikely to be perf bottlenecks.
function to_intermediate(two_body_tensor)
    index_order = find_index_order(two_body_tensor)
    if index_order == :intermediate
        return two_body_tensor
    end
    if index_order == :chemist  # doing the "wrong" permutation takes chem to intermediate.
        return phys_to_chem(two_body_tensor)
    end
    if index_order == :physicist  # doing the "wrong" permutation takes phys to intermediate.
        return chem_to_phys(two_body_tensor)
    end
    throw(ArgumentError("Unable to permute `two_body_tensor` to intermediate' index order"))
end

"""
    to_index_order(two_body_tensor, index_order)

Permute the rank-four tensor `two_body_tensor` so that the indices
are in the order specified by `index_order`, which must be one of
`:chemist`, `:physicist`, or `:intermediate`.

`two_body_tensor` must be in one of the three index orders in order that
the permutation algorithm be successful. This function is relatively inefficient.
It first performs sequences of symmetry tests to detemine the index order of the
input tensor. It then performs the correct sequence of permutations to bring the
tensor to the desired input order. If performance is important, and you know the
index order of the input tensor, it is better to apply just the single required
permutation.

# Throws:
- `ArgumentError` if the algorithm is unable to obtain the desired index order.
"""
function to_index_order(two_body_tensor, index_order)
    if index_order == :chemist
        return to_chem(two_body_tensor)
    elseif index_order == :physicist
        return to_phys(two_body_tensor)
    elseif index_order == :intermediate
        return to_intermediate(two_body_tensor)
    else
        throw(ArgumentError("`index_order` must be one of `:chemist`, `:physicist`, or `:intermediate`"))
    end
end

"""
    check_two_body_symmetries(two_body_tensor; chemist=true)

Return `true` if the rank-4 tensor `two_body_tensor` has the required symmetries for coefficents of the
two-electron terms.  If `chemist` is `true`, assume the input is in chemists' order,
otherwise in physicists' order.

If `two_body_tensor` is a correct tensor of indices, with the correct index order, it must pass the
tests. If `two_body_tensor` is a correct tensor of indicies, but the flag `chemist` is incorrect, it will
fail the tests, unless the tensor has accidental symmetries.
This test may be used with care to discriminiate between the orderings.

References: HJO Molecular Electronic-Structure Theory (1.4.17), (1.4.38)

See `phys_to_chem`, `chem_to_phys`.
"""
function check_two_body_symmetries(two_body_tensor_in; chemist=true)
    two_body_tensor = ZChop.zchop(two_body_tensor_in)
    if ! chemist
        ## Safest is to use property that if a tensor satisfies phys ordered symmetries,
        ## then it satisfies chem ordered symmetries after transformation
        two_body_tensor = phys_to_chem(two_body_tensor)
    end
    t1 = similar(two_body_tensor) # `t1` will hold the transformed tensor
    for test! in _chem_tests
        test!(t1, two_body_tensor)  # write transformed tensor to `t1`
        if ! (t1 ≈ two_body_tensor)
            return false
        end
    end
    return true
end

"""
    find_index_order(two_body_tensor)

Return the index convention of rank-four `two_body_tensor`.

The index convention is determined by checking symmetries of the tensor.
If the indexing convention can be determined, then one of `:chemist`,
`:physicist`, or `:intermediate` is returned. The `:intermediate` indexing
may be obtained by applying `chem_to_phys` to the physicists' convention or
`phys_to_chem` to the chemists' convention. If the tests for each of these
conventions fail, then `:unknown` is returned.

See also: `chem_to_phys`, `phys_to_chem`, `check_two_body_symmetries`.

# Note
The first of `:chemist`, `:physicist`, and `:intermediate`, in that order, to pass the tests
is returned. If `two_body_tensor` has accidental symmetries, it may in fact satisfy more
than one set of symmetry tests. For example, if all elements have the same value, then the
symmetries for all three index orders are satisfied.
"""
function find_index_order(two_body_tensor)
    if check_two_body_symmetries(two_body_tensor)
        return :chemist
    else
        transformed_tensor = phys_to_chem(two_body_tensor)
        if check_two_body_symmetries(transformed_tensor)
            return :physicist
        end
    end
    transformed_tensor = phys_to_chem(transformed_tensor)
    if check_two_body_symmetries(transformed_tensor)
        return :intermediate
    end
    return :unknown
end

## Tests for two-electron coefficient symmetries in chemists' order.
## These anonymous functions each take two arguments. `t` is the tensor
## to be tested and `t1` is a tensor of the same shape that will be overwritten
## with the transformed tensor. The calling code should then compare `t` and
## `t1` for approximate equality.
const _chem_tests = [
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[q, p, r, s]),  # (1.4.38) in HJO book
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[p, q, s, r]),  # (1.4.38)
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[q, p, s, r]),  # (1.4.38)
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[r, s, p, q]),  # (1.4.17)
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[r, s, q, p]),  # 1 and 4
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[s, r, p, q]),  # 2 and 4
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[s, r, q, p]),  # 3 and 4
    ]


function data_header_dict()
    d = Dict()
    d["version"] = PkgVersion.Version(ElectronicStructure)
    d["date"] = Dates.now()
    return d
end
