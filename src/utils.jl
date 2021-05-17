"""
    non_zero_elements(t)

Return an iterator over non-zero elements of `t` as tuples whose first element
is a tuple of indices and whose second element is the value of `t` at those indices.
"""
function non_zero_elements(t)
    ((ind, t[ind...]) for ind in Iterators.product(axes(t)...) if !iszero(t[ind...]))
#   ((ind, t[reverse(ind)...]) for ind in Iterators.product(axes(t)...) if !iszero(t[reverse(ind)...]))
end

"""
    non_zero_elements_python(t)

Like `non_zero_elements`, but use Pythons order for product of iterators. Also use zero-based indexing in
reported indices.
"""
function non_zero_elements_python(t)
    ((reverse(ind .- 1), t[reverse(ind)...]) for ind in Iterators.product(axes(t)...) if !iszero(t[reverse(ind)...]))
end

"""
    phys_to_chem(two_body)

Convert the rank-four tensor `two_body` representing two-body integrals from physicist
index order to chemist index order: i,j,k,l -> i,l,j,k

See `chem_to_phys`, `test_two_body_symmetries`.
"""
function phys_to_chem(two_body)
    @tullio two_body_out[i,l,j,k] := two_body[i,j,k,l]
    return two_body_out
end

"""
    chem_to_phys(two_body)

Convert the rank-four tensor `two_body` representing two-body integrals from chemist
index order to physicist index order: i,j,k,l -> i,k,l,j

Denote `chem_to_phys` by `g` and `phys_to_chem` by `h`. The elements `g`, `h`, `I` form
a group with `gh = hg = I`, `g^2=h`, and `h^2=g`.

See `phys_to_chem`, `test_two_body_symmetries`.
"""
function chem_to_phys(two_body)
    @tullio two_body_out[i,k,l,j] := two_body[i,j,k,l]

    ## Wrong, these are backwards, I think.
#    @tullio two_body_out[i,l,j,k] := two_body[i,j,k,l]
#    @tullio two_body_out[i,j,k,l] := two_body[i,k,l,j]  # Same as above
    return two_body_out
end

## Checking broken physicist symmetry functions
function try_symmetries(two_body; chemist=true)
    t1 = similar(two_body)
    if chemist
        tests = _chem_tests
    else
        tests = _phys_tests
    end
    for tst in tests
        tst(t1, two_body)
        println(t1 ≈ two_body)
    end
end

"""
    check_two_body_symmetries(t; chemist=true)

Return `true` if the rank-4 tensor `t` has the required symmetries for coefficents of the two-electron terms.
If `chemist` is `true`, assume the input is in chemists' order, otherwise in physicists' order.

If `t` is a correct tensor of indices, it must pass the tests. If `t` is a correct tensor of indicies,
but the flag `chemist` is incorrect, it must fail the tests. So this test may be used to discriminiate
between the orderings.

References: HJO (1.4.17), (1.4.38)

See `phys_to_chem`, `chem_to_phys`.
"""
function check_two_body_symmetries(two_body_tensor_in; chemist=true)
    two_body_tensor = ZChop.zchop(two_body_tensor_in)
    if ! chemist
        ## Safest is to use property that if a tensor satisfies phys ordered symmetries,
        ## then it satisfies chem ordered symmetries under transformation
        two_body_tensor = phys_to_chem(two_body_tensor)
    end
    tests = _chem_tests
    t1 = similar(two_body_tensor)
    for test! in tests
        test!(t1, two_body_tensor)
        if ! (t1 ≈ two_body_tensor)
            return false
        end
    end
    return true
end

## Tests for two-electron coefficient symmetries in chemists' order
const _chem_tests = [
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[q, p, r, s]),  # (1.4.38) in HJO book
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[p, q, s, r]),  # (1.4.38)
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[q, p, s, r]),  # (1.4.38)
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[r, s, p, q]),  # (1.4.17)
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[r, s, q, p]),  # 1 and 4
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[s, r, p, q]),  # 2 and 4
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[s, r, q, p]),  # 3 and 4
    ]

## !!!!!!!!! WRONG probably, apparently.
## Tests for two-electron coefficient symmetries in physicists' order
## pqrs = rqps = psrq = srqp = qpsr = rspq = spqr = qrsp
const _phys_tests = [
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[r, q, p, s]),
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[p, s, r, q]),
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[s, r, q, p]),
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[q, p, s, r]),
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[r, s, p, q]),
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[s, p, q, r]),
        ((t1, t) ->  @tullio t1[p, q, r, s] = t[q, r, s, p])
    ]
