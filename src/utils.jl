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
index order to chemist index order: i,j,k,l -> i,k,l,j

See `chem_to_phys`, `test_hijkl_symmetries`.
"""
function phys_to_chem(two_body)
    @tullio two_body_out[i,k,l,j] := two_body[i,j,k,l]
    return two_body_out
end

"""
    chem_to_phys(two_body)

Convert the rank-four tensor `two_body` representing two-body integrals from chemist
index order to physicist index order: i,j,k,l -> i,l,j,k

See `phys_to_chem`, `test_hijkl_symmetries`.
"""
function chem_to_phys(two_body)
    @tullio two_body_out[i,l,j,k] := two_body[i,j,k,l]
    return two_body_out
end

## TODO: refactor quite a bit.
## 1. allow this to return nothing or throw an error.
## 2. Make each tullio transform into a function that overwrites an input.
##    We don't want to allocate so many arrays.
"""
    test_hijkl_symmetries(t; chemist=true)

Test that rank-4 tensor `t` has the required symmetries for coefficents of the two-electron terms.
If `chemist` is `true`, assume the input is in chemists' order, otherwise in physicists' order.

If `t` is a correct tensor of indices, it must pass the tests. If `t` is a correct tensor of indicies,
but the flag `chemist` is incorrect, it must fail the tests. So this test may be used to discriminiate
between the orderings.

References: HJO (1.4.17), (1.4.38)

See `phys_to_chem`, `chem_to_phys`.
"""
function test_hijkl_symmetries(t_in::Array{T}; chemist=true) where T
    t = ZChop.zchop(t_in)

    if chemist  # chemists' ordering
        @tullio t1[p, q, r, s] := t[q, p, r, s]  # (1.4.38) in HJO book
        @tullio t2[p, q, r, s] := t[p, q, s, r]  # (1.4.38)
        @tullio t3[p, q, r, s] := t[q, p, s, r]  # (1.4.38)
        @tullio t4[p, q, r, s] := t[r, s, p, q]  # (1.4.17)
        @tullio t5[p, q, r, s] := t[r, s, q, p]  # 1 and 4
        @tullio t6[p, q, r, s] := t[s, r, p, q]  # 2 and 4
        @tullio t7[p, q, r, s] := t[s, r, q, p]  # 3 and 4
        test_set = (t1, t2, t3, t4, t5, t6, t7)
    else  # physicists' ordering
        ## Agrees with following from OpenFermion
        ## pqrs = rqps = psrq = srqp = qpsr = rspq = spqr = qrsp
        @tullio t1[p, q, r, s] := t[r, q, p, s]
        @tullio t2[p, q, r, s] := t[p, s, r, q]
        @tullio t3[p, q, r, s] := t[s, r, q, p]
        @tullio t4[p, q, r, s] := t[q, p, s, r]
        @tullio t5[p, q, r, s] := t[r, s, p, q]
        @tullio t6[p, q, r, s] := t[s, p, q, r]
        @tullio t7[p, q, r, s] := t[q, r, s, p]
        test_set = (t1, t2, t3, t4, t5, t6, t7)
    end

    fail_count = 0
    for (i, ts) in enumerate(test_set)
        if ! isapprox(ts, t)
            println("** Failed symmetry test number $i.")
            fail_count += 1
        end
    end
    println(fail_count, "/", length(test_set), " tests failed.")
    return nothing
end
