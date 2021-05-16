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
"""
function phys_to_chem(two_body)
    @tullio two_body_out[i,k,l,j] := two_body[i,j,k,l]
    return two_body_out
end

function chem_to_phys(two_body)
    @tullio two_body_out[i,l,j,k] := two_body[i,j,k,l]
    return two_body_out
end

"""
    test_hijkl_symmetries(t; chemist=false)

Test that rank-4 tensor `t` has the required symmetries for coefficents of the two body terms.
If `chemist` is `true`, assume the input is in chemists' order, otherwise in physicists' order.

If `t` is a correct tensor of indices, it must pass the tests. If `t` is a correct tensor of indicies,
but the flag `chemist` is incorrect, it must fail the tests. So this test may be used to discriminiate
between the orderings.

References: HJO (1.4.17), (1.4.38)
"""
function test_hijkl_symmetries(t_in::Array{T}; chemist=false) where T
    t = ZChop.zchop(t_in)
    ## These are some of the required symmetries. This may be all.
    if chemist  # chemist ordering
        @tullio t1[p, q, r, s] := t[q, p, r, s]  # (1.4.38) in HJO book
        @tullio t2[p, q, r, s] := t[p, q, s, r]  # (1.4.38)
        @tullio t3[p, q, r, s] := t[q, p, s, r]  # (1.4.38)
        @tullio t4[p, q, r, s] := t[r, s, p, q]  # (1.4.17)
        @tullio t5[p, q, r, s] := t[r, s, q, p]  # 1 and 4
        @tullio t6[p, q, r, s] := t[s, r, p, q]  # 2 and 4
        @tullio t7[p, q, r, s] := t[s, r, q, p]  # 3 and 4
        test_set = (t1, t2, t3, t4, t5, t6, t7)
    else  # physicist ordering
        @tullio t1[i, j, k, l] := t[j, i, l, k]  # swap a^ s and a s
        @tullio t2[i, j, k, l] := t[j, i, k, l]  # Swap only a^ s; should have -1 phase, but does not!
        @tullio t3[i, j, k, l] := t[i, j, l, k]  # Swap only a s;  should have -1 phase, but does not!
        @tullio t4[i, j, k, l] := conj(t[k, l, i, j])  # h.c. of i,j,k,l
        @tullio t5[i, j, k, l] := conj(t[l, k, j, i])  # h.c. and 1
        @tullio t6[i, j, k, l] := conj(t[k, l, j, i])  # h.c. and 2
        @tullio t7[i, j, k, l] := conj(t[l, k, i, j])  # h.c. and 3
        test_set = (t1, t2, t3, t4, t5, t6, t7)
    end
         ## I'm not sure of the explanations for these symmetries.
        # @tullio t1[i, j, k, l] := t[k, l, i, j]  # swap a^ s and a s. HJO (1.4.17)
        # @tullio t2[i, j, k, l] := t[j, l, i, k]  # Swap only a^ s; should have -1 phase, but does not!
        # @tullio t3[i, j, k, l] := t[i, k, j, l]  # Swap only a s;  should have -1 phase, but does not!
        # @tullio t4[i, j, k, l] := conj(t[l, k, j, i])  # h.c.
        # @tullio t5[i, j, k, l] := conj(t[j, i, l, k])  # h.c. and 1
        # @tullio t6[i, j, k, l] := conj(t[k, i, l, j])  # h.c. and 2
        # @tullio t7[i, j, k, l] := conj(t[l, j, k, i])  # h.c. and 3
        # test_set = (t1, t2, t3, t4, t5, t6, t7)
    ## Here are 8-point symmetry from OF interaction operator.
    ## I don't know what they mean. They correspond to neither our chem nor phys symmetries.
    ## pqrs = rqps = psrq = srqp = qpsr = rspq = spqr = qrsp.
    ## ijkl = kjil = ilkj = lkji = jilk = klij = lijk = jkli.

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
