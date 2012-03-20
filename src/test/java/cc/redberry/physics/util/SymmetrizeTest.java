/*
 * Redberry: symbolic tensor computations.
 *
 * Copyright (c) 2010-2012:
 *   Stanislav Poslavsky   <stvlpos@mail.ru>
 *   Bolotin Dmitriy       <bolotin.dmitriy@gmail.com>
 *
 * This file is part of Redberry.
 *
 * Redberry is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * Redberry is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Redberry. If not, see <http://www.gnu.org/licenses/>.
 */
package cc.redberry.physics.util;

import cc.redberry.core.context.CC;
import cc.redberry.core.indices.IndexType;
import cc.redberry.core.tensor.Sum;
import cc.redberry.core.tensor.Tensor;
import static cc.redberry.physics.TAssert.*;
import static org.junit.Assert.assertEquals;
import cc.redberry.transformation.Transformations;
import org.junit.Test;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class SymmetrizeTest {
    public SymmetrizeTest() {
    }

    @Test
    public void test1() {
        System.out.println(Symmetrize.INSTANCE.transform(CC.parse("A_m*B_n")));
    }

    @Test
    public void test2() {
        System.out.println(Symmetrize.INSTANCE.transform(CC.parse("d_a^b*d_c^d")));
    }

    @Test
    public void test3() {
        Tensor t = Transformations.expandBrackets(Symmetrize.INSTANCE.transform(CC.parse("d_p^a*d_q^b*d_r^c")));
        assertTrue(((Sum) t).size() == 6);
    }

    @Test
    public void test4() {
        Tensor t = Transformations.expandBracketsExceptSymbols(Symmetrize.INSTANCE.transform(CC.parse("6*(-(1/2)+l*b*b)*g_pq*g^ab*d_r^c")));
        assertTrue(((Sum) t).size() == 9);
    }

    @Test
    public void test5() {
        Tensor t = Transformations.expandBracketsExceptSymbols(Symmetrize.INSTANCE.transform(CC.parse("3*(-1+l)*n_p*n^a*d_q^b*d_r^c")));
        assertTrue(((Sum) t).size() == 18);
    }

    @Test
    public void test6() {
        Tensor t = Transformations.expandBracketsExceptSymbols(Symmetrize.INSTANCE.transform(CC.parse("6*((1/2)+l*b)*(n_p*n_q*g^ab*d_r^c+n^a*n^b*g_pq*d_r^c)")));
        assertTrue(((Sum) t).size() == 18);
    }

    @Test
    public void test7() {
        Tensor t = Transformations.expandBracketsExceptSymbols(Symmetrize.INSTANCE.transform(CC.parse("6*(-(1/4)+l*b*b)*n_p*g_qr*n^a*g^bc")));
        assertTrue(((Sum) t).size() == 9);
    }

    @Test
    public void test8() {
        Tensor toInv = CC.parse("d_p^a*d_q^b*d_r^c+"
                + "6*(-(1/2)+l*b*b)*g_pq*g^ab*d_r^c+"
                + "3*(-1+l)*n_p*n^a*d_q^b*d_r^c+"
                + "6*((1/2)+l*b)*(n_p*n_q*g^ab*d_r^c+n^a*n^b*g_pq*d_r^c)+"
                + "6*(-(1/4)+l*b*b)*n_p*g_qr*n^a*g^bc");
        Tensor t = Transformations.expandBracketsExceptSymbols(Symmetrize.INSTANCE.transform(toInv));
        assertTrue(((Sum) t).size() == 60);
    }

    @Test
    public void test9() {
        Tensor u = CC.parse("A_ab*B_mn");
        CC.addSymmetry("A_mn", IndexType.LatinLower, false, new int[]{1, 0});
        CC.addSymmetry("B_mn", IndexType.LatinLower, false, new int[]{1, 0});
        Tensor sym = Symmetrize.INSTANCE.transform(u);
        Tensor expected = _("(1/6)*(A_{ab}*B_{mn}+A_{am}*B_{bn}+A_{an}*B_{bm}+A_{bm}*B_{an}+A_{bn}*B_{am}+A_{mn}*B_{ab})");
        assertEquals(sym, expected);
    }

    @Test
    public void test10() {
        Tensor u = CC.parse("A_ab*B_m*C_n");
        CC.addSymmetry("A_mn", IndexType.LatinLower, false, new int[]{1, 0});
        Tensor sym = Symmetrize.INSTANCE.transform(u);
        Tensor expected = _("(1/12)*(A_{ab}*B_{m}*C_{n}+A_{ab}*B_{n}*C_{m}+A_{am}*B_{b}*C_{n}+A_{am}*B_{n}*C_{b}+A_{an}*B_{b}*C_{m}+A_{an}*B_{m}*C_{b}+A_{bm}*B_{a}*C_{n}+A_{bm}*B_{n}*C_{a}+A_{bn}*B_{a}*C_{m}+A_{bn}*B_{m}*C_{a}+A_{mn}*B_{a}*C_{b}+A_{mn}*B_{b}*C_{a})");
        assertEquals(sym, expected);
    }

    @Test
    public void test11() {
        Tensor u = CC.parse("A_a*B_b*C_c*D_d");
        Tensor sym = Symmetrize.INSTANCE.transform(u);
        int size = -1;
        for (Tensor t : sym)
            if (t instanceof Sum) {
                size = ((Sum) t).size();
                break;
            }
        assertEquals(size, 24);
    }
}
