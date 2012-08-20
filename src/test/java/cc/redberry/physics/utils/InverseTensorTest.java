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
package cc.redberry.physics.utils;

import cc.redberry.core.combinatorics.symmetries.Symmetries;
import cc.redberry.core.context.CC;
import cc.redberry.core.tensor.*;
import cc.redberry.core.tensor.Expression;
import cc.redberry.core.tensor.Tensor;
import cc.redberry.core.transformations.Symmetrize;
import cc.redberry.core.transformations.Transformation;
import org.junit.*;
/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class InverseTensorTest {


    @Test
    public void test1() {
        Transformation[] transformations = new Transformation[]{Tensors.parseExpression("k_a*k^a=1")};
        Expression toInverse = Tensors.parseExpression("D_mn = k_m*k_n-(1/a)*k_i*k^i*g_mn");
        Expression equation = Tensors.parseExpression("D_ab*K^ac=d_b^c");
        Tensor[] samples = {Tensors.parse("g_mn"), Tensors.parse("g^mn"), Tensors.parse("d_m^n"), Tensors.parse("k_m"), Tensors.parse("k^b")};
        InverseTensor it = new InverseTensor(toInverse, equation, samples, transformations);
        for (Expression eq : it.linearEquations)
            System.out.println(eq);
    }

    @Test
    public void test2() {
        Transformation[] transformations = new Transformation[]{Tensors.parseExpression("k_a*k^a=1")};
        Expression toInverse = Tensors.parseExpression("D_mn = k_m*k_n-(1/a)*k_i*k^i*g_mn");
        Expression equation = Tensors.parseExpression("D_ab*K^ac=d_b^c");
        Tensor[] samples = {Tensors.parse("g_mn"), Tensors.parse("g^mn"), Tensors.parse("d_m^n"), Tensors.parse("k_m"), Tensors.parse("k^b")};
        InverseTensor it = new InverseTensor(toInverse, equation, samples, transformations);
        for (Expression eq : it.linearEquations)
            System.out.println(eq);
    }

    @Test
    public void test3() {
        Transformation[] transformations = new Transformation[]{Tensors.parseExpression("k_a*k^a=1"), Tensors.parseExpression("d_a^a=4")};

        Tensor toInv = Tensors.parse("d_p^a*d_q^b*d_r^c+"
                + "6*(-(1/2)+l*b*b)*g_pq*g^ab*d_r^c+"
                + "3*(-1+l)*n_p*n^a*d_q^b*d_r^c+"
                + "6*((1/2)+l*b)*(n_p*n_q*g^ab*d_r^c+n^a*n^b*g_pq*d_r^c)+"
                + "6*(-(1/4)+l*b*b)*n_p*g_qr*n^a*g^bc");
        Expression toInverse = ExpressionFactory.FACTORY.create(Tensors.parseSimple("K^abc_pqr"),
                Symmetrize..transform(toInv));

        Tensor eqRhs = CC.parse("d_i^a*d_j^b*d_k^c");
        Expression equation = new Expression(CC.parse("K^abc_pqr*KINV^pqr_ijk"),
                Transformations.expandBracketsExceptSymbols(Symmetrize.INSTANCE.transform(eqRhs)));

        Symmetries symmetries = Symmetries.getFullSymmetriesForSortedIndices(3, 3);
        Tensor[] samples = {CC.parse("g_mn"), CC.parse("g^mn"), CC.parse("d_m^n"), CC.parse("n_m"), CC.parse("n^b")};
        InverseTensor it = new InverseTensor(toInverse, equation, symmetries, samples, transformations);
        System.out.println(it.inverse);

    }

    @Test
    public void test4() {
        SqrSubs sqrSubs = new SqrSubs(CC.parseSimple("n_a"));
        Transformation[] transformations = new Transformation[]{sqrSubs, new Expression("d_a^a=4").asSubstitution()};

        Tensor toInv = CC.parse("d_p^a*d_q^b*d_r^c+"
                + "6*(-(1/2)+l*b*b)*g_pq*g^ab*d_r^c+"
                + "3*(-1+l)*n_p*n^a*d_q^b*d_r^c+"
                + "6*((1/2)+l*b)*(n_p*n_q*g^ab*d_r^c+n^a*n^b*g_pq*d_r^c)+"
                + "6*(-(1/4)+l*b*b)*n_p*g_qr*n^a*g^bc");
        Expression toInverse = new Expression(CC.parseSimple("K^abc_pqr"),
                Symmetrize.INSTANCE.transform(toInv));

        Tensor eqRhs = CC.parse("d_i^a*d_j^b*d_k^c");
        Expression equation = new Expression(CC.parse("K^abc_pqr*KINV^pqr_ijk"),
                Transformations.expandBracketsExceptSymbols(Symmetrize.INSTANCE.transform(eqRhs)));

        Symmetries symmetries = Symmetries.getFullSymmetriesForSortedIndices(3, 3);
        Tensor[] samples = {CC.parse("g_mn"), CC.parse("g^mn"), CC.parse("d_m^n"), CC.parse("n_m"), CC.parse("n^b")};
        InverseTensor it = new InverseTensor(toInverse, equation, symmetries, samples, transformations);
    }

    @Test
    public void test5() {
        Transformation[] transformations = new Transformation[]{new Expression("d_a^a=4").asSubstitution()};

        Expression toInverse =
                new Expression("F_p^mn_q^rs = "
                + "d^s_q*d^r_p*g^mn+d^m_q*d^n_p*g^rs+(-1)*d^r_p*d^n_q*g^ms+(-1)*d^s_p*d^m_q*g^rn");
        Expression equation = new Expression("F_p^mn_q^rs*iF^p_mn^a_bc=d^a_q*d_b^r*d_c^s-(1/4)*d^r_q*d_b^a*d_c^s");

        Tensor[] samples = {CC.parse("g_mn"), CC.parse("g^mn"), CC.parse("d_m^n")};
        InverseTensor it = new InverseTensor(toInverse, equation, samples, transformations);
        System.out.println(it.inverse);
    }


}