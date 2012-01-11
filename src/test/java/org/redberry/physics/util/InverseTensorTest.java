/*
 *
 * Redberry: symbolic tensor computations library.
 * Copyright (C) 2010-2011  Stanislav Poslavsky <stvlpos@mail.ru>
 *
 * This program is free software: you can redistribute it and/or
 * modify it under the terms of the GNU General Public License as
 * published by the Free Software Foundation, either version 3 of the
 * License, or (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 */
package org.redberry.physics.util;

import org.junit.Test;
import redberry.core.context.CC;
import redberry.core.tensor.Expression;
import redberry.core.tensor.Tensor;
import redberry.core.transformation.Transformation;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class InverseTensorTest {
    public InverseTensorTest() {
    }

    @Test
    public void test1() {
        SqrSubs sqrSubs = new SqrSubs(CC.parseSimple("k_a"));
        Transformation[] transformations = new Transformation[]{sqrSubs};
        Expression toInverse = new Expression("D_mn = k_m*k_n-(1/a)*k_i*k^i*g_mn");
        Expression equation = new Expression("D_ab*K^ac=d_b^c");
        Tensor[] samples = {CC.parse("g_mn"), CC.parse("g^mn"), CC.parse("d_m^n"), CC.parse("k_m"), CC.parse("k^b")};
        InverseTensor it = new InverseTensor(toInverse, equation, samples,transformations);
        for(Expression eq:it.linearEquations)
            System.out.println(eq);
    }
    
    @Test
    public void test2() {
        SqrSubs sqrSubs = new SqrSubs(CC.parseSimple("k_a"));
        Transformation[] transformations = new Transformation[]{sqrSubs};
        Expression toInverse = new Expression("D_mn = k_m*k_n-(1/a)*k_i*k^i*g_mn");
        Expression equation = new Expression("D_ab*K^ac=d_b^c");
        Tensor[] samples = {CC.parse("g_mn"), CC.parse("g^mn"), CC.parse("d_m^n"), CC.parse("k_m"), CC.parse("k^b")};
        InverseTensor it = new InverseTensor(toInverse, equation, samples,transformations);
        for(Expression eq:it.linearEquations)
            System.out.println(eq);
    }
}