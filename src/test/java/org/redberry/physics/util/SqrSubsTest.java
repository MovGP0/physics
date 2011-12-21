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

import org.redberry.physics.util.SqrSubs;
import org.junit.Test;
import redberry.core.context.CC;
import redberry.core.tensor.SimpleTensor;
import redberry.core.tensor.Tensor;
import redberry.core.transformation.Transformation;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class SqrSubsTest {
    public SqrSubsTest() {
    }

    @Test
    public void test1() {
        SimpleTensor n = (SimpleTensor) CC.parse("n_{a}");
        Transformation tr = new SqrSubs(n);
        Tensor t = CC.parse("n_m*n^m*a*n_a*n^a*n_i*n^j*b");
        System.out.println(tr.transform(t));
    }

    @Test
    public void test2() {
        SimpleTensor n = (SimpleTensor) CC.parse("n_{a}");
        Transformation tr = new SqrSubs(n);
        Tensor t = CC.parse("n_m*n^m*n_a*n^a");
        System.out.println(tr.transform(t));
    }
}
