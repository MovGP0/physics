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
package org.redberry.physics;

import org.junit.Test;
import redberry.core.context.CC;
import redberry.core.tensor.Tensor;
import redberry.core.transformation.Transformations;
import redberry.core.transformation.Transformer;
import redberry.core.transformation.collect.CollectPowers;
import static org.redberry.physics.TAssert.*;

/**
 *
 * @author Dmitry Bolotin
 * @author Stanislav Poslavsky
 */
public class ToFourierTest {
    public ToFourierTest() {
    }

    @Test
    public void test1() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[f[x],x]");
        assertParity(fourier.transform(target), "f[0]");
    }

    @Test
    public void test2() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[a*f[x]*b,x]");
        assertParity(fourier.transform(target), "a*b*f[0]");
    }

    @Test
    public void test3() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[a*f[x]*g[x]*b,x]");
        assertParity(fourier.transform(target), "Integral[a*f[y]*g[-y]*b,y]");
    }

    @Test
    public void test4() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[a*f[x]*g[x]*f[x]*f[x]*X[y1]*b,x]");
        assertParity(fourier.transform(target), "Integral[a*X[y1]*f[y]*g[y0]*f[y2]*f[-y-y0-y2]*b,y,y2,y0]");
    }

    @Test(expected = UnsupportedOperationException.class)
    public void test5() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[a*f[x,y]*g[x]*f[x]*f[x]*X[y1]*b,x]");
        System.out.println(fourier.transform(target));
    }

    /**
     * Derivatives
     */
    @Test
    public void testD1() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[D[f[x],x],x]");
        assertParity(fourier.transform(target), "0");
    }

    @Test
    public void testD2() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[f[x]*D[f[x],x],x]");
        assertParity(fourier.transform(target), "Integral[f[y]*I*(-1)*y*f[(-1)*y],y]");
    }

    @Test
    public void testD3() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[D[f[x],x]*f[x],x]");
        assertOpposite(fourier.transform(target), "Integral[f[y]*I*(-1)*y*f[(-1)*y],y]");
    }

    @Test
    public void testD23() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target1 = CC.parse("Integral[D[f[x],x]*f[x],x]");
        Tensor target2 = CC.parse("Integral[f[x]*D[f[x],x],x]");
        assertOpposite(fourier.transform(target1), fourier.transform(target2));
    }

    @Test
    public void testD4() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[D[f[x],x,x]*f[x],x]");
        assertParity(fourier.transform(target), "Integral[-f[y]*y*y*f[(-1)*y],y]");
    }

    @Test
    public void testD5() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[f[x]*D[f[x],x,x],x]");
        assertParity(Transformations.calculateNumbers(fourier.transform(target)), "Integral[-f[y]*y*y*f[(-1)*y],y]");
    }

    @Test
    public void testD6() {
        ToFourier fourier = new ToFourier("x", "y");
        Tensor target = CC.parse("Integral[D[g[x],x,x]*D[f[x],x,x],x]");
        target = new Transformer(CollectPowers.INSTANCE).transform(
                Transformations.calculateNumbers(fourier.transform(target)));
        assertParity(target, "Integral[g[y]*f[(-1)*y]*Pow[y,4],y]");
    }

    @Test
    public void testD7() {
        ToFourier fourier = new ToFourier("x_m", "p_m");
        Tensor target = CC.parse("Integral[D[g[x_m],x_a,x_b]*D[f[x_h],x_c,x_d],x^c]");
        target = new Transformer(CollectPowers.INSTANCE).transform(
                Transformations.calculateNumbers(fourier.transform(target)));
        assertParity(target, "Integral[p^a*p^b*p^c*p^d*g[p_y]*f[-p_x],p_x]");
    }
}
